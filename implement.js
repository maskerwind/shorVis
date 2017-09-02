/* made by Tao */

////////////////////// DEFINATION OF CLASSES /////////////////////
/// 2D-Point
class Point{
	constructor(x, y) {
        /**
         * The point's x coordinate.
         */
        this.x = x;
        /**
         * The point's y coordinate.
         */
        this.y = y;
	}

    times(c) {
        return new Point(this.x * c, this.y * c);
    }

    plus(p) {
        return new Point(this.x + p.x, this.y + p.y);
    }

    distanceTo(other) {
	    let dx = this.x - other.x;
	    let dy = this.y - other.y;
	    return Math.sqrt(dx*dx + dy*dy);
	}
}

/// State On Bloch Sphere
class State{
	constructor(theta, phi) {
        this.theta = theta;
        this.phi = phi;
	}
}

/// Axis On Bloch Sphere
class Axis{
	constructor(theta, phi) {
        this.theta = theta;
        this.phi = phi;
	}
}

/// Config Of Everything
class Config {}

////////////////////// DEFINATION OF CONSTANTS /////////////////////
/// Thickness
Config.DEFAULT_STROKE_THICKNESS = 1;
Config.DIM_STROKE_THICKNESS = 0.8;
Config.CIRCUIT_STROKE_THICKNESS = 1.5;
Config.TAG_THICKNESS = 3;

/// Color
Config.DEFAULT_STROKE_COLOR = '#FFFFFF';
Config.DEFAULT_FILL_COLOR = '#000000';
Config.CIRCLE_FILL_COLOR = '#202020';
Config.ARROW_FILL_COLOR = 'blue';
Config.ARROW_ENTANGLED_COLOR1 = '#00FF00';
Config.ARROW_ENTANGLED_COLOR2 = '#DC143C';
Config.ARROW_ENTANGLED_COLORA = '#FFFF00';
Config.TAG_COLORA = '#DC143C';
Config.TAG_COLORB = '#DC143C';
Config.GATE_COLORB = '#3B3B3B';

/// Single Bloch Canvas Size
Config.QUBIT_SIZE = 120;

/// Registers
Config.RGE_1_SIZE = 4;
Config.RGE_2_SIZE = 7;
Config.RGE_SINGLE = 1;
Config.MAX_N = 90;

Config.LINE_SIZE = 5;

/// Aminitation Loop
Config.ROTATE_RATE = 40;
Config.SUPER_RATE = 1000;


/// Circuit
Config.GAP_REG1 = 30;
Config.GAP_REG2 = 17;
Config.GAP_U = 40;
Config.CIRCUITW = 640;

Config.CIRCUIT_LOC_INIT = 35;
Config.CIRCUIT_LOC_H = 85;
Config.CIRCUIT_LOC_QFT = 520;
Config.CIRCUIT_LOC_MB = 430;
Config.CIRCUIT_LOC_MA = 590;

////// canvas ///////
var canvas_w = 0;
var canvas_h = 0;

///// factor  ///////
var N = 0;
var a = 2;

///// reg //////
var sizeReg1 = 0;
var sizeReg2 = 0;
var statesReg1 = new Array(); 
var statesReg2 = new Array(); 

///////////////////// DEFINATION OF VARIABLES //////////////////////
var canvasReg1=document.getElementById("Reg1");
var ctxReg1=Reg1.getContext("2d");
var canvasReg2=document.getElementById("Reg2");
var ctxReg2=Reg2.getContext("2d");
var canvasCircuit=document.getElementById("Circuit");
var ctxCircuit=Circuit.getContext("2d");

/// Radius Of Bloch Sphere
var radius = Config.QUBIT_SIZE* 37/100;

/// Rotation Loop
var timerRandom;
var timerSingle;
var timerCollapseA;
var rotate_counts = 0;
var initTheta = Math.PI/2;
var initPhi = 0;
var random_states = randomState(Config.RGE_1_SIZE);

/// Circuit 
var unitaryTimes = 0;
var temp_phase = Math.PI/12;
var fx = new Array();
var x = new Array();

// Measurement Loop
var superState_counts = 0;
var superState = new Array();
////////////////////// FUNCTIONS //////////////////////////////////
////////////// Generate ////////////////
function coordinateSystem(unit) {
    return {
        dx: new Point(unit / 3, -unit / 3),
        dy: new Point(unit, 0),
        dz: new Point(0, unit)
    };
}

function randomState(size = Config.RGE_SINGLE){
	var states = new Array();
	for(var i = 0; i < size; i++){
		thetaRom = Math.random()*2*Math.PI;
		phiRom = Math.random()*2*Math.PI;
		states[i] = new State(thetaRom, phiRom);
	}
	return states;
}

function initState(size, states){
	//var states = new Array();
	for(var i = 0; i < size; i++){
		thetaRom = 0;
		phiRom = 0;
		states[i] = new State(thetaRom, phiRom);
	}
	//return states;
}

////////////// Strock Registers////////////////
function strokeCircle(ctx, center, radius, color = Config.DEFAULT_STROKE_COLOR, thickness = Config.DEFAULT_STROKE_THICKNESS) {
	ctx.save();
    ctx.beginPath();
    ctx.arc(center.x, center.y, Math.max(radius - 0.5, 0), 0, 2 * Math.PI);
    ctx.closePath();
    ctx.strokeStyle = color;
    ctx.lineWidth = thickness;
    ctx.stroke();
    ctx.restore();
}

function strokeEllipse(ctx, center, radius, ratioX = 1 ,ratioY = 1 ,color = Config.DEFAULT_STROKE_COLOR, thickness = Config.DEFAULT_STROKE_THICKNESS){
   ctx.save();
   ctx.scale(ratioX, ratioY); //进行缩放（均匀压缩）
   ctx.moveTo((center.x + radius) / ratioX, center.y / ratioY);
   //从椭圆的左端点开始逆时针绘制
   ctx.beginPath();
   ctx.arc(center.x / ratioX, center.y / ratioY, radius, 0, 2 * Math.PI);
   ctx.closePath();
   ctx.strokeStyle = color;
   ctx.lineWidth = thickness
   ctx.stroke();
   ctx.restore();
}

function strokeLine(ctx, p1, p2, color = Config.DEFAULT_STROKE_COLOR, thickness = 1) {
	ctx.save();
    ctx.beginPath();
    ctx.moveTo(p1.x, p1.y);
    ctx.lineTo(p2.x, p2.y);
    ctx.strokeStyle = color;
    ctx.lineWidth = thickness;
    ctx.stroke();
    ctx.restore();
}

function strokeState(ctx, center, radius, theta, phi, color = Config.ARROW_FILL_COLOR){
	let u = radius;
	let {dx, dy, dz} = coordinateSystem(u);

	let px = -Math.sin(theta)*Math.cos(phi);
	let py = Math.sin(theta)*Math.sin(phi);
	let pz = -Math.cos(theta);

	let p = center.plus(dx.times(px)).plus(dy.times(py)).plus(dz.times(pz));
	let r = 5;
	// inidicator line
	strokeLine(ctx, center, p, Config.DEFAULT_STROKE_COLOR, 1.5);
	// inidicator circle head
	strokeCircle(ctx, p, r, Config.DEFAULT_STROKE_COLOR);
	fillCircle(ctx, p, r, color);
}

function strokeState_Cartesian(ctx, center, radius, px, py, pz, color = Config.ARROW_FILL_COLOR){
	let u = radius;
	let {dx, dy, dz} = coordinateSystem(u);
	let p = center.plus(dx.times(px)).plus(dy.times(py)).plus(dz.times(pz));
	//let p = center.plus(px).plus(py).plus(pz);
	let r = 5;
	// inidicator line
	strokeLine(ctx, center, p, Config.DEFAULT_STROKE_COLOR, 1.5);
	// inidicator circle head
	strokeCircle(ctx, p, r, Config.DEFAULT_STROKE_COLOR);
	fillCircle(ctx, p, r, color);
}

function strockBloch(ctx, center, radius){
	strokeCircle(ctx, center, radius);
	fillCircle(ctx, center, radius, Config.CIRCLE_FILL_COLOR);
	strokeEllipse(ctx, center, radius, 1 , 1/3 , Config.DEFAULT_STROKE_COLOR, Config.DIM_STROKE_THICKNESS);
	strokeEllipse(ctx, center, radius, 1/3 , 1, Config.DEFAULT_STROKE_COLOR , Config.DIM_STROKE_THICKNESS);
	strockAxis(ctx, center, radius, Config.DEFAULT_STROKE_COLOR, Config.DIM_STROKE_THICKNESS);
}

function strockReg(ctx, num, radius, canvas){
	var heightSize = Math.ceil( num / Config.LINE_SIZE);
	var widthSize = 0;
	(num >= Config.LINE_SIZE) ? widthSize = Config.LINE_SIZE: widthSize =  num;

	canvas_h = Config.QUBIT_SIZE * heightSize;
	canvas_w = Config.QUBIT_SIZE * widthSize;

	canvas.setAttribute("height",canvas_h);
	canvas.setAttribute("width",canvas_w);
	canvas.style.display = "";
	
	for (var i = 0; i < heightSize; i++){
		for (var j =0 ; j < num - i*Config.LINE_SIZE ; j++){
			var center = new Point(Config.QUBIT_SIZE/2 + j % (Config.LINE_SIZE + 1) * Config.QUBIT_SIZE , canvas_h/(heightSize*2)*(2*i+1) );
			strockBloch(ctx, center, radius);
		}
	}
}

function strockAxis(ctx, center, radius, color = Config.DEFAULT_STROKE_COLOR, thickness = Config.DEFAULT_STROKE_THICKNESS){
	let u = radius;
	let {dx, dy, dz} = coordinateSystem(u);

	var blochAxis = new Array();
	blochAxis[0] = new State(Math.PI/2,0);
	blochAxis[1] = new State(Math.PI/2,Math.PI/2);
	blochAxis[2] = new State(0,0);

	// stroke the dirac
	ctx.fillStyle = 'white';
	ctx.fillText('|0>',center.x-5,center.y- 1.1*radius); 
	ctx.fillText('|1>',center.x-5,center.y+ 1.3*radius); 

	for (var i = 0; i < 3; i++){
		let px = -Math.sin(blochAxis[i].theta)*Math.cos(blochAxis[i].phi);
		let py = Math.sin(blochAxis[i].theta)*Math.sin(blochAxis[i].phi);
		let pz = -Math.cos(blochAxis[i].theta);
		let p = center.plus(dx.times(px)).plus(dy.times(py)).plus(dz.times(pz));
		strokeLine(ctx, center, p, color, thickness);
	}
}

function fillCircle(ctx, center, radius, color = Config.DEFAULT_FILL_COLOR) {
	ctx.save();
    ctx.beginPath();
    ctx.arc(center.x, center.y, Math.max(radius - 0.5, 0), 0, 2 * Math.PI);
    ctx.closePath();
    ctx.fillStyle = color;
    ctx.fill();
    ctx.restore();
}

function clearCanvas(ctx, canvas, color = Config.DEFAULT_FILL_COLOR) {
    ctx.fillStyle = color;
    ctx.fillRect(0, 0, canvas.getAttribute("width"), canvas.getAttribute("height"));
}

function ClearQubit(ctx, statePos){
	var x = Config.QUBIT_SIZE * (statePos % Config.LINE_SIZE);
	var y = Config.QUBIT_SIZE * parseInt(statePos / Config.LINE_SIZE);
	var center = new Point(x + Config.QUBIT_SIZE/2, y + Config.QUBIT_SIZE/2);
	ctx.fillStyle = Config.DEFAULT_FILL_COLOR;
	ctx.fillRect(x, y, Config.QUBIT_SIZE, Config.QUBIT_SIZE);
	strockBloch(ctx, center, radius);
	return center;
}

function EntangleSingleState(ctx, statePos){
	var center = ClearQubit(ctx, statePos);
	strokeState(ctx, center, radius, 0, 0, Config.ARROW_ENTANGLED_COLOR1);
	strokeState(ctx, center, radius, Math.PI, 0, Config.ARROW_ENTANGLED_COLOR2);
}


////////////// Strock Circuits ////////////////////
function strockCircuit(ctx, canvas){
	var heightSize = Config.GAP_REG1 * (sizeReg1 + 2) + Config.GAP_REG2 * (sizeReg2 + 2);
	var gap = (sizeReg1+1) * Config.GAP_REG1;
	var sizeM1 = Config.GAP_REG1*5/6;
	var sizeM2 = Config.GAP_REG2*5/6;
	let pointbegin,pointend;

	canvas.setAttribute("height", heightSize);
	canvas.setAttribute("width", Config.CIRCUITW);
	canvas.style.display = "";

	// The First Reg
	for (var i = 0; i < sizeReg1; i++){
		pointbegin = new Point(0, (i+1) * Config.GAP_REG1);
		pointend = new Point(Config.CIRCUITW, (i+1) * Config.GAP_REG1);
		strokeCircuitLine(ctx, pointbegin, pointend);
		strokeCircuitHGate(ctx, pointbegin);

		pointbegin = new Point(550, (i+1) * Config.GAP_REG1 - sizeM1/2);
		strokeCircuitM(ctx, pointbegin, sizeM1);
	}
	
	// The Second Reg
	for (var i = 0; i < sizeReg2; i++){
		pointbegin = new Point(0, gap + (i+1) * Config.GAP_REG2);
		pointend = new Point(Config.CIRCUITW-150, gap + (i+1) * Config.GAP_REG2);
		strokeCircuitLine(ctx, pointbegin, pointend);
		pointbegin = new Point(400, gap + (i+1) * Config.GAP_REG2 - sizeM2/2);
		strokeCircuitM(ctx, pointbegin, sizeM2);
	}

	// The U gate
	for (var i = 0; i < sizeReg1; i++){
		pointbegin = new Point(60 + (i+1)*(300/sizeReg1), gap + Config.GAP_REG2);
		strokeCircuitU(ctx, pointbegin, i);

	}

	// The reverse QFT Gate
	pointbegin = new Point(400, Config.GAP_REG1/2);
	strokeCircuitRQFT(ctx, pointbegin);	
}

function strokeCircuitLine(ctx, beginPoint, endPoint, color = Config.DEFAULT_STROKE_COLOR, thickness = Config.CIRCUIT_STROKE_THICKNESS){
	ctx.font="18px Georgia white";
	ctx.fillStyle = 'white';
	ctx.fillText("|0〉",beginPoint.x,beginPoint.y+5);
	ctx.save();
	ctx.beginPath();
	ctx.moveTo(beginPoint.x+25, beginPoint.y);
	ctx.lineTo(endPoint.x, endPoint.y);
	ctx.closePath();

	ctx.strokeStyle = color;
	ctx.lineWidth = thickness;
	ctx.stroke();
    ctx.restore();
}

function strokeCircuitHGate(ctx, beginPoint, color = Config.DEFAULT_STROKE_COLOR, thickness = Config.CIRCUIT_STROKE_THICKNESS){
	var x = beginPoint.x + 50;
	var y = beginPoint.y - 15;
	var size = 25;
	ctx.clearRect(x, y, size,size);
	ctx.save();

	ctx.beginPath();

	ctx.rect(x, y, size, size);
	ctx.fillStyle = Config.GATE_COLORB;
	ctx.fill();

	ctx.font="13px Georgia";
	ctx.fillStyle = 'white';
	ctx.fillText("H", x+7, y+18);
	ctx.closePath();

	ctx.strokeStyle = color;
	ctx.lineWidth = thickness;
	ctx.stroke();
	ctx.restore();
}

function strokeCircuitU(ctx, beginPoint, count,  color = Config.DEFAULT_STROKE_COLOR, thickness = Config.CIRCUIT_STROKE_THICKNESS){

	var x = beginPoint.x;
	var y = beginPoint.y - 15;
	var width = 200 / sizeReg1;
	var height =  sizeReg2 * 23;
	var pointCenter = new Point(x + width/2, y - (Config.GAP_REG1 + 2 + Config.GAP_REG1*count));
	var pointRadius = 4;

	ctx.clearRect(x, y, width, height);
	ctx.save();
	ctx.beginPath();
	ctx.rect(x, y, width, height);
	ctx.fillStyle = Config.GATE_COLORB;
	ctx.fill();
	
	ctx.moveTo(pointCenter.x, y);
	ctx.lineTo(pointCenter.x, pointCenter.y);

	ctx.font= 20 - sizeReg1 + "px Georgia";
	ctx.fillStyle = 'white';
	ctx.fillText("U^" + count, x + 2, y + 45);
	ctx.closePath();
	ctx.strokeStyle = color;
	ctx.lineWidth = thickness;
	ctx.stroke();
	ctx.restore();
	strokeCircle(ctx, pointCenter, pointRadius);
	fillCircle(ctx, pointCenter, pointRadius, '#FFFFFF');
}

function strokeCircuitRQFT(ctx, beginPoint, color = Config.DEFAULT_STROKE_COLOR, thickness = Config.CIRCUIT_STROKE_THICKNESS){
	var x = beginPoint.x;
	var y = beginPoint.y;
	var width = 100;
	var height =  sizeReg1 * Config.GAP_REG1;

	ctx.clearRect(x, y, width, height);
	ctx.save();
	ctx.beginPath();
	ctx.rect(x, y, width, height);
	ctx.fillStyle = Config.GATE_COLORB;
	ctx.fill();

	ctx.font= "20px Georgia";
	ctx.fillStyle = 'white';
	ctx.fillText("Reverse ", x + 15, height/2);
	ctx.fillText("QFT", x + 25, height/2+30);
	ctx.closePath();
	ctx.strokeStyle = color;
	ctx.lineWidth = thickness;
	ctx.stroke();
	ctx.restore();
}

function strokeCircuitM(ctx, beginPoint, size, color = Config.DEFAULT_STROKE_COLOR, thickness = Config.CIRCUIT_STROKE_THICKNESS) {
	var x = beginPoint.x;
	var y = beginPoint.y;
	ctx.clearRect(x, y, size, size);
	ctx.save();
	
	ctx.beginPath();
	ctx.rect(x, y, size, size);
	
	ctx.fillStyle = Config.GATE_COLORB;
	ctx.fill();
	
	ctx.strokeStyle = color;
	ctx.lineWidth = thickness;
	ctx.stroke();

	ctx.beginPath();
	ctx.arc(x+size/2, y+size, size*3/5, 125/100*Math.PI,175/100*Math.PI);
	ctx.lineWidth = 1;
	ctx.stroke();

	ctx.beginPath();
    ctx.moveTo(x+ size/2, y + 4/5*size);
    ctx.lineTo(x+ size*5/6, y + 1/5*size);
    ctx.lineWidth = 1;
	ctx.stroke();
	ctx.restore();
}


////////////// Math Functions /////////////////////
function dotProduct(v1, v2) {
    return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
}

function Matrix_Multiply(matrix, v) {
    return [dotProduct(matrix[0], v), dotProduct(matrix[1], v), dotProduct(matrix[2], v)];
}

function State_Rotation(state , axis , beta) {
	let p = [Math.sin(state.theta)*Math.cos(state.phi), Math.sin(state.theta)*Math.sin(state.phi), Math.cos(state.theta)],
		v = [Math.sin(axis.theta)*Math.cos(axis.phi), Math.sin(axis.theta)*Math.sin(axis.phi), Math.cos(axis.theta)];
    let ca = Math.cos(beta), sa = Math.sin(beta), t=1-ca, x=v[0], y=v[1], z=v[2];
    let r = [
        [ca + x*x*t, x*y*t - z*sa, x*z*t + y*sa],
        [x*y*t + z*sa, ca + y*y*t, y*z*t - x*sa],
        [z*x*t - y*sa, z*y*t + x*sa, ca + z*z*t]
    ];
    return Matrix_Multiply(r, p);
}

function Check_N(){
	var N=document.getElementById("num_N").value;
	if(isNaN(Number(N)) || N % 2 == 0 || N < 0){
		alert('Please input a suitable number');
		return false;
	}
	else if(N > 90){
		alert('N too large');
		return false;
	}
	return N;
}

function GCD(a,b){
	let d;
	if(b != 0){
		while (a % b != 0){
			d = a % b;
			a = b;
			b = d;
		}
	}
	return b;
}


////////////// Auxiliary Functions ////////////////
function Refresh_Registers(){
	clearCanvas(ctxReg1, canvasReg1);
	clearCanvas(ctxReg2, canvasReg2);
	strockReg(ctxReg1, sizeReg1, radius, canvasReg1);
	strockReg(ctxReg2, sizeReg2, radius, canvasReg2);
}

function Refresh_Tag(locA, locB){
	strockCircuit(ctxCircuit,canvasCircuit);

	var beginPoint = new Point(locA, Config.GAP_REG1 - 10);
	var endPoint = new Point(locA, sizeReg1*Config.GAP_REG1 + 10 );
	strokeLine(ctxCircuit, beginPoint, endPoint, Config.TAG_COLORA, Config.TAG_THICKNESS);

	beginPoint = new Point(locB, (sizeReg1+1)*Config.GAP_REG1 + Config.GAP_REG2 - 10);
	endPoint = new Point(locB, (sizeReg1+1)*Config.GAP_REG1 + Config.GAP_REG2*sizeReg2 + 10 );
	strokeLine(ctxCircuit, beginPoint, endPoint, Config.TAG_COLORA, Config.TAG_THICKNESS);
}

function PhaseShift_Block(ctx, statesReg, statePos, phase){
	var P_axis = new Axis(0, 0);
	var H_final_state_matrix = State_Rotation(statesReg[statePos], P_axis, phase);
	var center = ClearQubit(ctx, statePos);
	strokeState_Cartesian(ctx, center, radius, -H_final_state_matrix[0], H_final_state_matrix[1], -H_final_state_matrix[2],Config.ARROW_ENTANGLED_COLORA);
}

function CollapseDraw(ctx, sizeReg, collapseState){
	var statePos = sizeReg - 1;
	for(var i = collapseState.length - 1; i >= 0; i--){
		var center = ClearQubit(ctx, statePos);
		if(collapseState.charAt(i) == '0'){
			strokeState(ctx, center, radius, 0, 0, Config.ARROW_ENTANGLED_COLOR1);
		}
		else{
			strokeState(ctx, center, radius, Math.PI, 0, Config.ARROW_ENTANGLED_COLOR1);
		}
		statePos--;
	}
	while(statePos >= 0){
		var center = ClearQubit(ctx, statePos);
		strokeState(ctx, center, radius, 0, 0, Config.ARROW_ENTANGLED_COLOR1);
		statePos--;
	}
}

function CollapseA(ctx, sizeReg, superState){
	if(superState_counts == superState.length){
		superState_counts = 0;
		return;
	}
	CollapseDraw(ctx, sizeReg, superState[superState_counts].toString(2));
	var value_current = document.getElementById("value_current");
	value_current.style.display = "";
	value_current.innerHTML ="  [CURRENT]:  " + superState[superState_counts] + "〉";

	superState_counts ++;
	timerCollapseA = setTimeout("CollapseA(ctxReg1, sizeReg1, superState)",Config.SUPER_RATE);
}


////////////// Animation Functions////////////////
function Show_Registers(){
	if(!Check_N()){
		return;
	}
	N = Check_N();
	 // sizeReg2 = l; sizeReg1 = 2*l
	sizeReg2 = Math.ceil(Math.log2(Number(N)));
	sizeReg1 = 2*sizeReg2;
	unitaryTimes = 0;

	var size_Reg1 = document.getElementById("size_Reg1");
	var size_Reg2 = document.getElementById("size_Reg2");
	var value_Reg1 = document.getElementById("value_Reg1");
	var value_Reg2 = document.getElementById("value_Reg2");

	var description = document.getElementById("des");

	size_Reg1.innerHTML = "Size: "+ sizeReg1 +" qubits";
	size_Reg2.innerHTML = "Size:  "+ sizeReg2 +" qubits";
	value_Reg1.innerHTML = "Value: NaN";
	value_Reg2.innerHTML = "Value: NaN";

	description.innerHTML = "<p style = 'margin-left: 10px; color: white' >The first step is to initialize registers 1 and 2. The size of the second register {t} are supposed to satisfy the function:</p><p style = 'margin-left: 50px'>N<sup>2</sup>≤t≤2N<sup>2</sup></p>";
	strockCircuit(ctxCircuit,canvasCircuit);
	Refresh_Registers();
}

function Animation_Initiate(){
	var heightSize1 = Math.ceil( sizeReg1/ Config.LINE_SIZE);
	var heightSize2 = Math.ceil( sizeReg2/ Config.LINE_SIZE);

	var value_Reg1 = document.getElementById("value_Reg1");
	var value_Reg2 = document.getElementById("value_Reg2");
	var description = document.getElementById("des");

	initState(sizeReg1,statesReg1);
	initState(sizeReg2 - 1,statesReg2);
	statesReg2[sizeReg2-1] = new State(Math.PI, 0);

	value_Reg1.innerHTML = "Value: |0〉<sup>⊗" + sizeReg1 + "</sup>";
	value_Reg2.innerHTML = "Value: |1〉";
	description.innerHTML = "<p style = 'margin-left: 10px; color: white' >Next we need to set the value of each qubit to |0〉</p>";
	// stroke the bloch sphere first
	Refresh_Registers();

	// stroke the states
	for (var i = 0; i < heightSize1; i++){
		for (var j =0 ; j < sizeReg1 - i*Config.LINE_SIZE ; j++){
			var center = new Point(Config.QUBIT_SIZE/2 + j % (Config.LINE_SIZE + 1) * Config.QUBIT_SIZE , Config.QUBIT_SIZE * heightSize1/(heightSize1*2)*(2*i+1) );
			strokeState(ctxReg1, center, radius, statesReg1[i*Config.LINE_SIZE + j].theta , statesReg1[i*Config.LINE_SIZE + j].phi);
		}
	}
	for (var i = 0; i < heightSize2; i++){
		for (var j =0 ; j < sizeReg2 - i*Config.LINE_SIZE ; j++){
			var center = new Point(Config.QUBIT_SIZE/2 + j % (Config.LINE_SIZE + 1) * Config.QUBIT_SIZE , Config.QUBIT_SIZE * heightSize2/(heightSize2*2)*(2*i+1) );
			strokeState(ctxReg2, center, radius, statesReg2[i*Config.LINE_SIZE + j].theta , statesReg2[i*Config.LINE_SIZE + j].phi);
		}
	}
	// stroke the tag on circuit
	Refresh_Tag(Config.CIRCUIT_LOC_INIT, Config.CIRCUIT_LOC_INIT);
}

function Animation_Hardmard(){
	// exit
	if (rotate_counts==Config.ROTATE_RATE + 1) {
		var value_Reg1 = document.getElementById("value_Reg1");
		var description = document.getElementById("des");
		value_Reg1.innerHTML = "Value: <sup>1</sup>&frasl;<sub>2<sup>"+ sizeReg1 + "</sup></sub>( 0〉+1〉+...+2<sup>" + sizeReg1 + "-1</sup>〉 )";
		description.innerHTML = "<p style = 'margin-left: 10px; color: white' >One application of the Hadamard gate to either a 0 or 1 qubit will produce a quantum state that, if observed, will be a 0 or 1 with equal probability (as seen in the first two operations).</p><p style = 'margin-left: 10px;color: white' > The law of hadamard transformation follows the following formular:</p> <p><img style='filter: invert(100%)' src='image/hadamard.svg'/> </p><p style = 'margin-left: 10px; color: white' >This would be like taking a fair coin that is showing heads, flipping it twice, and it always landing on heads after the second flip.</p> <p style = 'margin-left: 10px; color: white' >The Hadamard operation is a 180 degree rotation around the diagonal X+Z axis of the Bloch sphere.</p>";
		rotate_counts = 0;
		
		// Change the final states
		for(var i = 0; i < sizeReg1; i++){
			statesReg1[i].theta = Math.PI / 2;
		}

		// Refresh tag on circuit
		Refresh_Tag(Config.CIRCUIT_LOC_H, Config.CIRCUIT_LOC_INIT);
		return;
	}
	var heightSize1 = Math.ceil( sizeReg1/ Config.LINE_SIZE);
	var H_axis = new Axis(Math.PI/4, 0);
	var beta = Math.PI;
	var temp_beta = beta / Config.ROTATE_RATE * rotate_counts;

	
	// Refresh the first register
	clearCanvas(ctxReg1, canvasReg1);
	strockReg(ctxReg1, sizeReg1, radius, canvasReg1);

	for (var i = 0; i < heightSize1; i++){
		for (var j =0 ; j < sizeReg1 - i*Config.LINE_SIZE ; j++){
			var center = new Point(Config.QUBIT_SIZE/2 + j % (Config.LINE_SIZE + 1) * Config.QUBIT_SIZE , Config.QUBIT_SIZE * heightSize1/(heightSize1*2)*(2*i+1) );
			var H_final_state_matrix = State_Rotation(statesReg1[i*Config.LINE_SIZE + j], H_axis, temp_beta);
			strokeState_Cartesian(ctxReg1, center, radius, -H_final_state_matrix[0], H_final_state_matrix[1], -H_final_state_matrix[2]);
		}
	}

	rotate_counts += 1;
	timerSingle = setTimeout("Animation_Hardmard()",Config.ROTATE_RATE);
}

function Animation_Unitary(){
	var value_Reg2 = document.getElementById("value_Reg2");
	value_Reg2.innerHTML = "Value: INVISIBLE";

	// Unitary Change On Bloch Sphere
	unitaryTimes += 1;
	//// Entangle Each Qubit in first reg with all Qubits in second reg
	var statePos1 = sizeReg1 - unitaryTimes;
	var statePos2 = sizeReg2 - 1;
	//Unitary_Block
	PhaseShift_Block(ctxReg1, statesReg1 ,statePos1, temp_phase * Math.pow(2, unitaryTimes - 1));
	for (var statePos2 = 0; statePos2 < sizeReg2; statePos2++){
		EntangleSingleState(ctxReg2, statePos2);
	}

	// Refresh Tag on Circuit
	var loc = 60 + unitaryTimes*(300/sizeReg1) + 200/sizeReg1 + 50/sizeReg1 ;
	Refresh_Tag(loc, loc);

	//// number
	if( unitaryTimes == sizeReg1){
		for(var i = 0; i < Math.pow(2,sizeReg2); i++){
			fx[i] = Math.pow(a,i) % N;
		}
		for(var j = 0; j < Math.pow(2,sizeReg1); j++){
			x[j]  = fx[ j % Math.pow(2,sizeReg2) ];
		}
		unitaryTimes = 0;
		value_Reg2.innerHTML = "Value: <sup>1</sup>&frasl;<sub>2<sup>"+ sizeReg2 + "</sup></sub>( "+ fx[0] + "〉+|" + fx[1] + "〉+|" + fx[2] + "〉+...+|" + fx[fx.length-1] + "〉 )";
		description.innerHTML = "<p style = 'margin-left: 10px; color: white' ></p>"
		var description = document.getElementById("des");
		return;
	}
}

function Animation_MeasureB(){
	var collapseState = fx[parseInt(Math.random()* Math.pow(2, sizeReg2))];
	var i = 0;
	for(var j = 0; j < x.length ; j++){
		if(x[j] == collapseState){
			superState[i] = j;
			i++;
		}
	}

	var value_Reg1 = document.getElementById("value_Reg1");
	var value_Reg2 = document.getElementById("value_Reg2");

	value_Reg1.innerHTML = "Value: <sup>1</sup>&frasl;<sub>2<sup>"+ sizeReg1 + "-1</sup></sub>( "+ superState[0] + "〉+" + superState[1] + "〉+...+" + superState[superState.length-1] + "〉 )";
	value_Reg2.innerHTML = "Value: " + collapseState + "〉";

	CollapseDraw(ctxReg2, sizeReg2 , collapseState.toString(2));
	CollapseA(ctxReg1, sizeReg1, superState);

	var loc = 360 + 250/sizeReg1;
	Refresh_Tag(loc, Config.CIRCUIT_LOC_MB);
}

function Animation_InverseQFT(){
	Refresh_Tag(Config.CIRCUIT_LOC_QFT, Config.CIRCUIT_LOC_MB);
}

function Animation_MeasureA(){
	var l = superState.length;
	var v = Math.pow(2,sizeReg1);
	let quotient, temp;
	let r, p1, q1, p, q; 
	var result = new Array();
	var count = 0;

	Refresh_Tag(Config.CIRCUIT_LOC_MA, Config.CIRCUIT_LOC_MB);

	while (l>1){
		result[count] = Math.floor(v/l);
		count ++;
		temp = l;
		l = v % l;
		v = temp;
	}
	result[count] = v;

	r = result[0];
	if(r % 2 != 0){
		alert("failed!");
		return;
	}
	else{
		p1= Math.pow(a, r/2)+1;
		q1= Math.pow(a, r/2)-1;
		p = Math.max(GCD(p1,N), GCD(q1,N));
		q = N / p;
	}

	var ele_p = document.getElementById('result_p');
	var ele_q = document.getElementById('result_q');
	ele_p.innerHTML ="p: " + p;
	ele_p.style.display = "";
	
	ele_q.innerHTML ="q: " + q;
	ele_q.style.display = "";
}

