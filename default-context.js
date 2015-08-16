var DefaultBezierContext = exports.DefaultBezierContext = function(){
	this.strands = [];
	this.lastx = 0;
	this.lasty = 0;
}
DefaultBezierContext.prototype.moveTo = function(x, y) {
	this.lastx = x;
	this.lasty = y;
}
DefaultBezierContext.prototype.lineTo = function(x, y) {
	this.strands.push({
		order: 1,
		start: {x: this.lastx, y: this.lasty}, 
		end: {x: x, y: y}
	});
	this.lastx = x;
	this.lasty = y;
}
DefaultBezierContext.prototype.cubicTo = function(x1, y1, x2, y2, x, y) {
	this.strands.push({
		order: 3,
		start: {x: this.lastx, y: this.lasty},
		c1: {x: x1, y: y1},
		c2: {x: x2, y: y2},
		end: {x: x, y: y}
	});
	this.lastx = x;
	this.lasty = y;
}