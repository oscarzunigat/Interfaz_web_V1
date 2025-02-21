from flask import Blueprint, render_template, request

pagina1_bp = Blueprint('pagina1', __name__)

@pagina1_bp.route('/pagina1')
def pagina1():
    return render_template('pagina1.html')

@pagina1_bp.route('/suma1', methods=['POST'])
def suma1():
    num1 = request.form.get('num1', type=int)
    num2 = request.form.get('num2', type=int)
    resultado = num1 + num2
    return f'La suma de {num1} y {num2} es {resultado}.'
