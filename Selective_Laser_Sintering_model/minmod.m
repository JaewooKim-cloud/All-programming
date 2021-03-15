function answer=minmod(a,b)
    answer=0.5*(sign(a)+sign(b))*min(abs(a),abs(b));
end