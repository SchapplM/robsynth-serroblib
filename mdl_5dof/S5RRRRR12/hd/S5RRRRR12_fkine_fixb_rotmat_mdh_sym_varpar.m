% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d1,d2,d3,d4,d5]';
% 
% Output:
% T_c_mdh [4x4x(5+1)]
%   homogenous transformation matrices for each (body) frame (MDH)
%   1:  mdh base (link 0) -> mdh base link 0 (unit matrix, no information)
%   ...
%   6:  mdh base (link 0) -> mdh frame (6-1), link (6-1)
%   ...
%   5+1:  mdh base (link 0) -> mdh frame (5)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S5RRRRR12_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR12_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RRRRR12_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:46:24
% EndTime: 2019-12-31 22:46:25
% DurationCPUTime: 0.22s
% Computational Cost: add. (315->63), mult. (825->91), div. (0->0), fcn. (1108->14), ass. (0->55)
t43 = cos(pkin(5));
t47 = sin(qJ(2));
t52 = cos(qJ(1));
t69 = t52 * t47;
t48 = sin(qJ(1));
t51 = cos(qJ(2));
t71 = t48 * t51;
t28 = -t43 * t71 - t69;
t40 = sin(pkin(6));
t42 = cos(pkin(6));
t41 = sin(pkin(5));
t73 = t48 * t41;
t19 = -t28 * t40 + t42 * t73;
t75 = t41 * t51;
t25 = -t40 * t75 + t43 * t42;
t78 = cos(qJ(3));
t76 = t41 * t47;
t72 = t48 * t47;
t70 = t52 * t41;
t68 = t52 * t51;
t67 = pkin(7) + 0;
t64 = t40 * t78;
t63 = t42 * t78;
t62 = t43 * pkin(8) + t67;
t61 = t52 * pkin(1) + pkin(8) * t73 + 0;
t60 = t41 * t64;
t26 = t43 * t68 - t72;
t18 = -t26 * t40 - t42 * t70;
t59 = t48 * pkin(1) - pkin(8) * t70 + 0;
t29 = -t43 * t72 + t68;
t58 = t29 * pkin(2) + t19 * pkin(9) + t61;
t57 = pkin(2) * t76 + t25 * pkin(9) + t62;
t46 = sin(qJ(3));
t11 = -t28 * t63 + t29 * t46 - t48 * t60;
t12 = t29 * t78 + (t28 * t42 + t40 * t73) * t46;
t56 = t12 * pkin(3) + t11 * pkin(10) + t58;
t16 = -t43 * t64 + t46 * t76 - t63 * t75;
t17 = t43 * t40 * t46 + (t42 * t46 * t51 + t78 * t47) * t41;
t55 = t17 * pkin(3) + t16 * pkin(10) + t57;
t27 = t43 * t69 + t71;
t54 = t27 * pkin(2) + t18 * pkin(9) + t59;
t10 = t27 * t78 + (t26 * t42 - t40 * t70) * t46;
t9 = -t26 * t63 + t27 * t46 + t52 * t60;
t53 = t10 * pkin(3) + t9 * pkin(10) + t54;
t50 = cos(qJ(4));
t49 = cos(qJ(5));
t45 = sin(qJ(4));
t44 = sin(qJ(5));
t8 = t17 * t50 + t25 * t45;
t7 = t17 * t45 - t25 * t50;
t4 = t12 * t50 + t19 * t45;
t3 = t12 * t45 - t19 * t50;
t2 = t10 * t50 + t18 * t45;
t1 = t10 * t45 - t18 * t50;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t52, -t48, 0, 0; t48, t52, 0, 0; 0, 0, 1, t67; 0, 0, 0, 1; t29, t28, t73, t61; t27, t26, -t70, t59; t76, t75, t43, t62; 0, 0, 0, 1; t12, -t11, t19, t58; t10, -t9, t18, t54; t17, -t16, t25, t57; 0, 0, 0, 1; t4, -t3, t11, t56; t2, -t1, t9, t53; t8, -t7, t16, t55; 0, 0, 0, 1; t11 * t44 + t4 * t49, t11 * t49 - t4 * t44, t3, t4 * pkin(4) + t3 * pkin(11) + t56; t2 * t49 + t9 * t44, -t2 * t44 + t9 * t49, t1, t2 * pkin(4) + t1 * pkin(11) + t53; t16 * t44 + t8 * t49, t16 * t49 - t8 * t44, t7, t8 * pkin(4) + t7 * pkin(11) + t55; 0, 0, 0, 1;];
T_ges = t5;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
Tc_stack = NaN(3*size(T_c_mdh,3),4);
% Zusätzliche Ausgabe: Als 2D-array gestapelt, ohne Zeile mit 0001
for i = 1:size(T_c_mdh,3), Tc_stack((i-1)*3+1:3*i,1:4) = T_c_mdh(1:3,1:4,i); end
