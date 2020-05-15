% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-12-05 15:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S5PPRRR4_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR4_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PPRRR4_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:18:26
% EndTime: 2019-12-05 15:18:26
% DurationCPUTime: 0.24s
% Computational Cost: add. (315->63), mult. (825->94), div. (0->0), fcn. (1108->14), ass. (0->54)
t42 = sin(pkin(6));
t46 = cos(pkin(6));
t47 = cos(pkin(5));
t43 = sin(pkin(5));
t44 = cos(pkin(11));
t72 = t43 * t44;
t25 = -t42 * t72 + t47 * t46;
t40 = sin(pkin(11));
t45 = cos(pkin(10));
t41 = sin(pkin(10));
t74 = t41 * t47;
t28 = -t45 * t40 - t44 * t74;
t75 = t41 * t43;
t19 = -t28 * t42 + t46 * t75;
t77 = cos(qJ(3));
t73 = t43 * t40;
t71 = t45 * t43;
t70 = t45 * t47;
t68 = qJ(2) * t43;
t65 = qJ(1) + 0;
t64 = t42 * t77;
t63 = t46 * t77;
t62 = t45 * pkin(1) + t41 * t68 + 0;
t61 = t47 * qJ(2) + t65;
t60 = t43 * t64;
t26 = -t41 * t40 + t44 * t70;
t18 = -t26 * t42 - t46 * t71;
t59 = t41 * pkin(1) - t45 * t68 + 0;
t29 = -t40 * t74 + t45 * t44;
t58 = t29 * pkin(2) + t19 * pkin(7) + t62;
t57 = pkin(2) * t73 + t25 * pkin(7) + t61;
t50 = sin(qJ(3));
t10 = t29 * t77 + (t28 * t46 + t42 * t75) * t50;
t9 = -t28 * t63 + t29 * t50 - t41 * t60;
t56 = t10 * pkin(3) + t9 * pkin(8) + t58;
t16 = -t47 * t64 + t50 * t73 - t63 * t72;
t17 = t47 * t42 * t50 + (t44 * t46 * t50 + t77 * t40) * t43;
t55 = t17 * pkin(3) + t16 * pkin(8) + t57;
t27 = t40 * t70 + t41 * t44;
t54 = t27 * pkin(2) + t18 * pkin(7) + t59;
t7 = -t26 * t63 + t27 * t50 + t45 * t60;
t8 = t27 * t77 + (t26 * t46 - t42 * t71) * t50;
t53 = t8 * pkin(3) + t7 * pkin(8) + t54;
t52 = cos(qJ(4));
t51 = cos(qJ(5));
t49 = sin(qJ(4));
t48 = sin(qJ(5));
t12 = t17 * t52 + t25 * t49;
t11 = t17 * t49 - t25 * t52;
t4 = t10 * t52 + t19 * t49;
t3 = t10 * t49 - t19 * t52;
t2 = t18 * t49 + t8 * t52;
t1 = -t18 * t52 + t8 * t49;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t45, -t41, 0, 0; t41, t45, 0, 0; 0, 0, 1, t65; 0, 0, 0, 1; t29, t28, t75, t62; t27, t26, -t71, t59; t73, t72, t47, t61; 0, 0, 0, 1; t10, -t9, t19, t58; t8, -t7, t18, t54; t17, -t16, t25, t57; 0, 0, 0, 1; t4, -t3, t9, t56; t2, -t1, t7, t53; t12, -t11, t16, t55; 0, 0, 0, 1; t4 * t51 + t9 * t48, -t4 * t48 + t9 * t51, t3, t4 * pkin(4) + t3 * pkin(9) + t56; t2 * t51 + t7 * t48, -t2 * t48 + t7 * t51, t1, t2 * pkin(4) + t1 * pkin(9) + t53; t12 * t51 + t16 * t48, -t12 * t48 + t16 * t51, t11, t12 * pkin(4) + t11 * pkin(9) + t55; 0, 0, 0, 1;];
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
