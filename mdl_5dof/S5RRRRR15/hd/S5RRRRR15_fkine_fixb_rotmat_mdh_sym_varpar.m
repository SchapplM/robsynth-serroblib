% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRRRR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d1,d2,d3,d4,d5]';
% 
% Output:
% Tc_mdh [4x4x(5+1)]
%   homogenous transformation matrices for each (body) frame (MDH)
%   1:  mdh base (link 0) -> mdh base link 0 (unit matrix, no information)
%   ...
%   6:  mdh base (link 0) -> mdh frame (6-1), link (6-1)
%   ...
%   5+1:  mdh base (link 0) -> mdh frame (5)
% T_c_stack [(5+1)*3 x 4]
%   stacked matrices from Tc_mdh into one 2D array, last row left out.
%   Last row only contains [0 0 0 1].

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 22:28
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [Tc_mdh, Tc_stack] = S5RRRRR15_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR15_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RRRRR15_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 22:23:44
% EndTime: 2024-09-27 22:23:44
% DurationCPUTime: 0.09s
% Computational Cost: add. (277->98), mult. (296->139), div. (0->0), fcn. (365->26), ass. (0->79)
t46 = sin(qJ(1));
t80 = -t46 / 0.2e1;
t51 = cos(qJ(1));
t79 = t51 / 0.2e1;
t52 = pkin(8) + pkin(9);
t38 = sin(pkin(6));
t78 = pkin(11) * t38;
t50 = cos(qJ(2));
t29 = t50 * pkin(2) + pkin(1);
t37 = qJ(2) + qJ(3);
t34 = qJ(4) + t37;
t27 = cos(t34);
t40 = cos(pkin(6));
t77 = t27 * t40;
t41 = cos(pkin(5));
t76 = t38 * t41;
t39 = sin(pkin(5));
t45 = sin(qJ(2));
t75 = t39 * t45;
t26 = sin(t34);
t74 = t46 * t26;
t73 = t46 * t27;
t19 = t46 * t39;
t42 = sin(qJ(5));
t72 = t46 * t42;
t71 = t46 * t45;
t70 = t46 * t50;
t47 = cos(qJ(5));
t69 = t47 * t46;
t68 = t51 * t26;
t67 = t51 * t27;
t66 = t51 * t39;
t65 = t51 * t42;
t64 = t51 * t45;
t63 = t51 * t47;
t62 = t51 * t50;
t36 = pkin(10) + t52;
t61 = pkin(7) + 0;
t60 = t41 * t65;
t59 = t38 * t19;
t58 = t41 * t72;
t57 = t41 * t69;
t56 = t38 * t66;
t55 = t41 * t63;
t31 = pkin(5) - t37;
t30 = pkin(5) + t37;
t43 = sin(qJ(4));
t48 = cos(qJ(4));
t12 = -t43 * pkin(4) + t48 * t78;
t44 = sin(qJ(3));
t49 = cos(qJ(3));
t9 = pkin(4) * t48 + t43 * t78 + pkin(3);
t3 = t12 * t44 + t9 * t49 + pkin(2);
t4 = t12 * t49 - t44 * t9;
t54 = t3 * t45 - t4 * t50;
t53 = pkin(3) * t44 * t50 + (t49 * pkin(3) + pkin(2)) * t45;
t33 = cos(t37);
t32 = sin(t37);
t25 = -qJ(4) + t31;
t24 = qJ(4) + t30;
t23 = cos(t31);
t22 = cos(t30);
t21 = sin(t31);
t20 = sin(t30);
t18 = cos(t24);
t17 = sin(t25);
t16 = t40 * pkin(11) + t36;
t15 = cos(t25) / 0.2e1;
t14 = sin(t24) / 0.2e1;
t13 = pkin(3) * t33 + t29;
t11 = t23 + t22;
t10 = t20 - t21;
t8 = t41 * t45 * pkin(2) - t39 * t52;
t7 = t15 + t18 / 0.2e1;
t6 = t14 - t17 / 0.2e1;
t5 = -t39 * t36 + t53 * t41;
t2 = t3 * t50 + t4 * t45 + pkin(1);
t1 = t39 * t16 - t54 * t41;
t28 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; t51, -t46, 0, 0; t46, t51, 0, 0; 0, 0, 1, t61; -t41 * t71 + t62, -t41 * t70 - t64, t19, t51 * pkin(1) + pkin(8) * t19 + 0; t41 * t64 + t70, t41 * t62 - t71, -t66, t46 * pkin(1) - pkin(8) * t66 + 0; t75, t39 * t50, t41, t41 * pkin(8) + t61; t10 * t80 + t51 * t33, t11 * t80 - t51 * t32, t19, t51 * t29 - t8 * t46 + 0; t10 * t79 + t46 * t33, t11 * t79 - t46 * t32, -t66, t46 * t29 + t51 * t8 + 0; t23 / 0.2e1 - t22 / 0.2e1, t20 / 0.2e1 + t21 / 0.2e1, t41, pkin(2) * t75 + t41 * t52 + t61; -t46 * t6 + t67, -t46 * t7 - t68, t19, t51 * t13 - t46 * t5 + 0; t51 * t6 + t73, t51 * t7 - t74, -t66, t46 * t13 + t51 * t5 + 0; t15 - t18 / 0.2e1, t14 + t17 / 0.2e1, t41, t41 * t36 + t53 * t39 + t61; (-t40 * t58 + t63) * t27 + (-t40 * t65 - t57) * t26 + t42 * t59, (-t40 * t57 - t65) * t27 + (-t40 * t63 + t58) * t26 + t47 * t59, t40 * t19 + (t41 * t73 + t68) * t38, t1 * t46 + t2 * t51 + 0; (t40 * t60 + t69) * t27 + (-t40 * t72 + t55) * t26 - t42 * t56, (t40 * t55 - t72) * t27 + (-t40 * t69 - t60) * t26 - t47 * t56, -t40 * t66 + (-t41 * t67 + t74) * t38, -t1 * t51 + t2 * t46 + 0; t42 * t76 + (t26 * t47 + t42 * t77) * t39, t47 * t76 + (-t26 * t42 + t47 * t77) * t39, -t39 * t27 * t38 + t41 * t40, t16 * t41 + t54 * t39 + t61;];
Tc_stack = t28;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), Tc_mdh = NaN(4,4,5+1);               % numerisch
else,                         Tc_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  Tc_mdh(:,:,i) = [Tc_stack((i-1)*3+1 : 3*i, :);[0 0 0 1]];
end
