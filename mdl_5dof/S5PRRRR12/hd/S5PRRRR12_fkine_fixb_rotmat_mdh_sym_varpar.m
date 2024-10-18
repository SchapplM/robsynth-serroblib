% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PRRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d2,d3,d4,d5,theta1]';
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
% Datum: 2024-09-28 18:09
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [Tc_mdh, Tc_stack] = S5PRRRR12_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR12_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR12_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-28 18:06:58
% EndTime: 2024-09-28 18:06:58
% DurationCPUTime: 0.03s
% Computational Cost: add. (277->117), mult. (334->178), div. (0->0), fcn. (399->26), ass. (0->79)
t57 = pkin(7) + pkin(8);
t44 = sin(pkin(6));
t83 = pkin(10) * t44;
t56 = cos(qJ(2));
t32 = t56 * pkin(2) + pkin(1);
t42 = qJ(2) + qJ(3);
t39 = qJ(4) + t42;
t29 = sin(t39);
t43 = sin(pkin(11));
t82 = t43 * t29;
t45 = sin(pkin(5));
t24 = t43 * t45;
t48 = cos(pkin(5));
t81 = t43 * t48;
t80 = t44 * t45;
t52 = sin(qJ(2));
t79 = t45 * t52;
t46 = cos(pkin(11));
t78 = t46 * t29;
t77 = t46 * t45;
t76 = t46 * t48;
t47 = cos(pkin(6));
t49 = sin(qJ(5));
t75 = t47 * t49;
t53 = cos(qJ(5));
t74 = t47 * t53;
t73 = t48 * t49;
t72 = t48 * t52;
t71 = t48 * t53;
t70 = t48 * t56;
t41 = pkin(9) + t57;
t69 = t43 * t83;
t68 = t46 * t83;
t67 = t43 * pkin(1) + 0;
t66 = t46 * pkin(1) + 0;
t65 = t43 * t80;
t64 = t44 * t77;
t63 = t47 * t73;
t62 = t47 * t71;
t61 = qJ(1) + 0;
t34 = pkin(5) - t42;
t33 = pkin(5) + t42;
t51 = sin(qJ(3));
t55 = cos(qJ(3));
t60 = pkin(3) * t51 * t56 + (t55 * pkin(3) + pkin(2)) * t52;
t50 = sin(qJ(4));
t54 = cos(qJ(4));
t7 = -t43 * pkin(4) + t48 * t68;
t9 = pkin(4) * t76 + t69;
t59 = pkin(3) * t43 + t9 * t50 - t7 * t54;
t6 = t46 * pkin(4) + t48 * t69;
t8 = pkin(4) * t81 - t68;
t58 = t46 * pkin(3) - t8 * t50 + t6 * t54;
t36 = cos(t42);
t35 = sin(t42);
t30 = cos(t39);
t28 = -qJ(4) + t34;
t27 = qJ(4) + t33;
t26 = cos(t33);
t25 = sin(t34);
t23 = cos(t27);
t22 = sin(t28);
t20 = cos(t34) / 0.2e1;
t19 = sin(t33) / 0.2e1;
t18 = t47 * pkin(10) + t41;
t17 = cos(t28) / 0.2e1;
t16 = sin(t27) / 0.2e1;
t15 = pkin(3) * t36 + t32;
t14 = -t50 * pkin(4) + t54 * t83;
t13 = pkin(4) * t54 + t50 * t83 + pkin(3);
t12 = pkin(2) * t72 - t45 * t57;
t11 = t20 + t26 / 0.2e1;
t10 = t19 - t25 / 0.2e1;
t5 = t17 + t23 / 0.2e1;
t4 = t16 - t22 / 0.2e1;
t3 = -t45 * t41 + t60 * t48;
t2 = pkin(3) * t76 + t7 * t50 + t9 * t54;
t1 = -pkin(3) * t81 - t6 * t50 - t8 * t54;
t21 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; t46, -t43, 0, 0; t43, t46, 0, 0; 0, 0, 1, t61; -t43 * t72 + t46 * t56, -t43 * t70 - t46 * t52, t24, pkin(7) * t24 + t66; t43 * t56 + t46 * t72, -t43 * t52 + t46 * t70, -t77, -pkin(7) * t77 + t67; t79, t45 * t56, t48, t48 * pkin(7) + t61; -t43 * t10 + t46 * t36, -t43 * t11 - t46 * t35, t24, -t43 * t12 + t46 * t32 + 0; t46 * t10 + t43 * t36, t46 * t11 - t43 * t35, -t77, t46 * t12 + t43 * t32 + 0; t20 - t26 / 0.2e1, t19 + t25 / 0.2e1, t48, pkin(2) * t79 + t48 * t57 + t61; t46 * t30 - t43 * t4, -t43 * t5 - t78, t24, t46 * t15 - t43 * t3 + 0; t43 * t30 + t46 * t4, t46 * t5 - t82, -t77, t43 * t15 + t46 * t3 + 0; t17 - t23 / 0.2e1, t16 + t22 / 0.2e1, t48, t48 * t41 + t60 * t45 + t61; (-t43 * t63 + t46 * t53) * t30 + (-t43 * t71 - t46 * t75) * t29 + t49 * t65, (-t43 * t62 - t46 * t49) * t30 + (t43 * t73 - t46 * t74) * t29 + t53 * t65, t44 * t78 + (t30 * t44 * t48 + t45 * t47) * t43, (t46 * pkin(2) + t1 * t51 + t58 * t55) * t56 + (-pkin(2) * t81 + t1 * t55 - t58 * t51) * t52 + t18 * t24 + t66; (t43 * t53 + t46 * t63) * t30 + (-t43 * t75 + t46 * t71) * t29 - t49 * t64, (-t43 * t49 + t46 * t62) * t30 + (-t43 * t74 - t46 * t73) * t29 - t53 * t64, -t47 * t77 + (-t30 * t76 + t82) * t44, (t43 * pkin(2) + t2 * t51 + t59 * t55) * t56 + (pkin(2) * t76 + t2 * t55 - t59 * t51) * t52 - t18 * t77 + t67; t44 * t73 + (t29 * t53 + t30 * t75) * t45, t44 * t71 + (-t29 * t49 + t30 * t74) * t45, -t30 * t80 + t48 * t47, t18 * t48 + ((t13 * t55 + t14 * t51 + pkin(2)) * t52 - (-t51 * t13 + t14 * t55) * t56) * t45 + t61;];
Tc_stack = t21;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), Tc_mdh = NaN(4,4,5+1);               % numerisch
else,                         Tc_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  Tc_mdh(:,:,i) = [Tc_stack((i-1)*3+1 : 3*i, :);[0 0 0 1]];
end
