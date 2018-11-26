% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3,theta4]';
% 
% Output:
% T_c_mdh [4x4x(6+1)]
%   homogenous transformation matrices for each (body) frame (MDH)
%   1:  mdh base (link 0) -> mdh base link 0 (unit matrix, no information)
%   ...
%   7:  mdh base (link 0) -> mdh frame (7-1), link (7-1)
%   ...
%   6+1:  mdh base (link 0) -> mdh frame (6)

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 14:53
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6PRPPRR1_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR1_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPPRR1_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 14:53:27
% EndTime: 2018-11-23 14:53:27
% DurationCPUTime: 0.27s
% Computational Cost: add. (640->88), mult. (455->102), div. (0->0), fcn. (490->22), ass. (0->63)
t52 = pkin(6) - qJ(2);
t38 = cos(t52) / 0.2e1;
t51 = pkin(6) + qJ(2);
t47 = cos(t51);
t85 = t38 - t47 / 0.2e1;
t37 = sin(t51) / 0.2e1;
t44 = sin(t52);
t23 = t37 - t44 / 0.2e1;
t54 = sin(pkin(10));
t55 = sin(pkin(6));
t35 = t54 * t55;
t57 = cos(pkin(10));
t82 = t57 * t55;
t53 = sin(pkin(12));
t58 = cos(pkin(6));
t81 = t58 * t53;
t50 = qJ(2) + pkin(11);
t80 = t53 * t35;
t79 = t53 * t82;
t78 = qJ(1) + 0;
t60 = pkin(7) + qJ(3);
t17 = t23 * pkin(2) - t55 * t60;
t64 = cos(qJ(2));
t40 = t64 * pkin(2) + pkin(1);
t77 = t57 * t17 + t54 * t40 + 0;
t76 = pkin(6) - t50;
t75 = pkin(6) + t50;
t74 = -t54 * t17 + t57 * t40 + 0;
t73 = cos(t75);
t72 = sin(t76);
t71 = t85 * pkin(2) + t58 * t60 + t78;
t70 = cos(t76) / 0.2e1;
t69 = sin(t75) / 0.2e1;
t42 = sin(t50);
t65 = t73 / 0.2e1 + t70;
t13 = t57 * t42 + t54 * t65;
t20 = t69 - t72 / 0.2e1;
t46 = cos(t50);
t14 = -t54 * t20 + t57 * t46;
t56 = cos(pkin(12));
t39 = t56 * pkin(4) + pkin(3);
t59 = -pkin(8) - qJ(4);
t68 = pkin(4) * t80 - t13 * t59 + t14 * t39 + t74;
t21 = t72 / 0.2e1 + t69;
t22 = t70 - t73 / 0.2e1;
t67 = pkin(4) * t81 + t21 * t59 + t22 * t39 + t71;
t11 = t54 * t42 - t57 * t65;
t12 = t57 * t20 + t54 * t46;
t66 = -pkin(4) * t79 - t11 * t59 + t12 * t39 + t77;
t63 = cos(qJ(6));
t62 = sin(qJ(2));
t61 = sin(qJ(6));
t49 = pkin(12) + qJ(5);
t45 = cos(t49);
t41 = sin(t49);
t24 = t38 + t47 / 0.2e1;
t10 = t22 * t45 + t58 * t41;
t9 = t22 * t41 - t58 * t45;
t4 = t14 * t45 + t41 * t35;
t3 = t14 * t41 - t45 * t35;
t2 = t12 * t45 - t41 * t82;
t1 = t12 * t41 + t45 * t82;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t57, -t54, 0, 0; t54, t57, 0, 0; 0, 0, 1, t78; 0, 0, 0, 1; -t54 * t23 + t57 * t64, -t54 * t24 - t57 * t62, t35, t57 * pkin(1) + pkin(7) * t35 + 0; t57 * t23 + t54 * t64, t57 * t24 - t54 * t62, -t82, t54 * pkin(1) - pkin(7) * t82 + 0; t85, t37 + t44 / 0.2e1, t58, t58 * pkin(7) + t78; 0, 0, 0, 1; t14, -t13, t35, t74; t12, -t11, -t82, t77; t22, t21, t58, t71; 0, 0, 0, 1; t14 * t56 + t80, -t14 * t53 + t56 * t35, t13, t14 * pkin(3) + t13 * qJ(4) + t74; t12 * t56 - t79, -t12 * t53 - t56 * t82, t11, t12 * pkin(3) + t11 * qJ(4) + t77; t22 * t56 + t81, -t22 * t53 + t58 * t56, -t21, t22 * pkin(3) - t21 * qJ(4) + t71; 0, 0, 0, 1; t4, -t3, t13, t68; t2, -t1, t11, t66; t10, -t9, -t21, t67; 0, 0, 0, 1; t13 * t61 + t4 * t63, t13 * t63 - t4 * t61, t3, t4 * pkin(5) + t3 * pkin(9) + t68; t11 * t61 + t2 * t63, t11 * t63 - t2 * t61, t1, t2 * pkin(5) + t1 * pkin(9) + t66; t10 * t63 - t21 * t61, -t10 * t61 - t21 * t63, t9, t10 * pkin(5) + t9 * pkin(9) + t67; 0, 0, 0, 1;];
T_ges = t5;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
