% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
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
% Datum: 2018-11-23 15:03
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6PRPRRR1_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR1_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR1_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:02:58
% EndTime: 2018-11-23 15:02:59
% DurationCPUTime: 0.24s
% Computational Cost: add. (640->88), mult. (455->102), div. (0->0), fcn. (490->22), ass. (0->63)
t51 = pkin(6) - qJ(2);
t38 = cos(t51) / 0.2e1;
t50 = pkin(6) + qJ(2);
t45 = cos(t50);
t85 = t38 - t45 / 0.2e1;
t37 = sin(t50) / 0.2e1;
t43 = sin(t51);
t23 = t37 - t43 / 0.2e1;
t53 = sin(pkin(11));
t54 = sin(pkin(6));
t35 = t53 * t54;
t55 = cos(pkin(11));
t82 = t55 * t54;
t56 = cos(pkin(6));
t59 = sin(qJ(4));
t81 = t56 * t59;
t49 = qJ(2) + pkin(12);
t80 = t59 * t35;
t79 = t59 * t82;
t78 = qJ(1) + 0;
t57 = pkin(7) + qJ(3);
t16 = t23 * pkin(2) - t54 * t57;
t63 = cos(qJ(2));
t40 = t63 * pkin(2) + pkin(1);
t77 = t55 * t16 + t53 * t40 + 0;
t76 = pkin(6) - t49;
t75 = pkin(6) + t49;
t74 = -t53 * t16 + t55 * t40 + 0;
t73 = cos(t75);
t72 = sin(t76);
t71 = t85 * pkin(2) + t56 * t57 + t78;
t70 = cos(t76) / 0.2e1;
t69 = sin(t75) / 0.2e1;
t41 = sin(t49);
t65 = t73 / 0.2e1 + t70;
t13 = t55 * t41 + t53 * t65;
t20 = t69 - t72 / 0.2e1;
t44 = cos(t49);
t14 = -t53 * t20 + t55 * t44;
t62 = cos(qJ(4));
t39 = t62 * pkin(4) + pkin(3);
t64 = -pkin(9) - pkin(8);
t68 = pkin(4) * t80 - t13 * t64 + t14 * t39 + t74;
t21 = t72 / 0.2e1 + t69;
t22 = t70 - t73 / 0.2e1;
t67 = pkin(4) * t81 + t21 * t64 + t22 * t39 + t71;
t11 = t53 * t41 - t55 * t65;
t12 = t55 * t20 + t53 * t44;
t66 = -pkin(4) * t79 - t11 * t64 + t12 * t39 + t77;
t61 = cos(qJ(6));
t60 = sin(qJ(2));
t58 = sin(qJ(6));
t52 = qJ(4) + qJ(5);
t48 = cos(t52);
t47 = sin(t52);
t24 = t38 + t45 / 0.2e1;
t10 = t22 * t48 + t56 * t47;
t9 = t22 * t47 - t56 * t48;
t4 = t14 * t48 + t47 * t35;
t3 = t14 * t47 - t48 * t35;
t2 = t12 * t48 - t47 * t82;
t1 = t12 * t47 + t48 * t82;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t55, -t53, 0, 0; t53, t55, 0, 0; 0, 0, 1, t78; 0, 0, 0, 1; -t53 * t23 + t55 * t63, -t53 * t24 - t55 * t60, t35, t55 * pkin(1) + pkin(7) * t35 + 0; t55 * t23 + t53 * t63, t55 * t24 - t53 * t60, -t82, t53 * pkin(1) - pkin(7) * t82 + 0; t85, t37 + t43 / 0.2e1, t56, t56 * pkin(7) + t78; 0, 0, 0, 1; t14, -t13, t35, t74; t12, -t11, -t82, t77; t22, t21, t56, t71; 0, 0, 0, 1; t14 * t62 + t80, -t14 * t59 + t62 * t35, t13, t14 * pkin(3) + t13 * pkin(8) + t74; t12 * t62 - t79, -t12 * t59 - t62 * t82, t11, t12 * pkin(3) + t11 * pkin(8) + t77; t22 * t62 + t81, -t22 * t59 + t56 * t62, -t21, t22 * pkin(3) - t21 * pkin(8) + t71; 0, 0, 0, 1; t4, -t3, t13, t68; t2, -t1, t11, t66; t10, -t9, -t21, t67; 0, 0, 0, 1; t13 * t58 + t4 * t61, t13 * t61 - t4 * t58, t3, t4 * pkin(5) + t3 * pkin(10) + t68; t11 * t58 + t2 * t61, t11 * t61 - t2 * t58, t1, t2 * pkin(5) + t1 * pkin(10) + t66; t10 * t61 - t21 * t58, -t10 * t58 - t21 * t61, t9, t10 * pkin(5) + t9 * pkin(10) + t67; 0, 0, 0, 1;];
T_ges = t5;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
