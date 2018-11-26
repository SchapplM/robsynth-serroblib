% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
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
% Datum: 2018-11-23 17:02
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RRPRPR4_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR4_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR4_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:01:41
% EndTime: 2018-11-23 17:01:41
% DurationCPUTime: 0.25s
% Computational Cost: add. (640->88), mult. (455->102), div. (0->0), fcn. (490->22), ass. (0->63)
t52 = pkin(6) - qJ(2);
t37 = cos(t52) / 0.2e1;
t51 = pkin(6) + qJ(2);
t47 = cos(t51);
t85 = t37 - t47 / 0.2e1;
t36 = sin(t51) / 0.2e1;
t44 = sin(t52);
t23 = t36 - t44 / 0.2e1;
t54 = cos(pkin(6));
t58 = sin(qJ(4));
t82 = t54 * t58;
t53 = sin(pkin(6));
t60 = sin(qJ(1));
t38 = t60 * t53;
t64 = cos(qJ(1));
t81 = t64 * t53;
t80 = pkin(7) + 0;
t50 = qJ(2) + pkin(11);
t79 = t58 * t38;
t78 = t58 * t81;
t56 = pkin(8) + qJ(3);
t16 = t23 * pkin(2) - t53 * t56;
t63 = cos(qJ(2));
t40 = t63 * pkin(2) + pkin(1);
t77 = t64 * t16 + t60 * t40 + 0;
t76 = pkin(6) - t50;
t75 = pkin(6) + t50;
t74 = -t60 * t16 + t64 * t40 + 0;
t73 = cos(t75);
t72 = sin(t76);
t71 = t85 * pkin(2) + t54 * t56 + t80;
t70 = cos(t76) / 0.2e1;
t69 = sin(t75) / 0.2e1;
t42 = sin(t50);
t65 = t73 / 0.2e1 + t70;
t13 = t64 * t42 + t60 * t65;
t20 = t69 - t72 / 0.2e1;
t46 = cos(t50);
t14 = -t60 * t20 + t64 * t46;
t62 = cos(qJ(4));
t39 = t62 * pkin(4) + pkin(3);
t55 = -qJ(5) - pkin(9);
t68 = pkin(4) * t79 - t13 * t55 + t14 * t39 + t74;
t21 = t72 / 0.2e1 + t69;
t22 = t70 - t73 / 0.2e1;
t67 = pkin(4) * t82 + t21 * t55 + t22 * t39 + t71;
t11 = t60 * t42 - t64 * t65;
t12 = t64 * t20 + t60 * t46;
t66 = -pkin(4) * t78 - t11 * t55 + t12 * t39 + t77;
t61 = cos(qJ(6));
t59 = sin(qJ(2));
t57 = sin(qJ(6));
t49 = qJ(4) + pkin(12);
t45 = cos(t49);
t41 = sin(t49);
t24 = t37 + t47 / 0.2e1;
t10 = t22 * t45 + t54 * t41;
t9 = t22 * t41 - t54 * t45;
t4 = t14 * t45 + t41 * t38;
t3 = t14 * t41 - t45 * t38;
t2 = t12 * t45 - t41 * t81;
t1 = t12 * t41 + t45 * t81;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t64, -t60, 0, 0; t60, t64, 0, 0; 0, 0, 1, t80; 0, 0, 0, 1; -t60 * t23 + t64 * t63, -t60 * t24 - t64 * t59, t38, t64 * pkin(1) + pkin(8) * t38 + 0; t64 * t23 + t60 * t63, t64 * t24 - t60 * t59, -t81, t60 * pkin(1) - pkin(8) * t81 + 0; t85, t36 + t44 / 0.2e1, t54, t54 * pkin(8) + t80; 0, 0, 0, 1; t14, -t13, t38, t74; t12, -t11, -t81, t77; t22, t21, t54, t71; 0, 0, 0, 1; t14 * t62 + t79, -t14 * t58 + t62 * t38, t13, t14 * pkin(3) + t13 * pkin(9) + t74; t12 * t62 - t78, -t12 * t58 - t62 * t81, t11, t12 * pkin(3) + t11 * pkin(9) + t77; t22 * t62 + t82, -t22 * t58 + t54 * t62, -t21, t22 * pkin(3) - t21 * pkin(9) + t71; 0, 0, 0, 1; t4, -t3, t13, t68; t2, -t1, t11, t66; t10, -t9, -t21, t67; 0, 0, 0, 1; t13 * t57 + t4 * t61, t13 * t61 - t4 * t57, t3, t4 * pkin(5) + t3 * pkin(10) + t68; t11 * t57 + t2 * t61, t11 * t61 - t2 * t57, t1, t2 * pkin(5) + t1 * pkin(10) + t66; t10 * t61 - t21 * t57, -t10 * t57 - t21 * t61, t9, t10 * pkin(5) + t9 * pkin(10) + t67; 0, 0, 0, 1;];
T_ges = t5;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
