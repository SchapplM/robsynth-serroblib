% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PPRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2,theta5]';
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
% Datum: 2018-11-23 14:50
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6PPRRPR1_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR1_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRPR1_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 14:49:42
% EndTime: 2018-11-23 14:49:43
% DurationCPUTime: 0.32s
% Computational Cost: add. (1656->94), mult. (1692->118), div. (0->0), fcn. (1753->24), ass. (0->77)
t51 = pkin(6) - pkin(12);
t41 = cos(t51) / 0.2e1;
t50 = pkin(6) + pkin(12);
t44 = cos(t50);
t33 = t41 + t44 / 0.2e1;
t54 = sin(pkin(12));
t55 = sin(pkin(11));
t60 = cos(pkin(11));
t26 = -t55 * t33 - t60 * t54;
t56 = sin(pkin(7));
t61 = cos(pkin(7));
t57 = sin(pkin(6));
t97 = t55 * t57;
t80 = -t26 * t56 + t61 * t97;
t40 = sin(t50) / 0.2e1;
t43 = sin(t51);
t31 = t40 + t43 / 0.2e1;
t62 = cos(pkin(6));
t82 = -t31 * t56 + t62 * t61;
t100 = cos(qJ(4));
t96 = t60 * t57;
t94 = qJ(2) * t57;
t93 = pkin(7) - qJ(3);
t92 = pkin(7) + qJ(3);
t90 = qJ(1) + 0;
t53 = sin(pkin(13));
t89 = pkin(5) * t53 + pkin(9);
t88 = t60 * pkin(1) + t55 * t94 + 0;
t87 = t62 * qJ(2) + t90;
t86 = cos(t92);
t85 = sin(t93);
t84 = cos(t93) / 0.2e1;
t83 = sin(t92) / 0.2e1;
t24 = t60 * t33 - t55 * t54;
t81 = -t24 * t56 - t61 * t96;
t79 = t55 * pkin(1) - t60 * t94 + 0;
t32 = t40 - t43 / 0.2e1;
t59 = cos(pkin(12));
t27 = -t55 * t32 + t60 * t59;
t78 = t27 * pkin(2) + pkin(8) * t80 + t88;
t34 = t41 - t44 / 0.2e1;
t77 = t34 * pkin(2) + pkin(8) * t82 + t87;
t35 = t83 - t85 / 0.2e1;
t36 = t84 - t86 / 0.2e1;
t66 = cos(qJ(3));
t12 = t26 * t35 + t27 * t66 + t36 * t97;
t76 = t12 * pkin(3) + t78;
t15 = t31 * t35 + t34 * t66 + t62 * t36;
t75 = t15 * pkin(3) + t77;
t74 = t84 + t86 / 0.2e1;
t73 = t83 + t85 / 0.2e1;
t72 = t57 * t73;
t65 = sin(qJ(3));
t11 = -t26 * t74 + t27 * t65 - t55 * t72;
t71 = t11 * pkin(9) + t76;
t14 = -t31 * t74 + t34 * t65 - t62 * t73;
t70 = t14 * pkin(9) + t75;
t25 = t60 * t32 + t55 * t59;
t69 = t25 * pkin(2) + t81 * pkin(8) + t79;
t10 = t24 * t35 + t25 * t66 - t36 * t96;
t68 = t10 * pkin(3) + t69;
t9 = -t24 * t74 + t25 * t65 + t60 * t72;
t67 = t9 * pkin(9) + t68;
t64 = sin(qJ(4));
t63 = -pkin(10) - qJ(5);
t58 = cos(pkin(13));
t52 = pkin(13) + qJ(6);
t46 = cos(t52);
t45 = sin(t52);
t42 = t58 * pkin(5) + pkin(4);
t6 = t15 * t100 + t82 * t64;
t5 = -t82 * t100 + t15 * t64;
t4 = t12 * t100 + t80 * t64;
t3 = -t80 * t100 + t12 * t64;
t2 = t10 * t100 + t81 * t64;
t1 = t10 * t64 - t81 * t100;
t7 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t60, -t55, 0, 0; t55, t60, 0, 0; 0, 0, 1, t90; 0, 0, 0, 1; t27, t26, t97, t88; t25, t24, -t96, t79; t34, t31, t62, t87; 0, 0, 0, 1; t12, -t11, t80, t78; t10, -t9, t81, t69; t15, -t14, t82, t77; 0, 0, 0, 1; t4, -t3, t11, t71; t2, -t1, t9, t67; t6, -t5, t14, t70; 0, 0, 0, 1; t11 * t53 + t4 * t58, t11 * t58 - t4 * t53, t3, t4 * pkin(4) + t3 * qJ(5) + t71; t2 * t58 + t9 * t53, -t2 * t53 + t9 * t58, t1, t2 * pkin(4) + t1 * qJ(5) + t67; t14 * t53 + t6 * t58, t14 * t58 - t6 * t53, t5, t6 * pkin(4) + t5 * qJ(5) + t70; 0, 0, 0, 1; t11 * t45 + t4 * t46, t11 * t46 - t4 * t45, t3, t89 * t11 - t3 * t63 + t4 * t42 + t76; t2 * t46 + t9 * t45, -t2 * t45 + t9 * t46, t1, -t1 * t63 + t2 * t42 + t89 * t9 + t68; t14 * t45 + t6 * t46, t14 * t46 - t6 * t45, t5, t89 * t14 + t6 * t42 - t5 * t63 + t75; 0, 0, 0, 1;];
T_ges = t7;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
