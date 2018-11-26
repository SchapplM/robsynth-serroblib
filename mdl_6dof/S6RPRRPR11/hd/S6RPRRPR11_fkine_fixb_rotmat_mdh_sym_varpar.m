% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
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
% Datum: 2018-11-23 16:22
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RPRRPR11_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR11_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR11_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:22:04
% EndTime: 2018-11-23 16:22:04
% DurationCPUTime: 0.31s
% Computational Cost: add. (1656->94), mult. (1692->118), div. (0->0), fcn. (1753->24), ass. (0->77)
t51 = pkin(6) - pkin(12);
t41 = cos(t51) / 0.2e1;
t50 = pkin(6) + pkin(12);
t44 = cos(t50);
t33 = t41 + t44 / 0.2e1;
t54 = sin(pkin(12));
t64 = sin(qJ(1));
t66 = cos(qJ(1));
t26 = -t64 * t33 - t66 * t54;
t55 = sin(pkin(7));
t59 = cos(pkin(7));
t56 = sin(pkin(6));
t96 = t64 * t56;
t80 = -t26 * t55 + t59 * t96;
t40 = sin(t50) / 0.2e1;
t43 = sin(t51);
t31 = t40 + t43 / 0.2e1;
t60 = cos(pkin(6));
t82 = -t31 * t55 + t60 * t59;
t100 = cos(qJ(4));
t95 = t66 * t56;
t94 = qJ(2) * t56;
t93 = pkin(7) - qJ(3);
t92 = pkin(7) + qJ(3);
t91 = pkin(8) + 0;
t53 = sin(pkin(13));
t89 = pkin(5) * t53 + pkin(10);
t88 = t60 * qJ(2) + t91;
t87 = t66 * pkin(1) + t64 * t94 + 0;
t86 = cos(t92);
t85 = sin(t93);
t84 = cos(t93) / 0.2e1;
t83 = sin(t92) / 0.2e1;
t24 = t66 * t33 - t64 * t54;
t81 = -t24 * t55 - t59 * t95;
t79 = t64 * pkin(1) - t66 * t94 + 0;
t34 = t41 - t44 / 0.2e1;
t78 = t34 * pkin(2) + pkin(9) * t82 + t88;
t32 = t40 - t43 / 0.2e1;
t58 = cos(pkin(12));
t27 = -t64 * t32 + t66 * t58;
t77 = t27 * pkin(2) + pkin(9) * t80 + t87;
t35 = t83 - t85 / 0.2e1;
t36 = t84 - t86 / 0.2e1;
t65 = cos(qJ(3));
t15 = t31 * t35 + t34 * t65 + t60 * t36;
t76 = t15 * pkin(3) + t78;
t12 = t26 * t35 + t27 * t65 + t36 * t96;
t75 = t12 * pkin(3) + t77;
t74 = t84 + t86 / 0.2e1;
t73 = t83 + t85 / 0.2e1;
t72 = t56 * t73;
t63 = sin(qJ(3));
t14 = -t31 * t74 + t34 * t63 - t60 * t73;
t71 = t14 * pkin(10) + t76;
t11 = -t26 * t74 + t27 * t63 - t64 * t72;
t70 = t11 * pkin(10) + t75;
t25 = t66 * t32 + t64 * t58;
t69 = t25 * pkin(2) + t81 * pkin(9) + t79;
t10 = t24 * t35 + t25 * t65 - t36 * t95;
t68 = t10 * pkin(3) + t69;
t9 = -t24 * t74 + t25 * t63 + t66 * t72;
t67 = t9 * pkin(10) + t68;
t62 = sin(qJ(4));
t61 = -pkin(11) - qJ(5);
t57 = cos(pkin(13));
t52 = pkin(13) + qJ(6);
t46 = cos(t52);
t45 = sin(t52);
t42 = t57 * pkin(5) + pkin(4);
t6 = t15 * t100 + t82 * t62;
t5 = -t82 * t100 + t15 * t62;
t4 = t12 * t100 + t80 * t62;
t3 = -t80 * t100 + t12 * t62;
t2 = t10 * t100 + t81 * t62;
t1 = t10 * t62 - t81 * t100;
t7 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t66, -t64, 0, 0; t64, t66, 0, 0; 0, 0, 1, t91; 0, 0, 0, 1; t27, t26, t96, t87; t25, t24, -t95, t79; t34, t31, t60, t88; 0, 0, 0, 1; t12, -t11, t80, t77; t10, -t9, t81, t69; t15, -t14, t82, t78; 0, 0, 0, 1; t4, -t3, t11, t70; t2, -t1, t9, t67; t6, -t5, t14, t71; 0, 0, 0, 1; t11 * t53 + t4 * t57, t11 * t57 - t4 * t53, t3, t4 * pkin(4) + t3 * qJ(5) + t70; t2 * t57 + t9 * t53, -t2 * t53 + t9 * t57, t1, t2 * pkin(4) + t1 * qJ(5) + t67; t14 * t53 + t6 * t57, t14 * t57 - t6 * t53, t5, t6 * pkin(4) + t5 * qJ(5) + t71; 0, 0, 0, 1; t11 * t45 + t4 * t46, t11 * t46 - t4 * t45, t3, t89 * t11 - t3 * t61 + t4 * t42 + t75; t2 * t46 + t9 * t45, -t2 * t45 + t9 * t46, t1, -t1 * t61 + t2 * t42 + t89 * t9 + t68; t14 * t45 + t6 * t46, t14 * t46 - t6 * t45, t5, t89 * t14 + t6 * t42 - t5 * t61 + t76; 0, 0, 0, 1;];
T_ges = t7;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
