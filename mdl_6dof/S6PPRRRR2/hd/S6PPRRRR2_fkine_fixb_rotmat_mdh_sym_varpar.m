% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PPRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
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
% Datum: 2018-11-23 14:52
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6PPRRRR2_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR2_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR2_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 14:52:13
% EndTime: 2018-11-23 14:52:13
% DurationCPUTime: 0.32s
% Computational Cost: add. (1656->94), mult. (1692->118), div. (0->0), fcn. (1753->24), ass. (0->77)
t51 = pkin(6) - pkin(13);
t41 = cos(t51) / 0.2e1;
t50 = pkin(6) + pkin(13);
t44 = cos(t50);
t33 = t41 + t44 / 0.2e1;
t53 = sin(pkin(13));
t54 = sin(pkin(12));
t58 = cos(pkin(12));
t26 = -t54 * t33 - t58 * t53;
t55 = sin(pkin(7));
t59 = cos(pkin(7));
t56 = sin(pkin(6));
t97 = t54 * t56;
t80 = -t26 * t55 + t59 * t97;
t40 = sin(t50) / 0.2e1;
t43 = sin(t51);
t31 = t40 + t43 / 0.2e1;
t60 = cos(pkin(6));
t82 = -t31 * t55 + t60 * t59;
t100 = cos(qJ(4));
t96 = t58 * t56;
t94 = qJ(2) * t56;
t93 = pkin(7) - qJ(3);
t92 = pkin(7) + qJ(3);
t90 = qJ(1) + 0;
t61 = sin(qJ(5));
t89 = pkin(5) * t61 + pkin(9);
t88 = t58 * pkin(1) + t54 * t94 + 0;
t87 = t60 * qJ(2) + t90;
t86 = cos(t92);
t85 = sin(t93);
t84 = cos(t93) / 0.2e1;
t83 = sin(t92) / 0.2e1;
t24 = t58 * t33 - t54 * t53;
t81 = -t24 * t55 - t59 * t96;
t79 = t54 * pkin(1) - t58 * t94 + 0;
t32 = t40 - t43 / 0.2e1;
t57 = cos(pkin(13));
t27 = -t54 * t32 + t58 * t57;
t78 = t27 * pkin(2) + pkin(8) * t80 + t88;
t34 = t41 - t44 / 0.2e1;
t77 = t34 * pkin(2) + pkin(8) * t82 + t87;
t35 = t83 - t85 / 0.2e1;
t36 = t84 - t86 / 0.2e1;
t65 = cos(qJ(3));
t12 = t26 * t35 + t27 * t65 + t36 * t97;
t76 = t12 * pkin(3) + t78;
t15 = t31 * t35 + t34 * t65 + t60 * t36;
t75 = t15 * pkin(3) + t77;
t74 = t84 + t86 / 0.2e1;
t73 = t83 + t85 / 0.2e1;
t72 = t56 * t73;
t63 = sin(qJ(3));
t11 = -t26 * t74 + t27 * t63 - t54 * t72;
t71 = t11 * pkin(9) + t76;
t14 = -t31 * t74 + t34 * t63 - t60 * t73;
t70 = t14 * pkin(9) + t75;
t25 = t58 * t32 + t54 * t57;
t69 = t25 * pkin(2) + t81 * pkin(8) + t79;
t10 = t24 * t35 + t25 * t65 - t36 * t96;
t68 = t10 * pkin(3) + t69;
t9 = -t24 * t74 + t25 * t63 + t58 * t72;
t67 = t9 * pkin(9) + t68;
t66 = -pkin(11) - pkin(10);
t64 = cos(qJ(5));
t62 = sin(qJ(4));
t52 = qJ(5) + qJ(6);
t47 = cos(t52);
t46 = sin(t52);
t42 = t64 * pkin(5) + pkin(4);
t6 = t15 * t100 + t82 * t62;
t5 = -t82 * t100 + t15 * t62;
t4 = t12 * t100 + t80 * t62;
t3 = -t80 * t100 + t12 * t62;
t2 = t10 * t100 + t81 * t62;
t1 = t10 * t62 - t81 * t100;
t7 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t58, -t54, 0, 0; t54, t58, 0, 0; 0, 0, 1, t90; 0, 0, 0, 1; t27, t26, t97, t88; t25, t24, -t96, t79; t34, t31, t60, t87; 0, 0, 0, 1; t12, -t11, t80, t78; t10, -t9, t81, t69; t15, -t14, t82, t77; 0, 0, 0, 1; t4, -t3, t11, t71; t2, -t1, t9, t67; t6, -t5, t14, t70; 0, 0, 0, 1; t11 * t61 + t4 * t64, t11 * t64 - t4 * t61, t3, t4 * pkin(4) + t3 * pkin(10) + t71; t2 * t64 + t9 * t61, -t2 * t61 + t9 * t64, t1, t2 * pkin(4) + t1 * pkin(10) + t67; t14 * t61 + t6 * t64, t14 * t64 - t6 * t61, t5, t6 * pkin(4) + t5 * pkin(10) + t70; 0, 0, 0, 1; t11 * t46 + t4 * t47, t11 * t47 - t4 * t46, t3, t89 * t11 - t3 * t66 + t4 * t42 + t76; t2 * t47 + t9 * t46, -t2 * t46 + t9 * t47, t1, -t1 * t66 + t2 * t42 + t89 * t9 + t68; t14 * t46 + t6 * t47, t14 * t47 - t6 * t46, t5, t89 * t14 + t6 * t42 - t5 * t66 + t75; 0, 0, 0, 1;];
T_ges = t7;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
