% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2018-11-23 15:35
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6PRRRRR5_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR5_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR5_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:35:10
% EndTime: 2018-11-23 15:35:11
% DurationCPUTime: 0.31s
% Computational Cost: add. (1656->94), mult. (1692->117), div. (0->0), fcn. (1753->24), ass. (0->76)
t51 = pkin(6) - qJ(2);
t41 = cos(t51) / 0.2e1;
t50 = pkin(6) + qJ(2);
t44 = cos(t50);
t35 = t41 + t44 / 0.2e1;
t53 = sin(pkin(13));
t56 = cos(pkin(13));
t62 = sin(qJ(2));
t26 = -t53 * t35 - t56 * t62;
t54 = sin(pkin(7));
t57 = cos(pkin(7));
t55 = sin(pkin(6));
t96 = t53 * t55;
t80 = -t26 * t54 + t57 * t96;
t40 = sin(t50) / 0.2e1;
t43 = sin(t51);
t32 = t40 + t43 / 0.2e1;
t58 = cos(pkin(6));
t82 = -t32 * t54 + t58 * t57;
t99 = cos(qJ(4));
t95 = t56 * t55;
t93 = pkin(7) - qJ(3);
t92 = pkin(7) + qJ(3);
t90 = qJ(1) + 0;
t59 = sin(qJ(5));
t89 = pkin(5) * t59 + pkin(10);
t88 = t56 * pkin(1) + pkin(8) * t96 + 0;
t87 = t58 * pkin(8) + t90;
t86 = cos(t92);
t85 = sin(t93);
t84 = cos(t93) / 0.2e1;
t83 = sin(t92) / 0.2e1;
t24 = t56 * t35 - t53 * t62;
t81 = -t24 * t54 - t57 * t95;
t79 = t53 * pkin(1) - pkin(8) * t95 + 0;
t33 = t40 - t43 / 0.2e1;
t65 = cos(qJ(2));
t27 = -t53 * t33 + t56 * t65;
t78 = t27 * pkin(2) + pkin(9) * t80 + t88;
t36 = t41 - t44 / 0.2e1;
t77 = t36 * pkin(2) + pkin(9) * t82 + t87;
t31 = t83 - t85 / 0.2e1;
t34 = t84 - t86 / 0.2e1;
t64 = cos(qJ(3));
t12 = t26 * t31 + t27 * t64 + t34 * t96;
t76 = t12 * pkin(3) + t78;
t15 = t32 * t31 + t58 * t34 + t36 * t64;
t75 = t15 * pkin(3) + t77;
t74 = t84 + t86 / 0.2e1;
t73 = t83 + t85 / 0.2e1;
t72 = t55 * t73;
t61 = sin(qJ(3));
t11 = -t26 * t74 + t27 * t61 - t53 * t72;
t71 = t11 * pkin(10) + t76;
t14 = -t32 * t74 + t36 * t61 - t58 * t73;
t70 = t14 * pkin(10) + t75;
t25 = t56 * t33 + t53 * t65;
t69 = t25 * pkin(2) + t81 * pkin(9) + t79;
t10 = t24 * t31 + t25 * t64 - t34 * t95;
t68 = t10 * pkin(3) + t69;
t9 = -t24 * t74 + t25 * t61 + t56 * t72;
t67 = t9 * pkin(10) + t68;
t66 = -pkin(12) - pkin(11);
t63 = cos(qJ(5));
t60 = sin(qJ(4));
t52 = qJ(5) + qJ(6);
t46 = cos(t52);
t45 = sin(t52);
t42 = t63 * pkin(5) + pkin(4);
t6 = t15 * t99 + t82 * t60;
t5 = t15 * t60 - t82 * t99;
t4 = t12 * t99 + t80 * t60;
t3 = t12 * t60 - t80 * t99;
t2 = t10 * t99 + t81 * t60;
t1 = t10 * t60 - t81 * t99;
t7 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t56, -t53, 0, 0; t53, t56, 0, 0; 0, 0, 1, t90; 0, 0, 0, 1; t27, t26, t96, t88; t25, t24, -t95, t79; t36, t32, t58, t87; 0, 0, 0, 1; t12, -t11, t80, t78; t10, -t9, t81, t69; t15, -t14, t82, t77; 0, 0, 0, 1; t4, -t3, t11, t71; t2, -t1, t9, t67; t6, -t5, t14, t70; 0, 0, 0, 1; t11 * t59 + t4 * t63, t11 * t63 - t4 * t59, t3, t4 * pkin(4) + t3 * pkin(11) + t71; t2 * t63 + t9 * t59, -t2 * t59 + t9 * t63, t1, t2 * pkin(4) + t1 * pkin(11) + t67; t14 * t59 + t6 * t63, t14 * t63 - t6 * t59, t5, t6 * pkin(4) + t5 * pkin(11) + t70; 0, 0, 0, 1; t11 * t45 + t4 * t46, t11 * t46 - t4 * t45, t3, t89 * t11 - t3 * t66 + t4 * t42 + t76; t2 * t46 + t9 * t45, -t2 * t45 + t9 * t46, t1, -t1 * t66 + t2 * t42 + t89 * t9 + t68; t14 * t45 + t6 * t46, t14 * t46 - t6 * t45, t5, t89 * t14 + t6 * t42 - t5 * t66 + t75; 0, 0, 0, 1;];
T_ges = t7;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
