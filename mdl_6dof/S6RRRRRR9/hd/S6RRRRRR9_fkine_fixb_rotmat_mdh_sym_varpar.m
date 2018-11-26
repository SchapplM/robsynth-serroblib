% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2018-11-23 18:45
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RRRRRR9_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR9_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR9_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:44:56
% EndTime: 2018-11-23 18:44:56
% DurationCPUTime: 0.30s
% Computational Cost: add. (1656->94), mult. (1692->117), div. (0->0), fcn. (1753->24), ass. (0->76)
t51 = pkin(6) - qJ(2);
t41 = cos(t51) / 0.2e1;
t50 = pkin(6) + qJ(2);
t44 = cos(t50);
t35 = t41 + t44 / 0.2e1;
t60 = sin(qJ(2));
t61 = sin(qJ(1));
t65 = cos(qJ(1));
t26 = -t61 * t35 - t65 * t60;
t53 = sin(pkin(7));
t55 = cos(pkin(7));
t54 = sin(pkin(6));
t95 = t61 * t54;
t80 = -t26 * t53 + t55 * t95;
t40 = sin(t50) / 0.2e1;
t43 = sin(t51);
t32 = t40 + t43 / 0.2e1;
t56 = cos(pkin(6));
t82 = -t32 * t53 + t56 * t55;
t99 = cos(qJ(4));
t94 = t65 * t54;
t93 = pkin(7) - qJ(3);
t92 = pkin(7) + qJ(3);
t91 = pkin(8) + 0;
t57 = sin(qJ(5));
t89 = pkin(5) * t57 + pkin(11);
t88 = t56 * pkin(9) + t91;
t87 = t65 * pkin(1) + pkin(9) * t95 + 0;
t86 = cos(t92);
t85 = sin(t93);
t84 = cos(t93) / 0.2e1;
t83 = sin(t92) / 0.2e1;
t24 = t65 * t35 - t61 * t60;
t81 = -t24 * t53 - t55 * t94;
t79 = t61 * pkin(1) - pkin(9) * t94 + 0;
t36 = t41 - t44 / 0.2e1;
t78 = t36 * pkin(2) + pkin(10) * t82 + t88;
t33 = t40 - t43 / 0.2e1;
t64 = cos(qJ(2));
t27 = -t61 * t33 + t65 * t64;
t77 = t27 * pkin(2) + pkin(10) * t80 + t87;
t31 = t83 - t85 / 0.2e1;
t34 = t84 - t86 / 0.2e1;
t63 = cos(qJ(3));
t15 = t32 * t31 + t56 * t34 + t36 * t63;
t76 = t15 * pkin(3) + t78;
t12 = t26 * t31 + t27 * t63 + t34 * t95;
t75 = t12 * pkin(3) + t77;
t74 = t84 + t86 / 0.2e1;
t73 = t83 + t85 / 0.2e1;
t72 = t54 * t73;
t59 = sin(qJ(3));
t14 = -t32 * t74 + t36 * t59 - t56 * t73;
t71 = t14 * pkin(11) + t76;
t11 = -t26 * t74 + t27 * t59 - t61 * t72;
t70 = t11 * pkin(11) + t75;
t25 = t65 * t33 + t61 * t64;
t69 = t25 * pkin(2) + t81 * pkin(10) + t79;
t10 = t24 * t31 + t25 * t63 - t34 * t94;
t68 = t10 * pkin(3) + t69;
t9 = -t24 * t74 + t25 * t59 + t65 * t72;
t67 = t9 * pkin(11) + t68;
t66 = -pkin(13) - pkin(12);
t62 = cos(qJ(5));
t58 = sin(qJ(4));
t52 = qJ(5) + qJ(6);
t46 = cos(t52);
t45 = sin(t52);
t42 = t62 * pkin(5) + pkin(4);
t6 = t15 * t99 + t82 * t58;
t5 = t15 * t58 - t82 * t99;
t4 = t12 * t99 + t80 * t58;
t3 = t12 * t58 - t80 * t99;
t2 = t10 * t99 + t81 * t58;
t1 = t10 * t58 - t81 * t99;
t7 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t65, -t61, 0, 0; t61, t65, 0, 0; 0, 0, 1, t91; 0, 0, 0, 1; t27, t26, t95, t87; t25, t24, -t94, t79; t36, t32, t56, t88; 0, 0, 0, 1; t12, -t11, t80, t77; t10, -t9, t81, t69; t15, -t14, t82, t78; 0, 0, 0, 1; t4, -t3, t11, t70; t2, -t1, t9, t67; t6, -t5, t14, t71; 0, 0, 0, 1; t11 * t57 + t4 * t62, t11 * t62 - t4 * t57, t3, t4 * pkin(4) + t3 * pkin(12) + t70; t2 * t62 + t9 * t57, -t2 * t57 + t9 * t62, t1, t2 * pkin(4) + t1 * pkin(12) + t67; t14 * t57 + t6 * t62, t14 * t62 - t6 * t57, t5, t6 * pkin(4) + t5 * pkin(12) + t71; 0, 0, 0, 1; t11 * t45 + t4 * t46, t11 * t46 - t4 * t45, t3, t89 * t11 - t3 * t66 + t4 * t42 + t75; t2 * t46 + t9 * t45, -t2 * t45 + t9 * t46, t1, -t1 * t66 + t2 * t42 + t89 * t9 + t68; t14 * t45 + t6 * t46, t14 * t46 - t6 * t45, t5, t89 * t14 + t6 * t42 - t5 * t66 + t76; 0, 0, 0, 1;];
T_ges = t7;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
