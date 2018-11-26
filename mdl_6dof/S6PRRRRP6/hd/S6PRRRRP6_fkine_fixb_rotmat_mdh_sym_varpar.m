% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,theta1]';
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
% Datum: 2018-11-23 15:31
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6PRRRRP6_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP6_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP6_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:31:08
% EndTime: 2018-11-23 15:31:09
% DurationCPUTime: 0.31s
% Computational Cost: add. (1819->84), mult. (1878->100), div. (0->0), fcn. (1949->22), ass. (0->76)
t61 = sin(pkin(12));
t63 = sin(pkin(6));
t102 = t61 * t63;
t60 = pkin(6) - qJ(2);
t53 = cos(t60) / 0.2e1;
t59 = pkin(6) + qJ(2);
t55 = cos(t59);
t47 = t53 + t55 / 0.2e1;
t64 = cos(pkin(12));
t70 = sin(qJ(2));
t38 = -t61 * t47 - t64 * t70;
t62 = sin(pkin(7));
t65 = cos(pkin(7));
t87 = t65 * t102 - t38 * t62;
t52 = sin(t59) / 0.2e1;
t54 = sin(t60);
t44 = t52 + t54 / 0.2e1;
t66 = cos(pkin(6));
t89 = -t44 * t62 + t66 * t65;
t105 = cos(qJ(4));
t101 = t64 * t63;
t99 = pkin(7) - qJ(3);
t98 = pkin(7) + qJ(3);
t96 = qJ(1) + 0;
t95 = t64 * pkin(1) + pkin(8) * t102 + 0;
t94 = t66 * pkin(8) + t96;
t93 = cos(t98);
t92 = sin(t99);
t91 = cos(t99) / 0.2e1;
t90 = sin(t98) / 0.2e1;
t36 = t64 * t47 - t61 * t70;
t88 = -t65 * t101 - t36 * t62;
t86 = t61 * pkin(1) - pkin(8) * t101 + 0;
t45 = t52 - t54 / 0.2e1;
t73 = cos(qJ(2));
t39 = -t61 * t45 + t64 * t73;
t85 = t39 * pkin(2) + t87 * pkin(9) + t95;
t48 = t53 - t55 / 0.2e1;
t84 = t48 * pkin(2) + t89 * pkin(9) + t94;
t83 = t91 + t93 / 0.2e1;
t82 = t90 + t92 / 0.2e1;
t81 = t63 * t82;
t69 = sin(qJ(3));
t22 = -t38 * t83 + t39 * t69 - t61 * t81;
t43 = t90 - t92 / 0.2e1;
t46 = t91 - t93 / 0.2e1;
t72 = cos(qJ(3));
t23 = t46 * t102 + t38 * t43 + t39 * t72;
t80 = t23 * pkin(3) + t22 * pkin(10) + t85;
t26 = -t44 * t83 + t48 * t69 - t66 * t82;
t27 = t44 * t43 + t66 * t46 + t48 * t72;
t79 = t27 * pkin(3) + t26 * pkin(10) + t84;
t37 = t64 * t45 + t61 * t73;
t78 = t37 * pkin(2) + t88 * pkin(9) + t86;
t68 = sin(qJ(4));
t11 = -t87 * t105 + t23 * t68;
t12 = t23 * t105 + t87 * t68;
t77 = t12 * pkin(4) + t11 * pkin(11) + t80;
t14 = -t89 * t105 + t27 * t68;
t15 = t27 * t105 + t89 * t68;
t76 = t15 * pkin(4) + t14 * pkin(11) + t79;
t20 = -t36 * t83 + t37 * t69 + t64 * t81;
t21 = -t46 * t101 + t36 * t43 + t37 * t72;
t75 = t21 * pkin(3) + t20 * pkin(10) + t78;
t10 = t21 * t105 + t88 * t68;
t9 = -t88 * t105 + t21 * t68;
t74 = t10 * pkin(4) + t9 * pkin(11) + t75;
t71 = cos(qJ(5));
t67 = sin(qJ(5));
t6 = t15 * t71 + t26 * t67;
t5 = t15 * t67 - t26 * t71;
t4 = t12 * t71 + t22 * t67;
t3 = t12 * t67 - t22 * t71;
t2 = t10 * t71 + t20 * t67;
t1 = t10 * t67 - t20 * t71;
t7 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t64, -t61, 0, 0; t61, t64, 0, 0; 0, 0, 1, t96; 0, 0, 0, 1; t39, t38, t102, t95; t37, t36, -t101, t86; t48, t44, t66, t94; 0, 0, 0, 1; t23, -t22, t87, t85; t21, -t20, t88, t78; t27, -t26, t89, t84; 0, 0, 0, 1; t12, -t11, t22, t80; t10, -t9, t20, t75; t15, -t14, t26, t79; 0, 0, 0, 1; t4, -t3, t11, t77; t2, -t1, t9, t74; t6, -t5, t14, t76; 0, 0, 0, 1; t4, t11, t3, t4 * pkin(5) + t3 * qJ(6) + t77; t2, t9, t1, t2 * pkin(5) + t1 * qJ(6) + t74; t6, t14, t5, t6 * pkin(5) + t5 * qJ(6) + t76; 0, 0, 0, 1;];
T_ges = t7;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
