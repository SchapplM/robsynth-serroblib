% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5]';
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
% Datum: 2018-11-23 18:36
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RRRRRP12_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP12_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP12_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:35:55
% EndTime: 2018-11-23 18:35:55
% DurationCPUTime: 0.29s
% Computational Cost: add. (1819->84), mult. (1878->100), div. (0->0), fcn. (1949->22), ass. (0->76)
t62 = sin(pkin(6));
t69 = sin(qJ(1));
t101 = t69 * t62;
t60 = pkin(6) - qJ(2);
t53 = cos(t60) / 0.2e1;
t59 = pkin(6) + qJ(2);
t55 = cos(t59);
t47 = t53 + t55 / 0.2e1;
t68 = sin(qJ(2));
t73 = cos(qJ(1));
t38 = -t69 * t47 - t73 * t68;
t61 = sin(pkin(7));
t63 = cos(pkin(7));
t87 = t63 * t101 - t38 * t61;
t52 = sin(t59) / 0.2e1;
t54 = sin(t60);
t44 = t52 + t54 / 0.2e1;
t64 = cos(pkin(6));
t89 = -t44 * t61 + t64 * t63;
t105 = cos(qJ(4));
t100 = t73 * t62;
t99 = pkin(7) - qJ(3);
t98 = pkin(7) + qJ(3);
t97 = pkin(8) + 0;
t95 = t64 * pkin(9) + t97;
t94 = t73 * pkin(1) + pkin(9) * t101 + 0;
t93 = cos(t98);
t92 = sin(t99);
t91 = cos(t99) / 0.2e1;
t90 = sin(t98) / 0.2e1;
t36 = t73 * t47 - t69 * t68;
t88 = -t63 * t100 - t36 * t61;
t86 = t69 * pkin(1) - pkin(9) * t100 + 0;
t48 = t53 - t55 / 0.2e1;
t85 = t48 * pkin(2) + t89 * pkin(10) + t95;
t45 = t52 - t54 / 0.2e1;
t72 = cos(qJ(2));
t39 = -t69 * t45 + t73 * t72;
t84 = t39 * pkin(2) + t87 * pkin(10) + t94;
t83 = t91 + t93 / 0.2e1;
t82 = t90 + t92 / 0.2e1;
t81 = t62 * t82;
t67 = sin(qJ(3));
t26 = -t44 * t83 + t48 * t67 - t64 * t82;
t43 = t90 - t92 / 0.2e1;
t46 = t91 - t93 / 0.2e1;
t71 = cos(qJ(3));
t27 = t44 * t43 + t64 * t46 + t48 * t71;
t80 = t27 * pkin(3) + t26 * pkin(11) + t85;
t22 = -t38 * t83 + t39 * t67 - t69 * t81;
t23 = t46 * t101 + t38 * t43 + t39 * t71;
t79 = t23 * pkin(3) + t22 * pkin(11) + t84;
t37 = t73 * t45 + t69 * t72;
t78 = t37 * pkin(2) + t88 * pkin(10) + t86;
t66 = sin(qJ(4));
t14 = -t89 * t105 + t27 * t66;
t15 = t27 * t105 + t89 * t66;
t77 = t15 * pkin(4) + t14 * pkin(12) + t80;
t11 = -t87 * t105 + t23 * t66;
t12 = t23 * t105 + t87 * t66;
t76 = t12 * pkin(4) + t11 * pkin(12) + t79;
t20 = -t36 * t83 + t37 * t67 + t73 * t81;
t21 = -t46 * t100 + t36 * t43 + t37 * t71;
t75 = t21 * pkin(3) + t20 * pkin(11) + t78;
t10 = t21 * t105 + t88 * t66;
t9 = -t88 * t105 + t21 * t66;
t74 = t10 * pkin(4) + t9 * pkin(12) + t75;
t70 = cos(qJ(5));
t65 = sin(qJ(5));
t6 = t15 * t70 + t26 * t65;
t5 = t15 * t65 - t26 * t70;
t4 = t12 * t70 + t22 * t65;
t3 = t12 * t65 - t22 * t70;
t2 = t10 * t70 + t20 * t65;
t1 = t10 * t65 - t20 * t70;
t7 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t73, -t69, 0, 0; t69, t73, 0, 0; 0, 0, 1, t97; 0, 0, 0, 1; t39, t38, t101, t94; t37, t36, -t100, t86; t48, t44, t64, t95; 0, 0, 0, 1; t23, -t22, t87, t84; t21, -t20, t88, t78; t27, -t26, t89, t85; 0, 0, 0, 1; t12, -t11, t22, t79; t10, -t9, t20, t75; t15, -t14, t26, t80; 0, 0, 0, 1; t4, -t3, t11, t76; t2, -t1, t9, t74; t6, -t5, t14, t77; 0, 0, 0, 1; t4, t11, t3, t4 * pkin(5) + t3 * qJ(6) + t76; t2, t9, t1, t2 * pkin(5) + t1 * qJ(6) + t74; t6, t14, t5, t6 * pkin(5) + t5 * qJ(6) + t77; 0, 0, 0, 1;];
T_ges = t7;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
