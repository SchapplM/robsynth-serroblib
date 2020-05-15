% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,theta2]';
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
% Datum: 2018-11-23 16:31
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6RPRRRP12_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP12_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP12_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:31:12
% EndTime: 2018-11-23 16:31:12
% DurationCPUTime: 0.30s
% Computational Cost: add. (1819->84), mult. (1878->101), div. (0->0), fcn. (1949->22), ass. (0->77)
t63 = sin(pkin(6));
t70 = sin(qJ(1));
t102 = t70 * t63;
t60 = pkin(6) - pkin(12);
t53 = cos(t60) / 0.2e1;
t59 = pkin(6) + pkin(12);
t55 = cos(t59);
t45 = t53 + t55 / 0.2e1;
t61 = sin(pkin(12));
t73 = cos(qJ(1));
t38 = -t70 * t45 - t73 * t61;
t62 = sin(pkin(7));
t65 = cos(pkin(7));
t87 = t65 * t102 - t38 * t62;
t52 = sin(t59) / 0.2e1;
t54 = sin(t60);
t43 = t52 + t54 / 0.2e1;
t66 = cos(pkin(6));
t89 = -t43 * t62 + t66 * t65;
t106 = cos(qJ(4));
t101 = t73 * t63;
t100 = qJ(2) * t63;
t99 = pkin(7) - qJ(3);
t98 = pkin(7) + qJ(3);
t97 = pkin(8) + 0;
t95 = t66 * qJ(2) + t97;
t94 = t73 * pkin(1) + t70 * t100 + 0;
t93 = cos(t98);
t92 = sin(t99);
t91 = cos(t99) / 0.2e1;
t90 = sin(t98) / 0.2e1;
t36 = t73 * t45 - t70 * t61;
t88 = -t65 * t101 - t36 * t62;
t86 = t70 * pkin(1) - t73 * t100 + 0;
t46 = t53 - t55 / 0.2e1;
t85 = t46 * pkin(2) + t89 * pkin(9) + t95;
t44 = t52 - t54 / 0.2e1;
t64 = cos(pkin(12));
t39 = -t70 * t44 + t73 * t64;
t84 = t39 * pkin(2) + t87 * pkin(9) + t94;
t83 = t91 + t93 / 0.2e1;
t82 = t90 + t92 / 0.2e1;
t81 = t63 * t82;
t69 = sin(qJ(3));
t26 = -t43 * t83 + t46 * t69 - t66 * t82;
t47 = t90 - t92 / 0.2e1;
t48 = t91 - t93 / 0.2e1;
t72 = cos(qJ(3));
t27 = t43 * t47 + t46 * t72 + t66 * t48;
t80 = t27 * pkin(3) + t26 * pkin(10) + t85;
t22 = -t38 * t83 + t39 * t69 - t70 * t81;
t23 = t48 * t102 + t38 * t47 + t39 * t72;
t79 = t23 * pkin(3) + t22 * pkin(10) + t84;
t37 = t73 * t44 + t70 * t64;
t78 = t37 * pkin(2) + t88 * pkin(9) + t86;
t68 = sin(qJ(4));
t14 = -t89 * t106 + t27 * t68;
t15 = t27 * t106 + t89 * t68;
t77 = t15 * pkin(4) + t14 * pkin(11) + t80;
t11 = -t87 * t106 + t23 * t68;
t12 = t23 * t106 + t87 * t68;
t76 = t12 * pkin(4) + t11 * pkin(11) + t79;
t20 = -t36 * t83 + t37 * t69 + t73 * t81;
t21 = -t48 * t101 + t36 * t47 + t37 * t72;
t75 = t21 * pkin(3) + t20 * pkin(10) + t78;
t10 = t21 * t106 + t88 * t68;
t9 = -t88 * t106 + t21 * t68;
t74 = t10 * pkin(4) + t9 * pkin(11) + t75;
t71 = cos(qJ(5));
t67 = sin(qJ(5));
t6 = t15 * t71 + t26 * t67;
t5 = t15 * t67 - t26 * t71;
t4 = t12 * t71 + t22 * t67;
t3 = t12 * t67 - t22 * t71;
t2 = t10 * t71 + t20 * t67;
t1 = t10 * t67 - t20 * t71;
t7 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t73, -t70, 0, 0; t70, t73, 0, 0; 0, 0, 1, t97; 0, 0, 0, 1; t39, t38, t102, t94; t37, t36, -t101, t86; t46, t43, t66, t95; 0, 0, 0, 1; t23, -t22, t87, t84; t21, -t20, t88, t78; t27, -t26, t89, t85; 0, 0, 0, 1; t12, -t11, t22, t79; t10, -t9, t20, t75; t15, -t14, t26, t80; 0, 0, 0, 1; t4, -t3, t11, t76; t2, -t1, t9, t74; t6, -t5, t14, t77; 0, 0, 0, 1; t4, t11, t3, t4 * pkin(5) + t3 * qJ(6) + t76; t2, t9, t1, t2 * pkin(5) + t1 * qJ(6) + t74; t6, t14, t5, t6 * pkin(5) + t5 * qJ(6) + t77; 0, 0, 0, 1;];
T_ges = t7;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
Tc_stack = NaN(3*size(T_c_mdh,3),4);
% Zusätzliche Ausgabe: Als 2D-array gestapelt, ohne Zeile mit 0001
for i = 1:size(T_c_mdh,3), Tc_stack((i-1)*3+1:3*i,1:4) = T_c_mdh(1:3,1:4,i); end
