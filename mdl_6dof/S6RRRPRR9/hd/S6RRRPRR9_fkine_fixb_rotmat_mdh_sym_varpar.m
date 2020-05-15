% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2018-11-23 17:58
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6RRRPRR9_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR9_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR9_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:56:54
% EndTime: 2018-11-23 17:56:55
% DurationCPUTime: 0.29s
% Computational Cost: add. (1517->108), mult. (1294->132), div. (0->0), fcn. (1306->28), ass. (0->85)
t69 = pkin(7) + qJ(3);
t52 = sin(t69) / 0.2e1;
t70 = pkin(7) - qJ(3);
t59 = sin(t70);
t44 = t52 - t59 / 0.2e1;
t54 = cos(t70) / 0.2e1;
t62 = cos(t69);
t48 = t54 - t62 / 0.2e1;
t74 = sin(pkin(6));
t82 = sin(qJ(1));
t109 = t82 * t74;
t87 = cos(qJ(1));
t108 = t87 * t74;
t107 = pkin(8) + 0;
t68 = qJ(3) + pkin(13);
t106 = t82 * pkin(1) + 0;
t76 = cos(pkin(6));
t105 = t76 * pkin(9) + t107;
t104 = t87 * pkin(1) + pkin(9) * t109 + 0;
t103 = pkin(7) - t68;
t102 = pkin(7) + t68;
t71 = pkin(6) + qJ(2);
t53 = sin(t71) / 0.2e1;
t72 = pkin(6) - qJ(2);
t60 = sin(t72);
t45 = t53 + t60 / 0.2e1;
t73 = sin(pkin(7));
t75 = cos(pkin(7));
t31 = -t45 * t73 + t76 * t75;
t101 = cos(t102);
t100 = sin(t103);
t55 = cos(t72) / 0.2e1;
t64 = cos(t71);
t49 = t55 + t64 / 0.2e1;
t81 = sin(qJ(2));
t32 = t87 * t49 - t82 * t81;
t22 = -t75 * t108 - t32 * t73;
t34 = -t82 * t49 - t87 * t81;
t23 = t75 * t109 - t34 * t73;
t99 = -pkin(9) * t108 + t106;
t77 = pkin(10) + qJ(4);
t36 = t44 * pkin(3) - t73 * t77;
t37 = t48 * pkin(3) + t75 * t77;
t50 = t55 - t64 / 0.2e1;
t85 = cos(qJ(3));
t56 = t85 * pkin(3) + pkin(2);
t98 = t45 * t36 + t76 * t37 + t50 * t56 + t105;
t97 = cos(t103) / 0.2e1;
t96 = sin(t102) / 0.2e1;
t46 = t53 - t60 / 0.2e1;
t86 = cos(qJ(2));
t35 = -t82 * t46 + t87 * t86;
t95 = t37 * t109 + t34 * t36 + t35 * t56 + t104;
t33 = t87 * t46 + t82 * t86;
t94 = t32 * t36 + t33 * t56 + (-pkin(9) - t37) * t108 + t106;
t57 = sin(t68);
t89 = t100 / 0.2e1 + t96;
t90 = t101 / 0.2e1 + t97;
t14 = -t45 * t90 + t50 * t57 - t76 * t89;
t41 = t96 - t100 / 0.2e1;
t42 = t97 - t101 / 0.2e1;
t61 = cos(t68);
t15 = t45 * t41 + t76 * t42 + t50 * t61;
t93 = t15 * pkin(4) + t14 * pkin(11) + t98;
t88 = t74 * t89;
t11 = -t34 * t90 + t35 * t57 - t82 * t88;
t12 = t42 * t109 + t34 * t41 + t35 * t61;
t92 = t12 * pkin(4) + t11 * pkin(11) + t95;
t10 = -t42 * t108 + t32 * t41 + t33 * t61;
t9 = -t32 * t90 + t33 * t57 + t87 * t88;
t91 = t10 * pkin(4) + t9 * pkin(11) + t94;
t84 = cos(qJ(5));
t83 = cos(qJ(6));
t80 = sin(qJ(3));
t79 = sin(qJ(5));
t78 = sin(qJ(6));
t47 = t54 + t62 / 0.2e1;
t43 = t52 + t59 / 0.2e1;
t6 = t15 * t84 + t31 * t79;
t5 = t15 * t79 - t31 * t84;
t4 = t12 * t84 + t23 * t79;
t3 = t12 * t79 - t23 * t84;
t2 = t10 * t84 + t22 * t79;
t1 = t10 * t79 - t22 * t84;
t7 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t87, -t82, 0, 0; t82, t87, 0, 0; 0, 0, 1, t107; 0, 0, 0, 1; t35, t34, t109, t104; t33, t32, -t108, t99; t50, t45, t76, t105; 0, 0, 0, 1; t48 * t109 + t34 * t44 + t35 * t85, t43 * t109 + t34 * t47 - t35 * t80, t23, t35 * pkin(2) + t23 * pkin(10) + t104; -t48 * t108 + t32 * t44 + t33 * t85, -t43 * t108 + t32 * t47 - t33 * t80, t22, t33 * pkin(2) + t22 * pkin(10) + t99; t45 * t44 + t76 * t48 + t50 * t85, t76 * t43 + t45 * t47 - t50 * t80, t31, t50 * pkin(2) + t31 * pkin(10) + t105; 0, 0, 0, 1; t12, -t11, t23, t95; t10, -t9, t22, t94; t15, -t14, t31, t98; 0, 0, 0, 1; t4, -t3, t11, t92; t2, -t1, t9, t91; t6, -t5, t14, t93; 0, 0, 0, 1; t11 * t78 + t4 * t83, t11 * t83 - t4 * t78, t3, t4 * pkin(5) + t3 * pkin(12) + t92; t2 * t83 + t9 * t78, -t2 * t78 + t9 * t83, t1, t2 * pkin(5) + t1 * pkin(12) + t91; t14 * t78 + t6 * t83, t14 * t83 - t6 * t78, t5, t6 * pkin(5) + t5 * pkin(12) + t93; 0, 0, 0, 1;];
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
