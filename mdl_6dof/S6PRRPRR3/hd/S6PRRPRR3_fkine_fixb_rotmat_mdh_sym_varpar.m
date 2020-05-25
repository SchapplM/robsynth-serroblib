% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
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
% Datum: 2018-11-23 15:16
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6PRRPRR3_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR3_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR3_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:15:42
% EndTime: 2018-11-23 15:15:42
% DurationCPUTime: 0.30s
% Computational Cost: add. (1517->108), mult. (1294->132), div. (0->0), fcn. (1306->28), ass. (0->85)
t69 = pkin(7) + qJ(3);
t52 = sin(t69) / 0.2e1;
t70 = pkin(7) - qJ(3);
t59 = sin(t70);
t44 = t52 - t59 / 0.2e1;
t54 = cos(t70) / 0.2e1;
t62 = cos(t69);
t48 = t54 - t62 / 0.2e1;
t73 = sin(pkin(12));
t75 = sin(pkin(6));
t109 = t73 * t75;
t76 = cos(pkin(12));
t108 = t76 * t75;
t68 = qJ(3) + pkin(13);
t107 = t73 * pkin(1) + 0;
t106 = qJ(1) + 0;
t105 = t76 * pkin(1) + pkin(8) * t109 + 0;
t78 = cos(pkin(6));
t104 = t78 * pkin(8) + t106;
t103 = pkin(7) - t68;
t102 = pkin(7) + t68;
t71 = pkin(6) + qJ(2);
t53 = sin(t71) / 0.2e1;
t72 = pkin(6) - qJ(2);
t60 = sin(t72);
t45 = t53 + t60 / 0.2e1;
t74 = sin(pkin(7));
t77 = cos(pkin(7));
t31 = -t45 * t74 + t78 * t77;
t101 = cos(t102);
t100 = sin(t103);
t55 = cos(t72) / 0.2e1;
t64 = cos(t71);
t49 = t55 + t64 / 0.2e1;
t83 = sin(qJ(2));
t32 = t76 * t49 - t73 * t83;
t22 = -t77 * t108 - t32 * t74;
t34 = -t73 * t49 - t76 * t83;
t23 = t77 * t109 - t34 * t74;
t99 = -pkin(8) * t108 + t107;
t98 = cos(t103) / 0.2e1;
t97 = sin(t102) / 0.2e1;
t46 = t53 - t60 / 0.2e1;
t87 = cos(qJ(2));
t35 = -t73 * t46 + t76 * t87;
t79 = pkin(9) + qJ(4);
t37 = t44 * pkin(3) - t74 * t79;
t38 = t48 * pkin(3) + t77 * t79;
t86 = cos(qJ(3));
t56 = t86 * pkin(3) + pkin(2);
t96 = t38 * t109 + t34 * t37 + t35 * t56 + t105;
t50 = t55 - t64 / 0.2e1;
t95 = t45 * t37 + t78 * t38 + t50 * t56 + t104;
t33 = t76 * t46 + t73 * t87;
t94 = t32 * t37 + t33 * t56 + (-pkin(8) - t38) * t108 + t107;
t57 = sin(t68);
t89 = t100 / 0.2e1 + t97;
t88 = t75 * t89;
t90 = t101 / 0.2e1 + t98;
t11 = -t34 * t90 + t35 * t57 - t73 * t88;
t41 = t97 - t100 / 0.2e1;
t42 = t98 - t101 / 0.2e1;
t61 = cos(t68);
t12 = t42 * t109 + t34 * t41 + t35 * t61;
t93 = t12 * pkin(4) + t11 * pkin(10) + t96;
t14 = -t45 * t90 + t50 * t57 - t78 * t89;
t15 = t45 * t41 + t78 * t42 + t50 * t61;
t92 = t15 * pkin(4) + t14 * pkin(10) + t95;
t10 = -t42 * t108 + t32 * t41 + t33 * t61;
t9 = -t32 * t90 + t33 * t57 + t76 * t88;
t91 = t10 * pkin(4) + t9 * pkin(10) + t94;
t85 = cos(qJ(5));
t84 = cos(qJ(6));
t82 = sin(qJ(3));
t81 = sin(qJ(5));
t80 = sin(qJ(6));
t47 = t54 + t62 / 0.2e1;
t43 = t52 + t59 / 0.2e1;
t6 = t15 * t85 + t31 * t81;
t5 = t15 * t81 - t31 * t85;
t4 = t12 * t85 + t23 * t81;
t3 = t12 * t81 - t23 * t85;
t2 = t10 * t85 + t22 * t81;
t1 = t10 * t81 - t22 * t85;
t7 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t76, -t73, 0, 0; t73, t76, 0, 0; 0, 0, 1, t106; 0, 0, 0, 1; t35, t34, t109, t105; t33, t32, -t108, t99; t50, t45, t78, t104; 0, 0, 0, 1; t48 * t109 + t34 * t44 + t35 * t86, t43 * t109 + t34 * t47 - t35 * t82, t23, t35 * pkin(2) + t23 * pkin(9) + t105; -t48 * t108 + t32 * t44 + t33 * t86, -t43 * t108 + t32 * t47 - t33 * t82, t22, t33 * pkin(2) + t22 * pkin(9) + t99; t45 * t44 + t78 * t48 + t50 * t86, t78 * t43 + t45 * t47 - t50 * t82, t31, t50 * pkin(2) + t31 * pkin(9) + t104; 0, 0, 0, 1; t12, -t11, t23, t96; t10, -t9, t22, t94; t15, -t14, t31, t95; 0, 0, 0, 1; t4, -t3, t11, t93; t2, -t1, t9, t91; t6, -t5, t14, t92; 0, 0, 0, 1; t11 * t80 + t4 * t84, t11 * t84 - t4 * t80, t3, t4 * pkin(5) + t3 * pkin(11) + t93; t2 * t84 + t9 * t80, -t2 * t80 + t9 * t84, t1, t2 * pkin(5) + t1 * pkin(11) + t91; t14 * t80 + t6 * t84, t14 * t84 - t6 * t80, t5, t6 * pkin(5) + t5 * pkin(11) + t92; 0, 0, 0, 1;];
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
