% Calculate inertial parameters regressor of joint inertia matrix for
% S6RRPPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta4]';
% 
% Output:
% MM_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRPPRP5_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP5_inertiaJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP5_inertiaJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t72 = sin(qJ(2));
t67 = t72 ^ 2;
t73 = cos(qJ(2));
t68 = t73 ^ 2;
t119 = t67 + t68;
t115 = -0.2e1 * t72;
t110 = cos(qJ(5));
t69 = sin(pkin(9));
t109 = sin(qJ(5));
t70 = cos(pkin(9));
t87 = t109 * t70;
t40 = t110 * t69 + t87;
t37 = t40 ^ 2;
t89 = t110 * t70;
t42 = -t109 * t69 + t89;
t38 = t42 ^ 2;
t118 = t37 + t38;
t71 = -pkin(2) - qJ(4);
t111 = -pkin(8) + t71;
t46 = t111 * t69;
t22 = t109 * t46 - t111 * t89;
t24 = t110 * t46 + t111 * t87;
t84 = -t22 * t42 + t24 * t40;
t27 = t42 * t73;
t117 = t27 ^ 2;
t55 = t69 * pkin(4) + qJ(3);
t116 = 0.2e1 * t55;
t114 = 0.2e1 * t73;
t113 = 0.2e1 * qJ(3);
t112 = t72 * pkin(5);
t85 = -t72 * qJ(3) - pkin(1);
t36 = t71 * t73 + t85;
t60 = t72 * pkin(7);
t49 = t72 * pkin(3) + t60;
t45 = t70 * t49;
t13 = t72 * pkin(4) + t45 + (pkin(8) * t73 - t36) * t69;
t20 = t70 * t36 + t69 * t49;
t97 = t70 * t73;
t15 = -pkin(8) * t97 + t20;
t4 = t109 * t13 + t110 * t15;
t107 = t22 * t72;
t105 = t24 * t72;
t104 = t27 * t72;
t29 = t40 * t73;
t103 = t29 * t27;
t18 = t40 * t27;
t102 = t42 * t29;
t101 = t42 * t40;
t33 = t42 * t72;
t100 = t69 * t72;
t99 = t69 * t73;
t98 = t70 * t69;
t32 = t72 * t40;
t96 = t72 * t73;
t95 = t119 * pkin(7) ^ 2;
t62 = t73 * pkin(7);
t50 = t73 * pkin(3) + t62;
t64 = t69 ^ 2;
t65 = t70 ^ 2;
t52 = t64 + t65;
t94 = t72 * qJ(6);
t93 = t73 * qJ(3);
t92 = -0.2e1 * t96;
t91 = t22 ^ 2 + t24 ^ 2;
t90 = t69 * t97;
t34 = pkin(4) * t97 + t50;
t86 = -t22 * t29 - t24 * t27;
t83 = -t18 + t102;
t1 = t94 + t4;
t3 = -t109 * t15 + t110 * t13;
t2 = -t3 - t112;
t82 = t1 * t40 - t2 * t42;
t81 = t3 * t42 + t4 * t40;
t80 = -t72 * pkin(2) + t93;
t79 = pkin(5) * t42 + t40 * qJ(6);
t19 = -t69 * t36 + t45;
t5 = t19 * t70 + t20 * t69;
t78 = -t42 * t27 + t29 * t40;
t77 = t71 * t72 + t93;
t76 = 0.2e1 * t84;
t74 = qJ(3) ^ 2;
t57 = t70 * t72;
t54 = 0.2e1 * t96;
t48 = -t73 * pkin(2) + t85;
t47 = 0.2e1 * t119 * pkin(7);
t39 = t52 * t71;
t26 = t29 ^ 2;
t25 = t29 * t115;
t16 = t40 * pkin(5) - t42 * qJ(6) + t55;
t6 = t27 * pkin(5) + t29 * qJ(6) + t34;
t7 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t67, t54, 0, t68, 0, 0, pkin(1) * t114, pkin(1) * t115, t47, pkin(1) ^ 2 + t95, 0, 0, 0, t67, t54, t68, t47, t48 * t114, t48 * t115, t48 ^ 2 + t95, t64 * t68, 0.2e1 * t68 * t98, t69 * t92, t65 * t68, t70 * t92, t67, 0.2e1 * t19 * t72 + 0.2e1 * t50 * t97, -0.2e1 * t20 * t72 - 0.2e1 * t50 * t99 (t19 * t69 - t20 * t70) * t114, t19 ^ 2 + t20 ^ 2 + t50 ^ 2, t26, 0.2e1 * t103, t25, t117, -0.2e1 * t104, t67, 0.2e1 * t34 * t27 + 0.2e1 * t3 * t72, -0.2e1 * t34 * t29 - 0.2e1 * t4 * t72, -0.2e1 * t4 * t27 + 0.2e1 * t3 * t29, t3 ^ 2 + t34 ^ 2 + t4 ^ 2, t26, t25, -0.2e1 * t103, t67, 0.2e1 * t104, t117, -0.2e1 * t2 * t72 + 0.2e1 * t6 * t27, -0.2e1 * t1 * t27 - 0.2e1 * t2 * t29, 0.2e1 * t1 * t72 + 0.2e1 * t6 * t29, t1 ^ 2 + t2 ^ 2 + t6 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72, 0, t73, 0, -t60, -t62, 0, 0, 0, -t72, -t73, 0, 0, 0, t80, t60, t62, t80 * pkin(7), -t90 (t64 - t65) * t73, t57, t90, -t100, 0, t50 * t69 + t77 * t70, t50 * t70 - t69 * t77, -t5, t50 * qJ(3) + t5 * t71, -t102, t78, t33, t18, -t32, 0, t55 * t27 + t34 * t40 - t107, -t55 * t29 + t34 * t42 - t105, -t81 + t86, -t3 * t22 + t4 * t24 + t34 * t55, -t102, t33, -t78, 0, t32, t18, t16 * t27 + t6 * t40 - t107, -t82 + t86, t16 * t29 - t6 * t42 + t105, t1 * t24 + t6 * t16 + t2 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -0.2e1 * pkin(2), t113, pkin(2) ^ 2 + t74, t65, -0.2e1 * t98, 0, t64, 0, 0, t69 * t113, t70 * t113, -0.2e1 * t39, t52 * t71 ^ 2 + t74, t38, -0.2e1 * t101, 0, t37, 0, 0, t40 * t116, t42 * t116, -t76, t55 ^ 2 + t91, t38, 0, 0.2e1 * t101, 0, 0, t37, 0.2e1 * t16 * t40, -t76, -0.2e1 * t16 * t42, t16 ^ 2 + t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72, 0, 0, t60, 0, 0, 0, 0, 0, 0, t57, -t100, 0, t5, 0, 0, 0, 0, 0, 0, t33, -t32, t83, t81, 0, 0, 0, 0, 0, 0, t33, t83, t32, t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(2), 0, 0, 0, 0, 0, 0, 0, 0, -t52, t39, 0, 0, 0, 0, 0, 0, 0, 0, -t118, t84, 0, 0, 0, 0, 0, 0, 0, -t118, 0, t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, 0, 0, 0, 0, 0, 0, 0, 0, 0, t118, 0, 0, 0, 0, 0, 0, 0, 0, 0, t118; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t97, -t99, 0, t50, 0, 0, 0, 0, 0, 0, t27, -t29, 0, t34, 0, 0, 0, 0, 0, 0, t27, 0, t29, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, t70, 0, qJ(3), 0, 0, 0, 0, 0, 0, t40, t42, 0, t55, 0, 0, 0, 0, 0, 0, t40, 0, -t42, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29, 0, -t27, t72, t3, -t4, 0, 0, 0, -t29, 0, t72, t27, 0, t3 + 0.2e1 * t112, pkin(5) * t29 - t27 * qJ(6), 0.2e1 * t94 + t4, -t2 * pkin(5) + t1 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, 0, -t40, 0, -t22, -t24, 0, 0, 0, t42, 0, 0, t40, 0, -t22, -t79, t24, -t22 * pkin(5) + t24 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, -t40, 0, 0, 0, 0, 0, 0, 0, 0, t42, 0, t40, t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(5), 0, 0.2e1 * qJ(6), pkin(5) ^ 2 + qJ(6) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t72, -t29, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, 0, t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t7;
