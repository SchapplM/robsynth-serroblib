% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRPRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6]';
% 
% Output:
% MM_reg [((6+1)*6/2)x35]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 14:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRPRRR12_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR12_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR12_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 00:44:27
% EndTime: 2019-05-07 00:44:31
% DurationCPUTime: 1.21s
% Computational Cost: add. (1194->155), mult. (2711->278), div. (0->0), fcn. (3161->10), ass. (0->107)
t106 = cos(qJ(5));
t77 = sin(qJ(4));
t76 = sin(qJ(5));
t80 = cos(qJ(4));
t98 = t76 * t80;
t50 = t106 * t77 + t98;
t47 = t50 ^ 2;
t86 = t106 * t80;
t49 = t76 * t77 - t86;
t48 = t49 ^ 2;
t119 = -t47 - t48;
t118 = -2 * pkin(2);
t73 = sin(pkin(6));
t81 = cos(qJ(2));
t102 = t73 * t81;
t74 = cos(pkin(6));
t40 = -t77 * t102 + t74 * t80;
t96 = t80 * t102 + t74 * t77;
t23 = t106 * t96 + t76 * t40;
t117 = -0.2e1 * t23;
t116 = -0.2e1 * t49;
t115 = 0.2e1 * t50;
t114 = 0.2e1 * t73;
t113 = 2 * qJ(3);
t82 = -pkin(2) - pkin(9);
t78 = sin(qJ(2));
t112 = pkin(1) * t78;
t111 = pkin(1) * t81;
t62 = t73 * t78;
t56 = pkin(8) * t62;
t88 = -pkin(2) - t111;
t29 = pkin(3) * t62 + t56 + (-pkin(9) + t88) * t74;
t85 = -qJ(3) * t78 - pkin(1);
t33 = (t82 * t81 + t85) * t73;
t15 = t80 * t29 - t77 * t33;
t93 = pkin(4) * t62;
t11 = -t40 * pkin(10) + t15 + t93;
t16 = t77 * t29 + t80 * t33;
t12 = -t96 * pkin(10) + t16;
t6 = t106 * t11 - t76 * t12;
t4 = -pkin(5) * t62 - t6;
t79 = cos(qJ(6));
t110 = t4 * t79;
t109 = t76 * pkin(4);
t89 = t106 * pkin(4);
t67 = -t89 - pkin(5);
t108 = pkin(5) - t67;
t107 = -pkin(10) + t82;
t24 = t106 * t40 - t76 * t96;
t75 = sin(qJ(6));
t19 = t79 * t24 + t75 * t62;
t17 = t19 * t75;
t52 = t107 * t77;
t30 = -t107 * t86 + t76 * t52;
t105 = t30 * t79;
t43 = t49 * t75;
t104 = t49 * t79;
t70 = t73 ^ 2;
t103 = t70 * t81;
t101 = t74 * t81;
t21 = t75 * t23;
t41 = t75 * t50;
t66 = pkin(11) + t109;
t100 = t75 * t66;
t99 = t75 * t79;
t22 = t79 * t23;
t42 = t79 * t50;
t97 = t79 * t66;
t45 = pkin(8) * t102 + t74 * t112;
t64 = t77 * pkin(4) + qJ(3);
t95 = t49 * t115;
t94 = 0.2e1 * t62;
t92 = t50 * t62;
t91 = t77 * t62;
t90 = t82 * t62;
t68 = t74 * qJ(3);
t34 = -t68 - t45;
t87 = t106 * t12;
t32 = pkin(3) * t102 - t34;
t84 = pkin(5) * t49 - pkin(11) * t50;
t83 = -t49 * t67 - t50 * t66;
t7 = t76 * t11 + t87;
t20 = t96 * pkin(4) + t32;
t72 = t79 ^ 2;
t71 = t75 ^ 2;
t61 = t70 * t78 ^ 2;
t60 = 0.2e1 * t99;
t54 = t80 * t62;
t44 = pkin(1) * t101 - t56;
t39 = t49 * t99;
t38 = t49 * t62;
t36 = t88 * t74 + t56;
t35 = (-pkin(2) * t81 + t85) * t73;
t31 = t106 * t52 + t107 * t98;
t28 = t30 * t75;
t27 = (t71 - t72) * t49;
t26 = t50 * pkin(5) + t49 * pkin(11) + t64;
t18 = t75 * t24 - t79 * t62;
t14 = t75 * t26 + t79 * t31;
t13 = t79 * t26 - t75 * t31;
t9 = -t75 * t18 + t19 * t79;
t8 = t23 * pkin(5) - t24 * pkin(11) + t20;
t5 = pkin(11) * t62 + t7;
t3 = t4 * t75;
t2 = t79 * t5 + t75 * t8;
t1 = -t75 * t5 + t79 * t8;
t10 = [1, 0, 0, t61, 0.2e1 * t78 * t103, t74 * t94, t101 * t114, t74 ^ 2, 0.2e1 * pkin(1) * t103 + 0.2e1 * t44 * t74, -0.2e1 * t70 * t112 - 0.2e1 * t45 * t74 (-t34 * t81 + t36 * t78) * t114, 0.2e1 * t35 * t102 + 0.2e1 * t36 * t74, -0.2e1 * t34 * t74 - 0.2e1 * t35 * t62, t34 ^ 2 + t35 ^ 2 + t36 ^ 2, t40 ^ 2, -0.2e1 * t40 * t96, t40 * t94, -t96 * t94, t61, 0.2e1 * t15 * t62 + 0.2e1 * t32 * t96, -0.2e1 * t16 * t62 + 0.2e1 * t32 * t40, t24 ^ 2, t24 * t117, t24 * t94, t62 * t117, t61, 0.2e1 * t20 * t23 + 0.2e1 * t6 * t62, 0.2e1 * t20 * t24 - 0.2e1 * t7 * t62, t19 ^ 2, -0.2e1 * t19 * t18, 0.2e1 * t19 * t23, t18 * t117, t23 ^ 2, 0.2e1 * t1 * t23 + 0.2e1 * t4 * t18, 0.2e1 * t4 * t19 - 0.2e1 * t2 * t23; 0, 0, 0, 0, 0, t62, t102, t74, t44, -t45 (-pkin(2) * t78 + qJ(3) * t81) * t73, t56 + (t118 - t111) * t74, 0.2e1 * t68 + t45, -t36 * pkin(2) - t34 * qJ(3), t40 * t80, -t40 * t77 - t80 * t96, t54, -t91, 0, qJ(3) * t96 + t32 * t77 + t80 * t90, qJ(3) * t40 + t32 * t80 - t77 * t90, -t24 * t49, t49 * t23 - t24 * t50, -t38, -t92, 0, t20 * t50 + t64 * t23 - t30 * t62, -t20 * t49 + t64 * t24 - t31 * t62, -t19 * t104 (t18 * t79 + t17) * t49, t19 * t50 - t49 * t22, -t18 * t50 + t49 * t21, t23 * t50, t1 * t50 + t13 * t23 + t30 * t18 - t4 * t43, -t4 * t104 - t14 * t23 + t30 * t19 - t2 * t50; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, t118, t113, pkin(2) ^ 2 + qJ(3) ^ 2, t80 ^ 2, -0.2e1 * t80 * t77, 0, 0, 0, t77 * t113, t80 * t113, t48, t95, 0, 0, 0, t64 * t115, t64 * t116, t72 * t48, -0.2e1 * t48 * t99, t42 * t116, t75 * t95, t47, 0.2e1 * t13 * t50 - 0.2e1 * t30 * t43, -0.2e1 * t30 * t104 - 0.2e1 * t14 * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, t74, 0, t36, 0, 0, 0, 0, 0, t54, -t91, 0, 0, 0, 0, 0, -t38, -t92, 0, 0, 0, 0, 0, t49 * t18 - t23 * t41, t49 * t19 - t23 * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t119 * t75, t119 * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, -t96, t62, t15, -t16, 0, 0, t24, -t23, t62, t89 * t62 + t6, -t87 + (-t11 - t93) * t76, t17, t9, t21, t22, 0, -t23 * t100 + t67 * t18 - t110, t67 * t19 - t23 * t97 + t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80, -t77, 0, t80 * t82, -t77 * t82, 0, 0, -t49, -t50, 0, -t30, -t31, -t39, t27, t41, t42, 0, t83 * t75 - t105, t83 * t79 + t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80, -t77, 0, 0, 0, 0, 0, -t49, -t50, 0, 0, 0, 0, 0, -t104, t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t89, -0.2e1 * t109, t71, t60, 0, 0, 0, -0.2e1 * t67 * t79, 0.2e1 * t67 * t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, -t23, t62, t6, -t7, t17, t9, t21, t22, 0, -pkin(5) * t18 - pkin(11) * t21 - t110, -pkin(5) * t19 - pkin(11) * t22 + t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49, -t50, 0, -t30, -t31, -t39, t27, t41, t42, 0, t84 * t75 - t105, t84 * t79 + t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49, -t50, 0, 0, 0, 0, 0, -t104, t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t89, -t109, t71, t60, 0, 0, 0, t108 * t79, -t108 * t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t71, t60, 0, 0, 0, 0.2e1 * pkin(5) * t79, -0.2e1 * pkin(5) * t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, -t18, t23, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t104, t43, t50, t13, -t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41, -t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, t79, 0, -t100, -t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, t79, 0, -t75 * pkin(11), -t79 * pkin(11); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t10;
