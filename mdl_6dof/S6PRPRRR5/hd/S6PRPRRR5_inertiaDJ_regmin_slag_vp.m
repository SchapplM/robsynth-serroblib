% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6PRPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x28]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6PRPRRR5_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR5_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR5_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR5_inertiaDJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:43:35
% EndTime: 2019-03-08 20:43:38
% DurationCPUTime: 0.85s
% Computational Cost: add. (895->134), mult. (2173->251), div. (0->0), fcn. (2086->10), ass. (0->97)
t71 = sin(qJ(4));
t74 = cos(qJ(4));
t121 = cos(qJ(5));
t95 = t121 * qJD(4);
t120 = sin(qJ(5));
t97 = qJD(4) * t120;
t125 = -t71 * t95 - t74 * t97;
t94 = t121 * qJD(5);
t124 = t95 + t94;
t73 = cos(qJ(6));
t67 = t73 ^ 2;
t70 = sin(qJ(6));
t109 = t70 ^ 2 - t67;
t93 = qJD(6) * t109;
t42 = t120 * t71 - t121 * t74;
t40 = t42 ^ 2;
t123 = 2 * qJD(3);
t76 = -pkin(2) - pkin(8);
t122 = pkin(9) - t76;
t96 = qJD(5) * t120;
t28 = -t71 * t94 - t74 * t96 + t125;
t119 = t42 * t28;
t118 = t42 * t70;
t117 = t42 * t73;
t89 = t71 * t97;
t29 = t124 * t74 - t71 * t96 - t89;
t43 = t120 * t74 + t121 * t71;
t116 = t43 * t29;
t68 = sin(pkin(6));
t72 = sin(qJ(2));
t115 = t68 * t72;
t75 = cos(qJ(2));
t114 = t68 * t75;
t113 = t70 * t28;
t112 = t73 * t28;
t45 = t122 * t71;
t46 = t122 * t74;
t32 = -t120 * t46 - t121 * t45;
t14 = t32 * qJD(5) + t125 * t122;
t31 = -t120 * t45 + t121 * t46;
t64 = qJD(6) * t73;
t111 = t14 * t70 + t31 * t64;
t62 = -t121 * pkin(4) - pkin(5);
t91 = pkin(4) * t96;
t110 = t62 * t64 + t70 * t91;
t59 = t71 * pkin(4) + qJ(3);
t108 = qJD(2) * t75;
t107 = qJD(6) * t70;
t106 = t71 * qJD(4);
t105 = t74 * qJD(4);
t51 = pkin(4) * t105 + qJD(3);
t104 = qJ(3) * qJD(4);
t103 = pkin(5) * t107;
t102 = pkin(5) * t64;
t52 = qJD(2) * t115;
t101 = t68 * t108;
t100 = t70 * t64;
t99 = t43 ^ 2 + t40;
t98 = 0.4e1 * t70 * t117;
t92 = pkin(4) * t94;
t21 = t43 * pkin(5) + t42 * pkin(10) + t59;
t87 = t73 * t21 - t70 * t32;
t86 = t70 * t21 + t73 * t32;
t61 = t120 * pkin(4) + pkin(10);
t85 = t42 * t62 + t43 * t61;
t69 = cos(pkin(6));
t81 = t71 * t114 - t69 * t74;
t82 = t74 * t114 + t69 * t71;
t20 = -t120 * t82 - t121 * t81;
t84 = t73 * t115 - t70 * t20;
t83 = t70 * t115 + t73 * t20;
t17 = t42 * t64 - t113;
t18 = t42 * t107 + t112;
t16 = t70 * t29 + t43 * t64;
t15 = t43 * t107 - t73 * t29;
t80 = t62 * t107 - t73 * t91;
t79 = -0.2e1 * t116 + 0.2e1 * t119;
t78 = t81 * qJD(4) + t74 * t52;
t77 = t28 * t62 - t29 * t61 + (-t120 * t42 - t121 * t43) * qJD(5) * pkin(4);
t50 = 0.2e1 * t100;
t41 = -0.2e1 * t93;
t30 = t82 * qJD(4) - t71 * t52;
t22 = t31 * t107;
t19 = -t120 * t81 + t121 * t82;
t13 = -t122 * t89 + t124 * t46 - t45 * t96;
t11 = t29 * pkin(5) - t28 * pkin(10) + t51;
t10 = t70 * t112 + t42 * t93;
t9 = qJD(6) * t98 - t109 * t28;
t8 = t20 * qJD(5) - t120 * t30 - t121 * t78;
t7 = -t120 * t78 + t121 * t30 - t81 * t96 + t82 * t94;
t6 = t19 * t107 - t8 * t73;
t5 = t19 * t64 + t8 * t70;
t4 = -t83 * qJD(6) + t73 * t101 + t70 * t7;
t3 = -t84 * qJD(6) - t70 * t101 + t73 * t7;
t2 = -t86 * qJD(6) + t73 * t11 + t70 * t13;
t1 = -t87 * qJD(6) - t70 * t11 + t73 * t13;
t12 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t52, -t101, t52, t101 (qJD(3) * t72 + (-pkin(2) * t72 + qJ(3) * t75) * qJD(2)) * t68, 0, 0, 0, 0, 0 (t72 * t105 + t71 * t108) * t68 (-t72 * t106 + t74 * t108) * t68, 0, 0, 0, 0, 0 (t43 * t108 + t29 * t72) * t68 (-t42 * t108 + t28 * t72) * t68, 0, 0, 0, 0, 0, -t8 * t118 - t17 * t19 + t84 * t29 + t4 * t43, -t8 * t117 + t18 * t19 - t83 * t29 + t3 * t43; 0, 0, 0, 0, 0, t123, qJ(3) * t123, -0.2e1 * t71 * t105, 0.2e1 * (t71 ^ 2 - t74 ^ 2) * qJD(4), 0, 0, 0, 0.2e1 * qJD(3) * t71 + 0.2e1 * t74 * t104, 0.2e1 * qJD(3) * t74 - 0.2e1 * t71 * t104, -0.2e1 * t119, -0.2e1 * t28 * t43 + 0.2e1 * t42 * t29, 0, 0, 0, 0.2e1 * t59 * t29 + 0.2e1 * t51 * t43, 0.2e1 * t59 * t28 - 0.2e1 * t51 * t42, -0.2e1 * t40 * t100 - 0.2e1 * t67 * t119, t28 * t98 + 0.2e1 * t40 * t93, 0.2e1 * t43 * t112 + 0.2e1 * t15 * t42, -0.2e1 * t43 * t113 + 0.2e1 * t16 * t42, 0.2e1 * t116, -0.2e1 * t14 * t118 - 0.2e1 * t17 * t31 + 0.2e1 * t2 * t43 + 0.2e1 * t87 * t29, 0.2e1 * t1 * t43 - 0.2e1 * t14 * t117 + 0.2e1 * t18 * t31 - 0.2e1 * t86 * t29; 0, 0, 0, 0, 0, 0, t52, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t99 * t64 + t70 * t79, t99 * t107 + t73 * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t78, t30, 0, 0, 0, 0, 0, -t8, t7, 0, 0, 0, 0, 0, t6, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t106, -t105, 0, -t76 * t106, -t76 * t105, 0, 0, t28, -t29, 0, -t14, t13, t10, t9, t16, -t15, 0, t22 + (-t85 * qJD(6) - t14) * t73 + t77 * t70, t85 * t107 + t77 * t73 + t111; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t106, -t105, 0, 0, 0, 0, 0, t28, -t29, 0, 0, 0, 0, 0, t18, t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t91, -0.2e1 * t92, t50, t41, 0, 0, 0, 0.2e1 * t80, 0.2e1 * t110; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, t7, 0, 0, 0, 0, 0, t6, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, -t29, 0, -t14, t13, t10, t9, t16, -t15, 0, t22 + (-pkin(5) * t28 - pkin(10) * t29) * t70 + (-t14 + (pkin(5) * t42 - pkin(10) * t43) * qJD(6)) * t73, -t18 * pkin(5) + t15 * pkin(10) + t111; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, -t29, 0, 0, 0, 0, 0, t18, t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t91, -t92, t50, t41, 0, 0, 0, t80 - t103, -t102 + t110; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, t41, 0, 0, 0, -0.2e1 * t103, -0.2e1 * t102; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, t17, t29, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64, -t107, 0, -t61 * t64 - t70 * t92, t61 * t107 - t73 * t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64, -t107, 0, -pkin(10) * t64, pkin(10) * t107; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t12;
