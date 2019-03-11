% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x32]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPRRPR4_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR4_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR4_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR4_inertiaDJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:10:16
% EndTime: 2019-03-09 05:10:21
% DurationCPUTime: 1.51s
% Computational Cost: add. (3926->169), mult. (8818->302), div. (0->0), fcn. (9262->10), ass. (0->120)
t95 = sin(pkin(11));
t97 = cos(pkin(11));
t137 = t95 ^ 2 + t97 ^ 2;
t158 = t137 * qJD(5);
t160 = 0.2e1 * t158;
t154 = cos(qJ(4));
t128 = qJD(4) * t154;
t124 = pkin(3) * t128;
t79 = t124 + qJD(5);
t159 = t137 * t79;
t101 = sin(qJ(3));
t103 = cos(qJ(3));
t96 = sin(pkin(10));
t98 = cos(pkin(10));
t71 = t101 * t98 + t103 * t96;
t102 = cos(qJ(6));
t99 = sin(qJ(6));
t146 = t99 * t97;
t70 = t102 * t95 + t146;
t65 = t70 * qJD(6);
t110 = t101 * t96 - t103 * t98;
t145 = pkin(7) + qJ(2);
t76 = t145 * t96;
t78 = t145 * t98;
t112 = t101 * t78 + t103 * t76;
t157 = t110 * qJD(2) + t112 * qJD(3);
t86 = -t98 * pkin(2) - pkin(1);
t156 = 0.2e1 * t86;
t155 = t97 * pkin(5);
t90 = t97 * pkin(9);
t100 = sin(qJ(4));
t139 = t71 * qJD(3);
t105 = -t139 * pkin(8) - t157;
t111 = t101 * t76 - t103 * t78;
t107 = -t71 * qJD(2) + t111 * qJD(3);
t66 = t110 * qJD(3);
t106 = t66 * pkin(8) + t107;
t104 = t100 * t105 - t154 * t106;
t51 = -t71 * pkin(8) - t112;
t52 = -t110 * pkin(8) - t111;
t36 = t100 * t51 + t154 * t52;
t20 = t36 * qJD(4) + t104;
t153 = t20 * t97;
t55 = -t100 * t110 + t154 * t71;
t152 = t55 * t95;
t147 = t99 * t95;
t131 = qJD(6) * t147;
t134 = qJD(6) * t97;
t64 = -t102 * t134 + t131;
t85 = -pkin(4) - t155;
t151 = t85 * t64;
t150 = t85 * t65;
t54 = t100 * t71 + t154 * t110;
t39 = -t54 * qJD(4) - t100 * t139 - t154 * t66;
t149 = t95 * t39;
t148 = t97 * t39;
t113 = t102 * t97 - t147;
t133 = qJD(4) * t100;
t12 = pkin(5) * t149 + t52 * t128 + t51 * t133 + t104;
t35 = t100 * t52 - t154 * t51;
t30 = pkin(5) * t152 + t35;
t144 = -t113 * t12 + t30 * t65;
t143 = t12 * t70 - t30 * t64;
t19 = -t100 * t106 - t154 * t105 - t51 * t128 + t52 * t133;
t130 = t139 * pkin(3);
t40 = t55 * qJD(4) - t100 * t66 + t154 * t139;
t23 = t40 * pkin(4) - t39 * qJ(5) - t55 * qJD(5) + t130;
t8 = -t97 * t19 + t95 * t23;
t59 = t110 * pkin(3) + t86;
t34 = t54 * pkin(4) - t55 * qJ(5) + t59;
t25 = t95 * t34 + t97 * t36;
t132 = pkin(3) * t133;
t87 = -t154 * pkin(3) - pkin(4);
t74 = t87 - t155;
t142 = -t113 * t132 + t74 * t65;
t141 = t70 * t132 - t74 * t64;
t7 = t95 * t19 + t97 * t23;
t24 = t97 * t34 - t95 * t36;
t126 = t95 * t132;
t125 = t97 * t132;
t123 = 0.2e1 * (t96 ^ 2 + t98 ^ 2) * qJD(2);
t122 = t7 * t97 + t8 * t95;
t3 = -t7 * t95 + t8 * t97;
t121 = t20 * t55 + t35 * t39;
t120 = -t24 * t95 + t25 * t97;
t28 = t70 * t40 - t64 * t54;
t13 = t54 * pkin(5) - t55 * t90 + t24;
t18 = -pkin(9) * t152 + t25;
t119 = t102 * t18 + t99 * t13;
t118 = t102 * t13 - t99 * t18;
t84 = t100 * pkin(3) + qJ(5);
t67 = (-pkin(9) - t84) * t95;
t68 = t97 * t84 + t90;
t117 = t102 * t68 + t99 * t67;
t116 = t102 * t67 - t99 * t68;
t75 = (-pkin(9) - qJ(5)) * t95;
t77 = t97 * qJ(5) + t90;
t115 = t102 * t77 + t99 * t75;
t114 = t102 * t75 - t99 * t77;
t109 = -pkin(4) * t39 - qJ(5) * t40 - qJD(5) * t54;
t108 = t55 * t132 + t39 * t87 - t40 * t84 - t54 * t79;
t56 = -0.2e1 * t70 * t64;
t49 = -t70 * qJD(5) - t115 * qJD(6);
t48 = -t113 * qJD(5) - t114 * qJD(6);
t46 = -0.2e1 * t113 * t64 - 0.2e1 * t70 * t65;
t42 = -t117 * qJD(6) - t70 * t79;
t41 = -t116 * qJD(6) - t113 * t79;
t38 = t113 * t55;
t37 = t70 * t55;
t29 = t113 * t40 - t65 * t54;
t17 = t20 * t95;
t15 = t39 * t146 - t55 * t131 + (t55 * t134 + t149) * t102;
t14 = t113 * t39 - t55 * t65;
t9 = t14 * t70 - t38 * t64;
t6 = -pkin(9) * t149 + t8;
t5 = t40 * pkin(5) - pkin(9) * t148 + t7;
t4 = t113 * t14 - t70 * t15 + t64 * t37 - t38 * t65;
t2 = -t119 * qJD(6) + t102 * t5 - t99 * t6;
t1 = -t118 * qJD(6) - t102 * t6 - t99 * t5;
t10 = [0, 0, 0, 0, 0, t123, qJ(2) * t123, -0.2e1 * t71 * t66, 0.2e1 * t66 * t110 - 0.2e1 * t71 * t139, 0, 0, 0, t139 * t156, -t66 * t156, 0.2e1 * t55 * t39, -0.2e1 * t39 * t54 - 0.2e1 * t55 * t40, 0, 0, 0, 0.2e1 * t54 * t130 + 0.2e1 * t59 * t40, 0.2e1 * t55 * t130 + 0.2e1 * t59 * t39, 0.2e1 * t121 * t95 + 0.2e1 * t24 * t40 + 0.2e1 * t7 * t54, 0.2e1 * t121 * t97 - 0.2e1 * t25 * t40 - 0.2e1 * t8 * t54, -0.2e1 * t122 * t55 + 0.2e1 * (-t24 * t97 - t25 * t95) * t39, 0.2e1 * t35 * t20 + 0.2e1 * t24 * t7 + 0.2e1 * t25 * t8, 0.2e1 * t38 * t14, -0.2e1 * t14 * t37 - 0.2e1 * t38 * t15, 0.2e1 * t14 * t54 + 0.2e1 * t38 * t40, -0.2e1 * t15 * t54 - 0.2e1 * t37 * t40, 0.2e1 * t54 * t40, 0.2e1 * t118 * t40 + 0.2e1 * t12 * t37 + 0.2e1 * t30 * t15 + 0.2e1 * t2 * t54, 0.2e1 * t1 * t54 - 0.2e1 * t119 * t40 + 0.2e1 * t12 * t38 + 0.2e1 * t30 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t139, -t66, 0, 0, 0, 0, 0, t40, t39, t97 * t40, -t95 * t40, -t137 * t39, t122, 0, 0, 0, 0, 0, t29, -t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t66, -t139, 0, t107, t157, 0, 0, t39, -t40, 0, -t20, t19, t108 * t95 - t153, t108 * t97 + t17, t3, t120 * t79 + t35 * t132 + t20 * t87 + t3 * t84, t9, t4, t28, t29, 0, t116 * t40 + t37 * t132 + t74 * t15 + t42 * t54 + t144, -t117 * t40 + t38 * t132 + t74 * t14 + t41 * t54 + t143; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t132, -0.2e1 * t124, -0.2e1 * t125, 0.2e1 * t126, 0.2e1 * t159, 0.2e1 * t87 * t132 + 0.2e1 * t159 * t84, t56, t46, 0, 0, 0, 0.2e1 * t142, 0.2e1 * t141; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, -t40, 0, -t20, t19, t109 * t95 - t153, t109 * t97 + t17, t3, -t20 * pkin(4) + t3 * qJ(5) + t120 * qJD(5), t9, t4, t28, t29, 0, t114 * t40 + t85 * t15 + t49 * t54 + t144, -t115 * t40 + t85 * t14 + t48 * t54 + t143; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t132, -t124, -t125, t126, t158 + t159, -pkin(4) * t132 + qJ(5) * t159 + t158 * t84, t56, t46, 0, 0, 0, t142 + t150, t141 - t151; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t160, qJ(5) * t160, t56, t46, 0, 0, 0, 0.2e1 * t150, -0.2e1 * t151; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t149, t148, 0, t20, 0, 0, 0, 0, 0, t15, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t132, 0, 0, 0, 0, 0, t65, -t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, -t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, -t15, t40, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t65, t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t64, -t65, 0, t42, t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t64, -t65, 0, t49, t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t10;
