% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6PRPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta4]';
% 
% Output:
% tauc_reg [6x24]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRPPRR3_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR3_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR3_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR3_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:23:02
% EndTime: 2019-03-08 19:23:06
% DurationCPUTime: 1.21s
% Computational Cost: add. (989->213), mult. (2268->332), div. (0->0), fcn. (1631->10), ass. (0->126)
t62 = sin(pkin(6));
t123 = qJD(1) * t62;
t71 = cos(qJ(2));
t100 = t71 * t123;
t37 = (qJD(3) + t100) * qJD(2);
t61 = sin(pkin(11));
t63 = cos(pkin(11));
t68 = sin(qJ(2));
t101 = t68 * t123;
t88 = qJD(2) * t101;
t18 = t61 * t37 - t63 * t88;
t24 = t61 * t100 - t63 * t101;
t92 = t61 * qJD(3) - t24;
t154 = t92 * qJD(2) + t18;
t31 = (-t61 * t71 + t63 * t68) * t62;
t70 = cos(qJ(5));
t111 = t70 * qJD(2);
t52 = qJD(6) + t111;
t67 = sin(qJ(5));
t59 = t67 ^ 2;
t66 = sin(qJ(6));
t153 = (qJD(2) * t59 - t52 * t70) * t66;
t72 = -pkin(2) - pkin(3);
t127 = t63 * qJ(3) + t61 * t72;
t42 = -pkin(8) + t127;
t73 = qJD(5) ^ 2;
t152 = -t42 * t73 + t154;
t74 = qJD(2) ^ 2;
t151 = (t73 + t74) * t61;
t150 = -qJD(6) + t52;
t85 = qJD(3) - t100;
t34 = t72 * qJD(2) + t85;
t46 = qJD(2) * qJ(3) + t101;
t17 = t61 * t34 + t63 * t46;
t15 = -qJD(2) * pkin(8) + t17;
t64 = cos(pkin(6));
t51 = -t64 * qJD(1) + qJD(4);
t10 = t70 * t15 + t67 * t51;
t19 = t63 * t37 + t61 * t88;
t4 = t10 * qJD(5) + t67 * t19;
t148 = t4 * t66;
t69 = cos(qJ(6));
t147 = t4 * t69;
t109 = qJD(5) * qJD(6);
t110 = qJD(2) * qJD(5);
t97 = t70 * t110;
t121 = qJD(2) * t67;
t99 = t66 * t121;
t22 = qJD(6) * t99 + (t109 - t97) * t69;
t146 = t22 * t66;
t145 = t22 * t70;
t116 = qJD(6) * t69;
t102 = t67 * t116;
t113 = t66 * qJD(5);
t23 = -t66 * t109 + (t70 * t113 + t102) * qJD(2);
t144 = t23 * t70;
t112 = t69 * qJD(5);
t39 = t99 + t112;
t143 = t39 * t52;
t40 = t69 * t121 - t113;
t142 = t40 * t52;
t141 = t46 * t71;
t139 = t62 * t74;
t137 = t66 * t52;
t136 = t66 * t70;
t135 = t67 * t22;
t134 = t67 * t39;
t133 = t67 * t40;
t132 = t69 * t52;
t131 = t69 * t70;
t130 = t73 * t67;
t129 = t73 * t70;
t86 = -pkin(5) * t67 + pkin(9) * t70;
t77 = t86 * qJD(5);
t128 = -t77 - t92;
t126 = -t70 ^ 2 + t59;
t124 = qJD(2) * pkin(2);
t120 = qJD(5) * t42;
t119 = qJD(5) * t67;
t118 = qJD(5) * t70;
t117 = qJD(6) * t66;
t115 = qJD(6) * t70;
t108 = t68 * t139;
t107 = t71 * t139;
t106 = t52 * t131;
t105 = t52 * t117;
t104 = t67 * t117;
t103 = t52 * t116;
t8 = qJD(5) * pkin(9) + t10;
t98 = t42 * t52 + t8;
t96 = t67 * t110;
t16 = t63 * t34 - t61 * t46;
t94 = -t61 * qJ(3) + t63 * t72;
t14 = qJD(2) * pkin(4) - t16;
t93 = qJD(2) * t14 - t19;
t30 = (t61 * t68 + t63 * t71) * t62;
t27 = qJD(1) * t30;
t91 = t63 * qJD(3) - t27;
t90 = t52 * t102;
t89 = 0.2e1 * t96;
t41 = pkin(4) - t94;
t87 = t70 * pkin(5) + t67 * pkin(9);
t11 = t87 * qJD(2) + t14;
t1 = t69 * t11 - t66 * t8;
t2 = t66 * t11 + t69 * t8;
t84 = t67 * t15 - t70 * t51;
t83 = t16 * t61 - t17 * t63;
t21 = t31 * t70 - t64 * t67;
t82 = t21 * t69 + t30 * t66;
t81 = -t21 * t66 + t30 * t69;
t20 = t31 * t67 + t64 * t70;
t78 = t69 * t59 * t110 + t52 * t104;
t76 = qJD(5) * (-qJD(2) * t41 - t14 - t91);
t3 = -t84 * qJD(5) + t70 * t19;
t7 = -qJD(5) * pkin(5) + t84;
t75 = -t7 * qJD(5) - qJD(6) * t11 - t91 * t52 - t3;
t45 = t86 * qJD(2);
t38 = t85 - t124;
t28 = t41 + t87;
t26 = qJD(2) * t30;
t25 = qJD(2) * t31;
t13 = qJD(2) * t77 + t18;
t12 = t69 * t13;
t6 = -t20 * qJD(5) + t26 * t70;
t5 = t21 * qJD(5) + t26 * t67;
t9 = [0, 0, -t108, -t107, -t108, t107 (t37 * t68 + (t141 + (t38 - t100) * t68) * qJD(2)) * t62, -t25 * qJD(2), t26 * qJD(2), t16 * t25 + t17 * t26 + t18 * t30 + t19 * t31, 0, 0, 0, 0, 0, -t5 * qJD(5) + (-t30 * t119 - t25 * t70) * qJD(2), -t6 * qJD(5) + (-t118 * t30 + t25 * t67) * qJD(2), 0, 0, 0, 0, 0 (-qJD(6) * t82 - t25 * t69 - t6 * t66) * t52 - t81 * t96 - t5 * t39 - t20 * t23 -(qJD(6) * t81 - t25 * t66 + t6 * t69) * t52 + t82 * t96 - t5 * t40 + t20 * t22; 0, 0, 0, 0, 0, 0.2e1 * qJD(2) * qJD(3), t37 * qJ(3) + t46 * qJD(3) + (-t141 + (-t38 - t124) * t68) * t123, t154, t91 * qJD(2) + t19, -t83 * qJD(3) + t19 * t127 + t16 * t24 - t17 * t27 - t18 * t94, t70 * t89, -0.2e1 * t126 * t110, -t129, t130, 0, t152 * t70 + t67 * t76, -t152 * t67 + t70 * t76, -t69 * t135 + (t112 * t70 - t104) * t40 (-t39 * t69 - t40 * t66) * t118 + (t146 - t23 * t69 + (t39 * t66 - t40 * t69) * qJD(6)) * t67, t145 + (-t106 + t133) * qJD(5) + t78, t90 + t144 + (-t134 - t153) * qJD(5) (-t52 - t111) * t119 (-t28 * t117 - t128 * t69) * t52 + (-t116 * t98 - t39 * t120 + t66 * t75 + t12) * t70 + (-t7 * t116 - t42 * t23 - t148 - t91 * t39 + (t42 * t137 - (-t42 * t136 + t69 * t28) * qJD(2) - t1) * qJD(5)) * t67 (-t28 * t116 + t128 * t66) * t52 + (-t40 * t120 + (qJD(6) * t98 - t13) * t66 + t75 * t69) * t70 + (t7 * t117 + t42 * t22 - t147 - t91 * t40 + (t42 * t132 + (t42 * t131 + t66 * t28) * qJD(2) + t2) * qJD(5)) * t67; 0, 0, 0, 0, 0, -t74 (-t46 + t101) * qJD(2), -t61 * t74, -t63 * t74, t83 * qJD(2) - t18 * t63 + t19 * t61, 0, 0, 0, 0, 0, -t151 * t70 + t63 * t89, t151 * t67 + 0.2e1 * t63 * t97, 0, 0, 0, 0, 0, t63 * t105 + ((t113 * t67 - t115 * t69) * t52 - t39 * t118 - t67 * t23) * t61 + (-(-t63 * t136 + t61 * t69) * t52 + (-(-t61 * t136 - t69 * t63) * qJD(5) + t63 * t39) * t67) * qJD(2), t63 * t103 + (-(-t112 * t67 - t115 * t66) * t52 - t40 * t118 + t135) * t61 + ((t63 * t131 + t61 * t66) * t52 + ((t61 * t131 - t66 * t63) * qJD(5) + t63 * t40) * t67) * qJD(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t130, -t129, 0, 0, 0, 0, 0, -t90 + t144 + (-t134 + t153) * qJD(5), -t145 + (-t106 - t133) * qJD(5) + t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t67 * t74 * t70, t126 * t74, 0, 0, 0, t93 * t67, t93 * t70, -t40 * t132 + t146 (t22 + t143) * t69 + (t23 + t142) * t66, t103 + (t106 + (-t40 - t113) * t67) * qJD(2), -t105 + (-t52 * t136 + (t39 - t112) * t67) * qJD(2), t52 * t121, pkin(5) * t23 - t147 - (t69 * t45 + t66 * t84) * t52 + t10 * t39 + (-pkin(9) * t132 + t7 * t66) * qJD(6) + (t1 * t67 + (pkin(9) * t119 + t7 * t70) * t66) * qJD(2), -pkin(5) * t22 + t148 + (t66 * t45 - t69 * t84) * t52 + t10 * t40 + (pkin(9) * t137 + t7 * t69) * qJD(6) + (t7 * t131 + (pkin(9) * t112 - t2) * t67) * qJD(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40 * t39, -t39 ^ 2 + t40 ^ 2, t22 - t143, t23 - t142, -t96, t150 * t2 - t66 * t3 + t7 * t40 + t12, t1 * t150 - t66 * t13 - t69 * t3 - t7 * t39;];
tauc_reg  = t9;
