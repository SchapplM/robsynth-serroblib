% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% tau_reg [5x18]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRRP5_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP5_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP5_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP5_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP5_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP5_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:40:59
% EndTime: 2019-12-31 18:41:01
% DurationCPUTime: 0.78s
% Computational Cost: add. (1285->183), mult. (2168->208), div. (0->0), fcn. (1205->12), ass. (0->121)
t82 = sin(qJ(4));
t85 = cos(qJ(4));
t103 = t85 * pkin(4) + t82 * qJ(5);
t40 = -pkin(3) - t103;
t80 = sin(pkin(8));
t164 = pkin(1) * t80;
t81 = cos(pkin(8));
t67 = t81 * pkin(1) + pkin(2);
t52 = t67 * qJD(1);
t170 = qJD(3) * t52 + qJDD(1) * t164;
t121 = qJD(1) * t164;
t172 = -qJD(3) * t121 + t67 * qJDD(1);
t83 = sin(qJ(3));
t86 = cos(qJ(3));
t165 = -t170 * t86 - t172 * t83;
t75 = qJDD(1) + qJDD(3);
t173 = t75 * pkin(7) + qJD(2) * qJD(4) - t165;
t29 = t86 * t121 + t83 * t52;
t76 = qJD(1) + qJD(3);
t25 = t76 * pkin(7) + t29;
t146 = t82 * t25;
t17 = t85 * qJD(2) - t146;
t171 = qJD(5) - t17;
t124 = t82 * qJDD(2) + t173 * t85;
t125 = qJDD(4) * qJ(5);
t3 = t125 + (qJD(5) - t146) * qJD(4) + t124;
t132 = qJD(4) * t85;
t111 = -t85 * qJDD(2) + t25 * t132 + t173 * t82;
t131 = qJDD(4) * pkin(4);
t167 = qJDD(5) - t131;
t4 = t111 + t167;
t169 = t3 * t85 + t4 * t82;
t77 = qJ(1) + pkin(8);
t73 = qJ(3) + t77;
t65 = sin(t73);
t66 = cos(t73);
t168 = -t66 * pkin(7) - t40 * t65;
t141 = g(1) * t66 + g(2) * t65;
t140 = t86 * t164 + t83 * t67;
t14 = -qJD(4) * pkin(4) + t171;
t128 = qJD(4) * qJ(5);
t18 = t82 * qJD(2) + t85 * t25;
t15 = t18 + t128;
t113 = -t83 * t164 + t86 * t67;
t101 = t14 * t82 + t15 * t85;
t130 = t15 * qJD(4);
t166 = -t82 * t130 + t14 * t132 + t169;
t88 = qJD(4) ^ 2;
t163 = pkin(7) * t88;
t60 = g(1) * t65;
t162 = g(2) * t66;
t161 = t75 * pkin(3);
t160 = t76 * pkin(3);
t28 = -t83 * t121 + t86 * t52;
t156 = t28 * t76;
t155 = t29 * t76;
t30 = t113 * qJD(3);
t154 = t30 * t76;
t31 = t140 * qJD(3);
t153 = t31 * t76;
t34 = pkin(7) + t140;
t152 = t34 * t88;
t151 = t40 * t75;
t150 = t40 * t76;
t149 = t65 * t82;
t148 = t66 * t82;
t147 = t76 * t82;
t145 = t85 * t75;
t129 = t82 * qJD(5);
t133 = qJD(4) * t82;
t35 = pkin(4) * t133 - t85 * t128 - t129;
t143 = t35 - t29;
t142 = g(1) * t149 - g(2) * t148;
t78 = t82 ^ 2;
t79 = t85 ^ 2;
t139 = t78 - t79;
t138 = t78 + t79;
t136 = pkin(7) * qJDD(4);
t134 = qJD(4) * t76;
t127 = qJDD(4) * t34;
t74 = t76 ^ 2;
t123 = t82 * t74 * t85;
t47 = t85 * t60;
t122 = t28 * t133 + t85 * t155 + t47;
t110 = -t170 * t83 + t172 * t86;
t10 = -t110 - t161;
t118 = -t10 - t162;
t117 = t138 * t75;
t116 = t17 + t146;
t26 = -t113 + t40;
t115 = t26 * t76 - t30;
t24 = -t28 - t160;
t112 = t10 * t82 + t24 * t132 - t142;
t109 = t65 * pkin(7) - t40 * t66;
t106 = -t161 + t163;
t84 = sin(qJ(1));
t87 = cos(qJ(1));
t104 = g(1) * t84 - g(2) * t87;
t102 = pkin(4) * t82 - qJ(5) * t85;
t5 = (t102 * qJD(4) - t129) * t76 + t151 - t110;
t100 = -t151 - t5 - t163;
t98 = t110 + t60 - t162;
t33 = -pkin(3) - t113;
t96 = t33 * t75 + t152 + t153;
t95 = g(1) * t148 + g(2) * t149 - g(3) * t85 - t111;
t94 = -t141 + t166;
t16 = t31 + t35;
t93 = -t16 * t76 - t26 * t75 - t152 - t5;
t92 = -t127 + (t33 * t76 - t30) * qJD(4);
t91 = t18 * qJD(4) + t95;
t90 = t141 + t165;
t64 = t82 * t75;
t42 = qJDD(4) * t85 - t88 * t82;
t41 = qJDD(4) * t82 + t88 * t85;
t36 = t102 * t76;
t32 = 0.2e1 * t132 * t147 + t78 * t75;
t22 = -0.2e1 * t139 * t134 + 0.2e1 * t82 * t145;
t19 = t24 * t133;
t13 = -t28 + t150;
t11 = t13 * t133;
t1 = [qJDD(1), t104, g(1) * t87 + g(2) * t84, (t104 + (t80 ^ 2 + t81 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), t75, t113 * t75 - t153 + t98, -t140 * t75 - t154 + t90, t32, t22, t41, t42, 0, t19 + t47 + t92 * t82 + (t118 - t96) * t85, t82 * t96 + t85 * t92 + t112, t11 + t47 + (t115 * qJD(4) - t127) * t82 + (t93 - t162) * t85, t34 * t117 + t138 * t154 + t94, (t127 + (-t115 - t13) * qJD(4)) * t85 + t93 * t82 + t142, t5 * t26 + t13 * t16 - g(1) * (-pkin(2) * sin(t77) - t84 * pkin(1) - t168) - g(2) * (pkin(2) * cos(t77) + t87 * pkin(1) + t109) + t101 * t30 + t166 * t34; 0, 0, 0, qJDD(2) - g(3), 0, 0, 0, 0, 0, 0, 0, 0, t42, -t41, t42, 0, t41, t101 * qJD(4) + t3 * t82 - t4 * t85 - g(3); 0, 0, 0, 0, t75, t98 + t155, t90 + t156, t32, t22, t41, t42, 0, t19 + (-pkin(3) * t134 - t136) * t82 + (-t106 + t118) * t85 + t122, (-t136 + (t28 - t160) * qJD(4)) * t85 + (t106 - t155) * t82 + t112, t11 + (t40 * t134 - t136) * t82 + (-t35 * t76 + t100 - t162) * t85 + t122, pkin(7) * t117 - t138 * t156 + t94, (t136 + (-t13 - t28 - t150) * qJD(4)) * t85 + (-t143 * t76 + t100) * t82 + t142, t5 * t40 - g(2) * t109 - t101 * t28 + t143 * t13 + ((t14 * t85 - t15 * t82) * qJD(4) + t169) * pkin(7) + t168 * g(1); 0, 0, 0, 0, 0, 0, 0, -t123, t139 * t74, t64, t145, qJDD(4), -t24 * t147 + t91, g(3) * t82 + t116 * qJD(4) + (-t24 * t76 + t141) * t85 - t124, 0.2e1 * t131 - qJDD(5) + (-t13 * t82 + t36 * t85) * t76 + t91, -t102 * t75, 0.2e1 * t125 + (t36 * t76 - g(3)) * t82 + (t13 * t76 - t141) * t85 + (0.2e1 * qJD(5) - t116) * qJD(4) + t124, -t4 * pkin(4) - g(3) * t103 + t3 * qJ(5) + t141 * t102 - t13 * t36 - t14 * t18 + t171 * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(4) - t123, t64, -t78 * t74 - t88, t13 * t147 - t130 + t167 - t95;];
tau_reg = t1;
