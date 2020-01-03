% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RPRRR8
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RPRRR8_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR8_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR8_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR8_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR8_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR8_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:06:06
% EndTime: 2019-12-31 19:06:10
% DurationCPUTime: 1.29s
% Computational Cost: add. (5707->216), mult. (7559->258), div. (0->0), fcn. (4000->8), ass. (0->141)
t118 = sin(qJ(5));
t112 = qJDD(4) + qJDD(5);
t115 = -qJD(1) + qJD(3);
t121 = cos(qJ(5));
t122 = cos(qJ(4));
t119 = sin(qJ(4));
t154 = t115 * t119;
t82 = -t121 * t122 * t115 + t118 * t154;
t84 = (t122 * t118 + t119 * t121) * t115;
t62 = t84 * t82;
t174 = -t62 + t112;
t177 = t118 * t174;
t176 = t121 * t174;
t151 = qJD(4) * t115;
t148 = t119 * t151;
t113 = qJDD(1) - qJDD(3);
t152 = t122 * t113;
t140 = -t148 - t152;
t147 = t122 * t151;
t153 = t119 * t113;
t89 = t147 - t153;
t49 = -t82 * qJD(5) + t118 * t140 + t121 * t89;
t114 = qJD(4) + qJD(5);
t76 = t114 * t82;
t175 = -t76 + t49;
t111 = t115 ^ 2;
t120 = sin(qJ(3));
t123 = cos(qJ(3));
t168 = sin(qJ(1));
t169 = cos(qJ(1));
t143 = t168 * g(1) - t169 * g(2);
t141 = -qJDD(2) + t143;
t171 = qJD(1) ^ 2;
t132 = t171 * qJ(2) + t141;
t170 = pkin(1) + pkin(2);
t126 = -t170 * qJDD(1) - t132;
t116 = qJDD(1) * qJ(2);
t136 = t169 * g(1) + t168 * g(2);
t133 = (2 * qJD(2) * qJD(1)) - t136;
t131 = t116 + t133;
t79 = -t170 * t171 + t131;
t56 = t120 * t126 + t123 * t79;
t52 = -t111 * pkin(3) - t113 * pkin(7) + t56;
t46 = -t122 * g(3) + t119 * t52;
t173 = -t46 + (-t89 + t147) * pkin(8);
t80 = t82 ^ 2;
t81 = t84 ^ 2;
t110 = t114 ^ 2;
t172 = t122 ^ 2;
t167 = t119 * g(3);
t157 = t111 * t122;
t99 = t119 * t157;
t150 = qJDD(4) + t99;
t125 = t150 * pkin(4) + t173;
t149 = pkin(8) * t154;
t97 = qJD(4) * pkin(4) - t149;
t32 = t167 + (-t97 - t149) * qJD(4) + (-pkin(4) * t157 - t113 * pkin(8) + t52) * t122;
t16 = t118 * t32 - t121 * t125;
t162 = t121 * t32;
t17 = t118 * t125 + t162;
t6 = t118 * t17 - t121 * t16;
t166 = t119 * t6;
t104 = t172 * t111;
t145 = t120 * t79 - t123 * t126;
t51 = t113 * pkin(3) - t111 * pkin(7) + t145;
t34 = -t140 * pkin(4) - pkin(8) * t104 + t97 * t154 + t51;
t165 = t118 * t34;
t59 = t62 + t112;
t164 = t118 * t59;
t163 = t119 * t150;
t161 = t121 * t34;
t160 = t121 * t59;
t159 = t122 * (qJDD(4) - t99);
t158 = qJDD(1) * pkin(1);
t156 = t114 * t118;
t155 = t114 * t121;
t7 = t118 * t16 + t121 * t17;
t146 = t118 * t89 - t121 * t140;
t47 = t122 * t52 + t167;
t29 = t119 * t46 + t122 * t47;
t144 = -pkin(3) * t51 + pkin(7) * t29;
t92 = -t123 * t111 + t120 * t113;
t93 = t120 * t111 + t123 * t113;
t139 = 0.2e1 * t148 + t152;
t117 = t119 ^ 2;
t103 = t117 * t111;
t124 = qJD(4) ^ 2;
t70 = -t159 - t119 * (-t103 - t124);
t88 = 0.2e1 * t147 - t153;
t138 = pkin(3) * t88 - pkin(7) * t70 - t119 * t51;
t69 = t122 * (-t104 - t124) - t163;
t137 = -pkin(3) * t139 + pkin(7) * t69 - t122 * t51;
t135 = (-qJD(5) + t114) * t84 - t146;
t91 = (-t117 - t172) * t113;
t94 = t103 + t104;
t134 = pkin(3) * t94 + pkin(7) * t91 + t29;
t42 = t76 + t49;
t24 = t118 * t135 - t121 * t42;
t25 = t118 * t42 + t121 * t135;
t10 = -t119 * t24 + t122 * t25;
t50 = -t80 - t81;
t130 = pkin(3) * t50 - pkin(7) * t10 - t119 * (-pkin(8) * t24 - t6) - t122 * (-pkin(4) * t50 + pkin(8) * t25 + t7);
t57 = -t110 - t80;
t35 = t118 * t57 + t176;
t36 = t121 * t57 - t177;
t20 = -t119 * t35 + t122 * t36;
t37 = (qJD(5) + t114) * t84 + t146;
t129 = pkin(3) * t37 - pkin(7) * t20 - t122 * (-pkin(4) * t37 + pkin(8) * t36 - t161) - t119 * (-pkin(8) * t35 + t165);
t71 = -t81 - t110;
t43 = t121 * t71 - t164;
t44 = -t118 * t71 - t160;
t26 = -t119 * t43 + t122 * t44;
t128 = pkin(3) * t175 - pkin(7) * t26 - t119 * (-pkin(8) * t43 + t161) - t122 * (-pkin(4) * t175 + pkin(8) * t44 + t165);
t2 = t122 * t7 - t166;
t127 = pkin(3) * t34 - pkin(7) * t2 + pkin(8) * t166 - t122 * (-pkin(4) * t34 + pkin(8) * t7);
t87 = t132 + t158;
t74 = -t81 + t110;
t73 = t80 - t110;
t68 = t163 + t122 * (-t103 + t124);
t67 = t119 * (t104 - t124) + t159;
t66 = (t89 + t147) * t119;
t65 = t139 * t122;
t64 = t120 * t91 + t123 * t94;
t63 = -t119 * t139 + t122 * t88;
t61 = t81 - t80;
t54 = t120 * t70 - t123 * t88;
t53 = t120 * t69 - t123 * t139;
t48 = -t84 * qJD(5) - t146;
t33 = t120 * t56 - t123 * t145;
t30 = (t119 * (t118 * t84 - t121 * t82) + t122 * (-t118 * t82 - t121 * t84)) * t114;
t28 = t119 * (t121 * t73 - t164) + t122 * (t118 * t73 + t160);
t27 = t119 * (-t118 * t74 + t176) + t122 * (t121 * t74 + t177);
t22 = t119 * (t121 * t49 - t84 * t156) + t122 * (t118 * t49 + t84 * t155);
t21 = t119 * (-t118 * t48 + t82 * t155) + t122 * (t121 * t48 + t82 * t156);
t18 = t120 * t29 - t123 * t51;
t14 = t120 * t26 - t123 * t175;
t12 = t120 * t20 - t123 * t37;
t9 = t119 * (-t118 * t175 - t121 * t37) + t122 * (-t118 * t37 + t121 * t175);
t8 = t120 * t10 - t123 * t50;
t1 = t120 * t2 - t123 * t34;
t3 = [0, 0, 0, 0, 0, qJDD(1), t143, t136, 0, 0, 0, 0, 0, qJDD(1), 0, 0, t141 + 0.2e1 * t158, 0, 0.2e1 * t116 + t133, pkin(1) * t87 + qJ(2) * (-t171 * pkin(1) + t131), 0, 0, 0, 0, 0, t113, qJ(2) * t92 + t170 * t93 + t145, qJ(2) * t93 - t170 * t92 + t56, 0, qJ(2) * (t120 * t145 + t123 * t56) - t170 * t33, -t66, -t63, -t68, t65, -t67, 0, qJ(2) * (t120 * t139 + t123 * t69) - t170 * t53 - t137, qJ(2) * (t120 * t88 + t123 * t70) - t170 * t54 + t138, qJ(2) * (-t120 * t94 + t123 * t91) - t170 * t64 - t134, qJ(2) * (t120 * t51 + t123 * t29) - t170 * t18 - t144, -t22, -t9, -t27, -t21, -t28, -t30, qJ(2) * (t120 * t37 + t123 * t20) - t170 * t12 + t129, qJ(2) * (t120 * t175 + t123 * t26) - t170 * t14 + t128, qJ(2) * (t123 * t10 + t120 * t50) - t170 * t8 + t130, qJ(2) * (t120 * t34 + t123 * t2) - t170 * t1 + t127; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(1), 0, -t171, -t87, 0, 0, 0, 0, 0, 0, -t93, t92, 0, t33, 0, 0, 0, 0, 0, 0, t53, t54, t64, t18, 0, 0, 0, 0, 0, 0, t12, t14, t8, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t113, -t145, -t56, 0, 0, t66, t63, t68, -t65, t67, 0, t137, -t138, t134, t144, t22, t9, t27, t21, t28, t30, -t129, -t128, -t130, -t127; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t99, t103 - t104, -t153, t99, -t152, qJDD(4), -t46, -t47, 0, 0, t62, t61, t42, -t62, t135, t112, pkin(4) * t35 - t16, -t162 - t118 * t173 + (-t118 * t150 + t43) * pkin(4), pkin(4) * t24, pkin(4) * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, t61, t42, -t62, t135, t112, -t16, -t17, 0, 0;];
tauJ_reg = t3;
