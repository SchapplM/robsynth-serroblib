% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5PRRRP6
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:53
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRRRP6_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP6_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP6_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP6_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP6_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP6_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:52:24
% EndTime: 2019-12-05 16:52:29
% DurationCPUTime: 2.49s
% Computational Cost: add. (2464->323), mult. (5433->402), div. (0->0), fcn. (3806->10), ass. (0->175)
t112 = qJD(3) + qJD(4);
t118 = sin(qJ(4));
t121 = cos(qJ(3));
t231 = cos(qJ(4));
t172 = qJDD(2) * t231;
t119 = sin(qJ(3));
t190 = t119 * qJDD(2);
t204 = t118 * t121;
t73 = t231 * t119 + t204;
t241 = t112 * t73;
t25 = qJD(2) * t241 + t118 * t190 - t121 * t172;
t66 = t73 * qJD(2);
t12 = -t66 * t112 + t25;
t115 = qJ(3) + qJ(4);
t109 = sin(t115);
t111 = qJDD(3) + qJDD(4);
t120 = sin(qJ(2));
t116 = sin(pkin(8));
t117 = cos(pkin(8));
t160 = g(1) * t117 + g(2) * t116;
t150 = t160 * t120;
t122 = cos(qJ(2));
t226 = g(3) * t122;
t138 = t150 - t226;
t205 = t118 * t119;
t146 = t231 * t121 - t205;
t144 = t146 * t122;
t232 = pkin(7) + pkin(6);
t80 = t232 * t119;
t81 = t232 * t121;
t151 = -t118 * t81 - t231 * t80;
t181 = qJD(3) * t232;
t75 = t119 * t181;
t221 = -qJD(1) * t144 + t151 * qJD(4) - t181 * t204 - t231 * t75;
t46 = -t118 * t80 + t231 * t81;
t242 = -t138 * t109 - t46 * t111 - t221 * t112;
t149 = t160 * t122;
t227 = g(3) * t120;
t139 = t149 + t227;
t113 = t119 ^ 2;
t114 = t121 ^ 2;
t200 = t113 + t114;
t167 = t122 * t200;
t194 = t120 * qJD(1);
t88 = qJD(2) * pkin(6) + t194;
t193 = t122 * qJD(1);
t218 = qJD(2) * pkin(2);
t89 = -t193 - t218;
t239 = t89 * t120 + t88 * t167;
t156 = t112 * t205;
t106 = t111 * qJ(5);
t107 = t112 * qJD(5);
t237 = t106 + t107;
t108 = t111 * pkin(4);
t236 = qJDD(5) - t108;
t124 = qJD(3) ^ 2;
t192 = qJD(1) * qJD(2);
t212 = qJDD(2) * pkin(2);
t165 = -t122 * qJDD(1) + t120 * t192;
t69 = t165 - t212;
t234 = -pkin(6) * t124 + t120 * (t160 + t192) + t212 - t226 - t69;
t233 = t66 ^ 2;
t228 = g(2) * t117;
t225 = t121 * pkin(3);
t179 = qJD(2) * t231;
t163 = t121 * t179;
t198 = qJD(2) * t119;
t180 = t118 * t198;
t64 = -t163 + t180;
t105 = pkin(2) + t225;
t71 = -t105 * qJD(2) - t193;
t26 = t64 * pkin(4) - t66 * qJ(5) + t71;
t224 = t26 * t64;
t223 = t66 * t64;
t222 = t71 * t64;
t178 = t231 * qJD(3);
t164 = t121 * t178;
t220 = t46 * qJD(4) - t118 * t75 + t232 * t164 - t73 * t193;
t177 = t231 * qJD(4);
t174 = pkin(7) * qJD(2) + t88;
t61 = t174 * t121;
t217 = t118 * t61;
t60 = t174 * t119;
t36 = -t231 * t60 - t217;
t219 = pkin(3) * t177 + qJD(5) - t36;
t184 = t231 * t61;
t52 = qJD(3) * pkin(3) - t60;
t33 = t118 * t52 + t184;
t216 = t33 * t112;
t211 = t109 * t120;
t110 = cos(t115);
t210 = t110 * t120;
t209 = t116 * t110;
t208 = t116 * t122;
t207 = t117 * t110;
t206 = t117 * t122;
t32 = t231 * t52 - t217;
t203 = qJD(5) - t32;
t202 = qJDD(1) - g(3);
t201 = t113 - t114;
t125 = qJD(2) ^ 2;
t199 = t124 + t125;
t197 = qJD(2) * t120;
t196 = qJD(3) * t119;
t195 = qJD(4) * t118;
t191 = qJD(2) * qJD(3);
t189 = t119 * qJDD(3);
t188 = t121 * qJDD(2);
t187 = t122 * qJDD(2);
t186 = t121 * t228;
t185 = pkin(3) * t196;
t27 = t64 ^ 2 - t233;
t182 = t119 * t125 * t121;
t176 = t119 * t191;
t175 = t121 * t191;
t58 = t109 * t206 - t209;
t59 = t116 * t109 + t110 * t206;
t173 = -t58 * pkin(4) + t59 * qJ(5);
t70 = qJDD(2) * pkin(6) + t120 * qJDD(1) + t122 * t192;
t171 = t200 * t70;
t31 = -t121 * t88 * qJD(3) + qJDD(3) * pkin(3) - t119 * t70 + (-t175 - t190) * pkin(7);
t34 = -t88 * t196 + t121 * t70 + (-t176 + t188) * pkin(7);
t4 = t118 * t31 + t52 * t177 - t61 * t195 + t231 * t34;
t170 = t118 * t34 + t61 * t177 + t52 * t195 - t231 * t31;
t169 = -t112 * t163 - t118 * t188 - t119 * t172;
t40 = -t121 * t177 + t156 - t164;
t7 = pkin(4) * t241 + t40 * qJ(5) - t73 * qJD(5) + t185;
t168 = -t7 + t194;
t166 = t200 * qJDD(2);
t162 = t119 * t175;
t35 = -t118 * t60 + t184;
t161 = -pkin(3) * t195 + t35;
t38 = t66 * pkin(4) + t64 * qJ(5);
t159 = g(1) * t116 - t228;
t158 = -t146 * t25 + t241 * t64;
t155 = pkin(4) * t110 + qJ(5) * t109;
t154 = -t111 * t146 + t112 * t241;
t148 = t185 - t194;
t15 = qJD(2) * t144 - t120 * t241;
t16 = t119 * t122 * t179 + (qJD(2) * t118 * t122 + (t178 + t177) * t120) * t121 - t120 * t156;
t24 = qJD(2) * t156 + t169;
t62 = t73 * t120;
t63 = t146 * t120;
t147 = -t15 * t64 + t16 * t66 - t62 * t24 - t63 * t25;
t142 = -t62 * t111 - t16 * t112 - t122 * t25 + t64 * t197;
t57 = -t117 * t109 + t110 * t208;
t141 = g(1) * t59 + g(2) * t57 + g(3) * t210 - t4;
t56 = t109 * t208 + t207;
t140 = g(1) * t58 + g(2) * t56 + g(3) * t211 - t170;
t136 = -t146 * t24 - t241 * t66 - t73 * t25 + t40 * t64;
t42 = pkin(3) * t176 - t105 * qJDD(2) + t165;
t135 = -pkin(6) * qJDD(3) + (t193 + t89 - t218) * qJD(3);
t134 = t63 * t111 + t15 * t112 - t122 * t24 - t66 * t197;
t133 = -t71 * t66 + t140;
t132 = t32 * t112 + t141;
t131 = -g(2) * (-t56 * pkin(4) + t57 * qJ(5)) - g(3) * (-pkin(4) * t211 + qJ(5) * t210);
t130 = -t110 * t226 + t151 * t111 - t220 * t112 + (g(1) * t207 + g(2) * t209) * t120;
t129 = -t89 * qJD(2) + t139 - t70;
t128 = -t26 * t66 + t140 - t236;
t126 = t151 * t24 + t220 * t66 - t221 * t64 - t46 * t25 - t139;
t104 = -t231 * pkin(3) - pkin(4);
t101 = t118 * pkin(3) + qJ(5);
t99 = t116 * t225;
t83 = t122 * t105;
t39 = -pkin(4) * t146 - t73 * qJ(5) - t105;
t37 = pkin(3) * t198 + t38;
t21 = t112 * qJ(5) + t33;
t20 = t73 * t111 - t40 * t112;
t19 = -t112 * pkin(4) + t203;
t11 = -t169 + (-t180 + t64) * t112;
t6 = -t24 * t73 - t66 * t40;
t3 = t25 * pkin(4) + t24 * qJ(5) - t66 * qJD(5) + t42;
t2 = t170 + t236;
t1 = t4 + t237;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t202, 0, 0, 0, 0, 0, 0, -t125 * t120 + t187, -qJDD(2) * t120 - t125 * t122, 0, -g(3) + (t120 ^ 2 + t122 ^ 2) * qJDD(1), 0, 0, 0, 0, 0, 0, (-0.2e1 * t176 + t188) * t122 + (-t199 * t121 - t189) * t120, (-qJDD(3) * t120 - 0.2e1 * t122 * t191) * t121 + (t199 * t120 - t187) * t119, t120 * t166 + t125 * t167, t239 * qJD(2) + t120 * t171 - t69 * t122 - g(3), 0, 0, 0, 0, 0, 0, t142, -t134, t147, -t42 * t122 + t33 * t15 - t32 * t16 + t170 * t62 + t71 * t197 + t4 * t63 - g(3), 0, 0, 0, 0, 0, 0, t142, t147, t134, t1 * t63 - t3 * t122 + t21 * t15 + t19 * t16 + t26 * t197 + t2 * t62 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t202 * t122 + t150, -t202 * t120 + t149, 0, 0, t113 * qJDD(2) + 0.2e1 * t162, 0.2e1 * t119 * t188 - 0.2e1 * t201 * t191, t124 * t121 + t189, t114 * qJDD(2) - 0.2e1 * t162, qJDD(3) * t121 - t124 * t119, 0, t135 * t119 + t234 * t121, -t234 * t119 + t135 * t121, -t227 + t171 + pkin(6) * t166 + (-t200 * t192 - t160) * t122, (-t69 + t138) * pkin(2) + (t171 - t139) * pkin(6) - t239 * qJD(1), t6, t136, t20, t158, -t154, 0, -t105 * t25 - t146 * t42 + t148 * t64 + t241 * t71 + t130, t105 * t24 + t148 * t66 - t71 * t40 + t42 * t73 + t242, t146 * t4 + t170 * t73 - t241 * t33 + t32 * t40 + t126, t4 * t46 - t170 * t151 - t42 * t105 - g(3) * (t120 * t232 + t83) + t148 * t71 + t221 * t33 - t220 * t32 + t160 * (t105 * t120 - t122 * t232), t6, t20, -t136, 0, t154, t158, -t146 * t3 - t168 * t64 + t241 * t26 + t39 * t25 + t130, t1 * t146 - t19 * t40 + t2 * t73 - t21 * t241 + t126, t168 * t66 + t39 * t24 + t26 * t40 - t3 * t73 - t242, -g(3) * t83 + t1 * t46 - t2 * t151 + t26 * t7 + t3 * t39 + t221 * t21 + t220 * t19 + (-g(3) * t155 - t160 * t232) * t122 + (-g(3) * t232 - t26 * qJD(1) + t160 * (t105 + t155)) * t120; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t182, t201 * t125, t190, t182, t188, qJDD(3), t129 * t119 - t159 * t121, t159 * t119 + t129 * t121, 0, 0, t223, -t27, t11, -t223, -t12, t111, t35 * t112 + (t231 * t111 - t112 * t195 - t64 * t198) * pkin(3) + t133, t36 * t112 + t222 + (-t111 * t118 - t112 * t177 - t66 * t198) * pkin(3) + t141, (t33 - t35) * t66 + (-t32 + t36) * t64 + (t231 * t24 - t118 * t25 + (t118 * t66 - t231 * t64) * qJD(4)) * pkin(3), -g(1) * t99 + t32 * t35 - t33 * t36 + (t186 - t170 * t231 + t4 * t118 + (-t32 * t118 + t33 * t231) * qJD(4) + (-t71 * qJD(2) + t139) * t119) * pkin(3), t223, t11, t27, t111, t12, -t223, -t104 * t111 + t112 * t161 - t37 * t64 + t128, -t101 * t25 - t104 * t24 + (-t161 + t21) * t66 + (t19 - t219) * t64, t101 * t111 + t219 * t112 + t37 * t66 - t141 - t224 + t237, t1 * t101 + t2 * t104 - t26 * t37 - t19 * t35 - g(1) * (t173 + t99) + t219 * t21 + (t119 * t139 + t19 * t195 + t186) * pkin(3) + t131; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t223, -t27, t11, -t223, -t12, t111, t133 + t216, t132 + t222, 0, 0, t223, t11, t27, t111, t12, -t223, -t38 * t64 + t108 + t128 + t216, pkin(4) * t24 - t25 * qJ(5) + (t21 - t33) * t66 + (t19 - t203) * t64, t38 * t66 + 0.2e1 * t106 + 0.2e1 * t107 - t132 - t224, -t2 * pkin(4) - g(1) * t173 + t1 * qJ(5) - t19 * t33 + t203 * t21 - t26 * t38 + t131; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t111 + t223, t11, -t112 ^ 2 - t233, -t21 * t112 - t128;];
tau_reg = t5;
