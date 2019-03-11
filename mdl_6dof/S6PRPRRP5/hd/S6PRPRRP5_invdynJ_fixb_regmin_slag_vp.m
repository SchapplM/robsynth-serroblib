% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6PRPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
% 
% Output:
% tau_reg [6x23]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PRPRRP5_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP5_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP5_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRP5_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP5_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP5_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:16:46
% EndTime: 2019-03-08 20:16:53
% DurationCPUTime: 2.50s
% Computational Cost: add. (2275->358), mult. (4745->486), div. (0->0), fcn. (3472->10), ass. (0->186)
t114 = sin(qJ(4));
t117 = cos(qJ(4));
t176 = t117 * qJDD(2);
t183 = qJD(2) * qJD(4);
t246 = -t114 * t183 + t176;
t113 = sin(qJ(5));
t232 = pkin(5) * t113;
t245 = pkin(8) + t232;
t116 = cos(qJ(5));
t118 = cos(qJ(2));
t109 = sin(pkin(6));
t196 = qJD(1) * t109;
t115 = sin(qJ(2));
t205 = t113 * t115;
t147 = pkin(4) * t117 + pkin(9) * t114;
t73 = t147 * qJD(4) + qJD(3);
t244 = -(-t114 * t205 + t116 * t118) * t196 + t116 * t73;
t209 = t109 * t118;
t108 = sin(pkin(10));
t110 = cos(pkin(10));
t111 = cos(pkin(6));
t206 = t111 * t118;
t59 = t108 * t115 - t110 * t206;
t61 = t108 * t206 + t110 * t115;
t129 = g(1) * t61 + g(2) * t59 - g(3) * t209;
t119 = -pkin(2) - pkin(8);
t185 = t116 * qJD(4);
t188 = qJD(5) * t116;
t202 = t115 * t116;
t81 = t114 * pkin(4) - t117 * pkin(9) + qJ(3);
t243 = (t113 * t118 + t114 * t202) * t196 - t117 * t119 * t185 - t113 * t73 - t81 * t188;
t208 = t111 * t114;
t163 = t118 * t196;
t144 = qJD(3) - t163;
t65 = t119 * qJD(2) + t144;
t242 = -qJD(1) * t208 + t117 * t65;
t210 = t109 * t117;
t30 = t108 * t210 + t61 * t114;
t32 = t110 * t210 - t59 * t114;
t64 = t111 * t117 - t114 * t209;
t35 = t109 * t202 - t64 * t113;
t207 = t111 * t115;
t60 = t108 * t118 + t110 * t207;
t62 = -t108 * t207 + t110 * t118;
t241 = -g(1) * (-t30 * t113 + t62 * t116) - g(2) * (t32 * t113 + t60 * t116) - g(3) * t35;
t187 = t113 * qJD(4);
t192 = qJD(2) * t117;
t77 = t116 * t192 + t187;
t216 = qJD(5) * t77;
t23 = -t116 * qJDD(4) + t246 * t113 + t216;
t148 = g(1) * t62 + g(2) * t60;
t199 = qJDD(1) - g(3);
t211 = t109 * t115;
t240 = -t199 * t211 + t148;
t162 = t114 * t185;
t189 = qJD(5) * t113;
t165 = t117 * t189;
t22 = (t162 + t165) * qJD(2) - qJD(5) * t185 - t113 * qJDD(4) - t116 * t176;
t164 = t115 * t196;
t194 = qJD(2) * qJ(3);
t80 = t164 + t194;
t239 = qJD(4) * (t164 - t80 - t194) - qJDD(4) * t119;
t130 = g(3) * t211 + t148;
t195 = qJD(1) * t117;
t94 = t111 * t195;
t39 = t114 * t65 + t94;
t27 = qJD(4) * pkin(9) + t39;
t186 = t114 * qJD(2);
t99 = qJD(5) + t186;
t238 = (t119 * t99 + t27) * qJD(5) + t130;
t237 = t77 ^ 2;
t45 = t81 * qJD(2) + t164;
t15 = -t113 * t27 + t116 * t45;
t7 = -t77 * qJ(6) + t15;
t4 = t99 * pkin(5) + t7;
t236 = -t7 + t4;
t204 = t113 * t119;
t156 = pkin(5) - t204;
t184 = t116 * qJD(6);
t203 = t114 * t116;
t225 = t113 * t81 + t119 * t203;
t234 = qJ(6) * t162 - t225 * qJD(5) + (qJ(6) * t189 + t156 * qJD(4) - t184) * t117 + t244;
t220 = qJ(6) * t117;
t168 = t116 * t220;
t233 = -qJD(5) * t168 + (-qJD(6) * t117 + (qJ(6) * qJD(4) - qJD(5) * t119) * t114) * t113 - t243;
t75 = t113 * t192 - t185;
t231 = t75 * t99;
t230 = t77 * t99;
t229 = qJ(6) + pkin(9);
t79 = t147 * qJD(2);
t228 = t113 * t79 + t116 * t242;
t151 = qJD(5) * t229;
t227 = t184 - t228 + (-qJ(6) * t186 - t151) * t113;
t67 = t116 * t79;
t226 = -t116 * t151 - t67 - (pkin(5) * t117 + qJ(6) * t203) * qJD(2) + (-qJD(6) + t242) * t113;
t224 = t113 * t99;
t223 = t116 * t77;
t221 = t22 * t113;
t219 = qJD(2) * t99;
t218 = qJD(4) * t75;
t217 = qJD(4) * t77;
t215 = qJD(5) * t99;
t214 = qJDD(2) * pkin(2);
t101 = t116 * pkin(5) + pkin(4);
t213 = t101 * t114;
t212 = t109 * t114;
t121 = qJD(2) ^ 2;
t201 = t118 * t121;
t200 = t80 * qJD(2);
t107 = t117 ^ 2;
t198 = t114 ^ 2 - t107;
t120 = qJD(4) ^ 2;
t197 = -t120 - t121;
t193 = qJD(2) * t109;
t191 = qJD(4) * t114;
t190 = qJD(4) * t117;
t182 = qJDD(1) * t109;
t181 = qJDD(1) * t111;
t180 = qJDD(2) * qJ(3);
t179 = qJDD(4) * t114;
t177 = t114 * qJDD(2);
t157 = t117 * t181;
t159 = t118 * t182;
t87 = qJD(2) * t164;
t137 = qJDD(3) + t87 - t159;
t41 = t119 * qJDD(2) + t137;
t12 = qJDD(4) * pkin(9) + qJD(4) * t242 + t114 * t41 + t157;
t160 = t115 * t182;
t20 = t160 + t81 * qJDD(2) + (t73 + t163) * qJD(2);
t175 = -t113 * t20 - t116 * t12 - t45 * t188;
t174 = -qJD(4) * t94 - t114 * t181 - t65 * t191;
t171 = g(3) * (pkin(2) * t209 + qJ(3) * t211);
t170 = t99 * t187;
t169 = t99 * t185;
t167 = t115 * t193;
t166 = t118 * t193;
t158 = t117 * t183;
t155 = -t119 + t232;
t154 = -t59 * pkin(2) + t60 * qJ(3);
t153 = -t61 * pkin(2) + t62 * qJ(3);
t150 = -t41 + t200;
t16 = t113 * t45 + t116 * t27;
t8 = -t75 * qJ(6) + t16;
t146 = t113 * t8 + t116 * t4;
t145 = t113 * t4 - t116 * t8;
t143 = t115 * (-qJD(2) * pkin(2) + t144) + t118 * t80;
t141 = qJDD(2) * t115 + t201;
t26 = -qJD(4) * pkin(4) - t242;
t72 = qJDD(5) + t158 + t177;
t139 = t113 * t72 + t99 * t188;
t138 = t116 * t72 - t99 * t189;
t36 = t109 * t205 + t64 * t116;
t63 = t117 * t209 + t208;
t135 = t27 * t189 + t175;
t29 = t108 * t212 - t61 * t117;
t31 = t110 * t212 + t59 * t117;
t134 = g(1) * t29 - g(2) * t31 + g(3) * t63;
t133 = g(1) * t30 - g(2) * t32 + g(3) * t64;
t13 = -qJDD(4) * pkin(4) - t117 * t41 - t174;
t128 = -pkin(9) * t72 + t99 * t26;
t127 = t129 + t159;
t126 = qJDD(3) - t127;
t19 = t116 * t20;
t125 = -t16 * qJD(5) - t113 * t12 + t19;
t124 = pkin(9) * t215 + t13 - t134;
t3 = t23 * pkin(5) + qJDD(6) + t13;
t42 = t160 + t180 + (qJD(3) + t163) * qJD(2);
t122 = t144 * qJD(2) - t119 * t120 - t130 + t180 + t42;
t104 = qJDD(4) * t117;
t84 = t229 * t116;
t83 = t229 * t113;
t71 = t75 ^ 2;
t70 = t116 * t81;
t58 = t141 * t109;
t57 = (-qJDD(2) * t118 + t115 * t121) * t109;
t46 = t137 - t214;
t34 = t64 * qJD(4) - t117 * t167;
t33 = -t63 * qJD(4) + t114 * t167;
t28 = -t113 * t220 + t225;
t24 = t156 * t114 - t168 + t70;
t21 = t75 * pkin(5) + qJD(6) + t26;
t10 = t35 * qJD(5) + t113 * t166 + t33 * t116;
t9 = -t36 * qJD(5) - t33 * t113 + t116 * t166;
t2 = -t23 * qJ(6) - t75 * qJD(6) - t135;
t1 = t72 * pkin(5) + t22 * qJ(6) - t77 * qJD(6) + t125;
t5 = [t199, 0, -t57, -t58, t57, t58, t111 ^ 2 * qJDD(1) - g(3) + (t143 * qJD(2) + t115 * t42 - t118 * t46) * t109, 0, 0, 0, 0, 0, -t34 * qJD(4) - t63 * qJDD(4) + (t114 * t141 + t115 * t158) * t109, -t33 * qJD(4) - t64 * qJDD(4) + (t246 * t115 + t117 * t201) * t109, 0, 0, 0, 0, 0, t63 * t23 + t34 * t75 + t35 * t72 + t9 * t99, -t10 * t99 - t63 * t22 + t34 * t77 - t36 * t72, -t10 * t75 + t35 * t22 - t36 * t23 - t9 * t77, t1 * t35 + t8 * t10 + t2 * t36 + t21 * t34 + t3 * t63 + t4 * t9 - g(3); 0, qJDD(2), t127, t240, t126 - 0.2e1 * t214, 0.2e1 * qJD(2) * qJD(3) + 0.2e1 * t180 - t240, -t46 * pkin(2) - g(1) * t153 - g(2) * t154 + t42 * qJ(3) + t80 * qJD(3) - t143 * t196 - t171, t107 * qJDD(2) - 0.2e1 * t114 * t158, -0.2e1 * t114 * t176 + 0.2e1 * t183 * t198, -t120 * t114 + t104, -t120 * t117 - t179, 0, t122 * t114 - t117 * t239, t114 * t239 + t122 * t117, -t77 * t165 + (-t117 * t22 - t191 * t77) * t116 (t113 * t77 + t116 * t75) * t191 + (t221 - t116 * t23 + (t113 * t75 - t223) * qJD(5)) * t117 (-t22 - t169) * t114 + (t138 + t217) * t117 (-t23 + t170) * t114 + (-t139 - t218) * t117, t72 * t114 + t190 * t99, t70 * t72 + t244 * t99 + (-t215 * t81 + t129) * t113 + (t119 * t218 + t19 + (-qJD(4) * t26 - qJD(5) * t45 - t119 * t72 - t12) * t113 - t238 * t116) * t114 + (t75 * t164 + t26 * t188 + t13 * t113 - t119 * t23 + (-t204 * t99 + t15) * qJD(4)) * t117, -t225 * t72 + t243 * t99 + t129 * t116 + ((-t26 * t116 + t119 * t77) * qJD(4) + t238 * t113 + t175) * t114 + (-t16 * qJD(4) + t13 * t116 + t119 * t22 + t164 * t77 - t189 * t26) * t117, t24 * t22 - t28 * t23 - t234 * t77 - t233 * t75 + t146 * t191 + (qJD(5) * t145 - t1 * t116 - t113 * t2 + t130) * t117, t2 * t28 + t1 * t24 - g(1) * (t213 * t62 - t245 * t61 + t153) - g(2) * (t213 * t60 - t245 * t59 + t154) - t171 + t233 * t8 + t234 * t4 - t21 * t155 * t191 + (pkin(5) * t188 * t21 + t148 * t229 + t155 * t3) * t117 + (-g(3) * t245 * t118 + (t21 * t195 - g(3) * (-t117 * t229 + t213)) * t115) * t109; 0, 0, 0, 0, qJDD(2), -t121, t126 - t200 + t87 - t214, 0, 0, 0, 0, 0, t114 * t197 + t104, t117 * t197 - t179, 0, 0, 0, 0, 0, -t116 * t219 + (-t23 - t170) * t117 + (-t139 + t218) * t114, t113 * t219 + (t22 - t169) * t117 + (-t138 + t217) * t114 (-t75 * t190 + qJD(2) * t77 + (-t23 + t216) * t114) * t116 + (t77 * t190 + qJD(2) * t75 + (qJD(5) * t75 - t22) * t114) * t113, -t146 * qJD(2) + (-qJD(4) * t145 - t3) * t117 + (qJD(4) * t21 - qJD(5) * t146 - t1 * t113 + t2 * t116) * t114 - t129; 0, 0, 0, 0, 0, 0, 0, t117 * t121 * t114, -t198 * t121, t176, -t177, qJDD(4), t39 * qJD(4) - t117 * t150 + t134 + t174, t150 * t114 + t133 - t157, t223 * t99 - t221 (-t22 - t231) * t116 + (-t23 - t230) * t113 (-t117 * t77 + t203 * t99) * qJD(2) + t139 (-t114 * t224 + t117 * t75) * qJD(2) + t138, -t99 * t192, -t15 * t192 - pkin(4) * t23 - t39 * t75 - t67 * t99 + (t242 * t99 + t128) * t113 - t124 * t116, pkin(4) * t22 + t113 * t124 + t116 * t128 + t16 * t192 + t228 * t99 - t39 * t77, -t83 * t22 - t84 * t23 - t226 * t77 - t227 * t75 + (-t4 * t99 + t2) * t116 + (-t8 * t99 - t1) * t113 - t133, t2 * t84 - t1 * t83 - t3 * t101 - g(1) * (-t29 * t101 + t229 * t30) - g(2) * (t31 * t101 - t229 * t32) - g(3) * (-t63 * t101 + t229 * t64) + t227 * t8 + t226 * t4 + (pkin(5) * t224 - t39) * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77 * t75, -t71 + t237, -t22 + t231, t230 - t23, t72, t16 * t99 - t26 * t77 + t125 + t241, t15 * t99 + t26 * t75 - g(1) * (-t62 * t113 - t30 * t116) - g(2) * (-t60 * t113 + t32 * t116) + g(3) * t36 + t135, pkin(5) * t22 - t236 * t75, t236 * t8 + (-t21 * t77 + t1 + t241) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t71 - t237, t4 * t77 + t8 * t75 - t134 + t3;];
tau_reg  = t5;
