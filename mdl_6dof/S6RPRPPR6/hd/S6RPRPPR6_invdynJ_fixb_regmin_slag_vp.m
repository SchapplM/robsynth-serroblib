% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPRPPR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4,theta5]';
% 
% Output:
% tau_reg [6x26]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPRPPR6_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR6_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR6_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPPR6_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR6_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR6_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:54:41
% EndTime: 2019-03-09 02:54:49
% DurationCPUTime: 3.50s
% Computational Cost: add. (4363->402), mult. (8860->534), div. (0->0), fcn. (6324->14), ass. (0->202)
t181 = sin(qJ(3));
t184 = cos(qJ(3));
t264 = sin(pkin(9));
t265 = cos(pkin(9));
t123 = -t264 * t181 + t265 * t184;
t171 = qJ(3) + pkin(9);
t160 = sin(t171);
t162 = cos(t171);
t182 = sin(qJ(1));
t185 = cos(qJ(1));
t283 = g(1) * t182 - g(2) * t185;
t194 = -g(3) * t160 + t162 * t283;
t186 = -pkin(1) - pkin(7);
t134 = qJDD(1) * t186 + qJDD(2);
t127 = t184 * t134;
t135 = qJD(1) * t186 + qJD(2);
t233 = t184 * qJDD(1);
t236 = qJD(1) * qJD(4);
t237 = qJD(1) * qJD(3);
t243 = qJD(3) * t181;
t56 = -t184 * t236 - t135 * t243 + qJDD(3) * pkin(3) + t127 + (t181 * t237 - t233) * qJ(4);
t242 = qJD(3) * t184;
t71 = (-qJ(4) * qJD(1) + t135) * t242 + (-qJ(4) * qJDD(1) + t134 - t236) * t181;
t28 = -t264 * t71 + t265 * t56;
t25 = -qJDD(3) * pkin(4) + qJDD(5) - t28;
t285 = t25 + t194;
t195 = t181 * t265 + t184 * t264;
t282 = t195 * qJD(1);
t110 = qJD(6) + t282;
t180 = sin(qJ(6));
t183 = cos(qJ(6));
t116 = t123 * qJD(1);
t177 = sin(pkin(10));
t178 = cos(pkin(10));
t91 = -t178 * qJD(3) + t177 * t116;
t93 = t177 * qJD(3) + t178 * t116;
t43 = t180 * t93 + t183 * t91;
t288 = t43 * t110;
t205 = t180 * t91 - t183 * t93;
t287 = t110 * t205;
t223 = qJD(3) * t264;
t224 = qJD(3) * t265;
t115 = -t181 * t224 - t184 * t223;
t245 = qJD(1) * t181;
t108 = -qJ(4) * t245 + t181 * t135;
t101 = t264 * t108;
t244 = qJD(1) * t184;
t109 = -qJ(4) * t244 + t184 * t135;
t104 = qJD(3) * pkin(3) + t109;
t62 = t104 * t265 - t101;
t48 = -qJD(3) * pkin(4) + qJD(5) - t62;
t192 = t48 * t115 + t25 * t123 + t283;
t188 = qJD(1) ^ 2;
t196 = -t188 * qJ(2) - t283;
t126 = t183 * t177 + t180 * t178;
t119 = t126 * qJD(6);
t125 = t180 * t177 - t183 * t178;
t221 = qJDD(1) * t264;
t222 = qJDD(1) * t265;
t83 = -t116 * qJD(3) - t181 * t222 - t184 * t221;
t81 = -qJDD(6) + t83;
t286 = -t119 * t110 + t125 * t81;
t216 = g(1) * t185 + g(2) * t182;
t284 = t160 * t216;
t281 = -qJD(6) + t110;
t172 = qJDD(1) * qJ(2);
t173 = qJD(1) * qJD(2);
t232 = 0.2e1 * t173;
t280 = 0.2e1 * t172 + t232 - t216;
t118 = t125 * qJD(6);
t267 = -t125 * t282 - t118;
t270 = t126 * t81;
t279 = -t110 * t267 + t270;
t84 = qJD(3) * t282 + t181 * t221 - t184 * t222;
t72 = -t178 * qJDD(3) - t177 * t84;
t73 = t177 * qJDD(3) - t178 * t84;
t9 = -qJD(6) * t205 + t180 * t73 + t183 * t72;
t111 = t282 ^ 2;
t277 = pkin(8) * t178;
t275 = g(3) * t162;
t274 = g(3) * t181;
t166 = t181 * pkin(3);
t146 = pkin(3) * t264 + qJ(5);
t273 = pkin(8) + t146;
t230 = t184 * t237;
t234 = t181 * qJDD(1);
t204 = qJDD(4) + t172 + t173 + (t230 + t234) * pkin(3);
t19 = -t83 * pkin(4) + t84 * qJ(5) - t116 * qJD(5) + t204;
t29 = t264 * t56 + t265 * t71;
t23 = qJDD(3) * qJ(5) + qJD(3) * qJD(5) + t29;
t7 = t177 * t19 + t178 * t23;
t227 = t265 * t108;
t63 = t264 * t104 + t227;
t49 = qJD(3) * qJ(5) + t63;
t129 = pkin(3) * t245 + qJD(1) * qJ(2) + qJD(4);
t64 = pkin(4) * t282 - t116 * qJ(5) + t129;
t27 = t177 * t64 + t178 * t49;
t114 = t181 * t223 - t184 * t224;
t238 = pkin(3) * t242 + qJD(2);
t42 = -t114 * pkin(4) - t115 * qJ(5) - t123 * qJD(5) + t238;
t251 = qJ(4) - t186;
t106 = -t184 * qJD(4) + t243 * t251;
t131 = t251 * t184;
t107 = -qJD(3) * t131 - t181 * qJD(4);
t70 = t106 * t264 + t107 * t265;
t21 = t177 * t42 + t178 * t70;
t75 = t109 * t265 - t101;
t76 = pkin(3) * t244 + t116 * pkin(4) + qJ(5) * t282;
t31 = t177 * t76 + t178 * t75;
t252 = qJ(2) + t166;
t82 = pkin(4) * t195 - t123 * qJ(5) + t252;
t130 = t251 * t181;
t89 = -t130 * t265 - t131 * t264;
t37 = t177 * t82 + t178 * t89;
t272 = t116 * t43;
t271 = t116 * t205;
t269 = t177 * t83;
t268 = t178 * t83;
t67 = t126 * t282;
t266 = t119 + t67;
t263 = pkin(1) * qJDD(1);
t262 = t110 * t114;
t261 = t282 * t177;
t260 = t115 * t177;
t259 = t123 * t177;
t258 = t123 * t178;
t170 = pkin(10) + qJ(6);
t159 = sin(t170);
t257 = t182 * t159;
t161 = cos(t170);
t256 = t182 * t161;
t255 = t185 * t159;
t254 = t185 * t161;
t250 = -qJD(5) + t48;
t249 = t185 * pkin(1) + t182 * qJ(2);
t176 = t184 ^ 2;
t247 = t181 ^ 2 - t176;
t187 = qJD(3) ^ 2;
t246 = -t187 - t188;
t241 = qJD(6) * t180;
t240 = qJD(6) * t183;
t239 = t129 * qJD(1);
t235 = qJDD(3) * t181;
t6 = -t177 * t23 + t178 * t19;
t2 = -t83 * pkin(5) - t73 * pkin(8) + t6;
t5 = -pkin(8) * t72 + t7;
t231 = -t180 * t5 + t183 * t2;
t229 = -t182 * pkin(1) + t185 * qJ(2);
t20 = -t177 * t70 + t178 * t42;
t26 = -t177 * t49 + t178 * t64;
t30 = -t177 * t75 + t178 * t76;
t36 = -t177 * t89 + t178 * t82;
t220 = qJDD(2) - t263;
t219 = -t110 * t67 + t286;
t151 = -pkin(3) * t265 - pkin(4);
t214 = -t7 * t177 - t6 * t178;
t213 = t180 * t2 + t183 * t5;
t212 = -t26 * t114 + t195 * t6;
t211 = t27 * t114 - t195 * t7;
t69 = -t265 * t106 + t107 * t264;
t210 = t160 * pkin(4) - t162 * qJ(5);
t12 = pkin(5) * t282 - t93 * pkin(8) + t26;
t15 = -pkin(8) * t91 + t27;
t3 = t183 * t12 - t180 * t15;
t4 = t180 * t12 + t183 * t15;
t208 = -t26 * t177 + t27 * t178;
t24 = pkin(5) * t195 - pkin(8) * t258 + t36;
t32 = -pkin(8) * t259 + t37;
t207 = t180 * t32 - t183 * t24;
t206 = t180 * t24 + t183 * t32;
t74 = t109 * t264 + t227;
t88 = -t130 * t264 + t265 * t131;
t179 = -qJ(4) - pkin(7);
t203 = t185 * t166 + t182 * t179 + t229;
t202 = t182 * t166 - t185 * t179 + t249;
t8 = -t180 * t72 + t183 * t73 - t91 * t240 - t241 * t93;
t120 = t273 * t177;
t201 = pkin(8) * t261 - t178 * qJD(5) + qJD(6) * t120 + t31;
t121 = t273 * t178;
t200 = t116 * pkin(5) + t177 * qJD(5) + qJD(6) * t121 + t277 * t282 + t30;
t199 = t110 * t126;
t198 = 0.2e1 * qJ(2) * t237 + qJDD(3) * t186;
t190 = -t63 * t114 + t62 * t115 + t28 * t123 + t195 * t29 - t283;
t189 = -t186 * t187 + t280;
t163 = qJDD(3) * t184;
t132 = -t178 * pkin(5) + t151;
t100 = t160 * t254 - t257;
t99 = t160 * t255 + t256;
t98 = t160 * t256 + t255;
t97 = -t160 * t257 + t254;
t78 = t125 * t123;
t77 = t126 * t123;
t57 = pkin(5) * t259 + t88;
t39 = -pkin(5) * t261 + t74;
t38 = pkin(5) * t260 + t69;
t35 = t91 * pkin(5) + t48;
t34 = t115 * t126 + t240 * t258 - t241 * t259;
t33 = -t115 * t125 - t123 * t119;
t13 = -pkin(8) * t260 + t21;
t11 = -t114 * pkin(5) - t115 * t277 + t20;
t10 = t72 * pkin(5) + t25;
t1 = [qJDD(1), t283, t216, qJDD(2) - t283 - 0.2e1 * t263, t280, -t220 * pkin(1) - g(1) * t229 - g(2) * t249 + (t232 + t172) * qJ(2), t176 * qJDD(1) - 0.2e1 * t181 * t230, -0.2e1 * t181 * t233 + 0.2e1 * t237 * t247, -t187 * t181 + t163, -t187 * t184 - t235, 0, t181 * t189 + t184 * t198, -t181 * t198 + t184 * t189, t69 * t116 - t282 * t70 + t89 * t83 - t88 * t84 - t190, -g(1) * t203 - g(2) * t202 + t129 * t238 + t204 * t252 - t28 * t88 + t29 * t89 - t62 * t69 + t63 * t70, t177 * t192 - t178 * t284 + t20 * t282 - t36 * t83 + t69 * t91 + t88 * t72 + t212, t177 * t284 + t178 * t192 - t21 * t282 + t37 * t83 + t69 * t93 + t88 * t73 + t211, -t20 * t93 - t21 * t91 - t36 * t73 - t37 * t72 + t216 * t162 + t214 * t123 + (-t177 * t27 - t178 * t26) * t115, t7 * t37 + t27 * t21 + t6 * t36 + t26 * t20 + t25 * t88 + t48 * t69 - g(1) * (t185 * t210 + t203) - g(2) * (t182 * t210 + t202) -t205 * t33 - t78 * t8, t205 * t34 - t33 * t43 - t77 * t8 + t78 * t9, t110 * t33 + t114 * t205 + t195 * t8 + t78 * t81, -t110 * t34 + t114 * t43 - t195 * t9 + t77 * t81, -t195 * t81 - t262 (t183 * t11 - t180 * t13) * t110 + t207 * t81 + t231 * t195 - t3 * t114 + t38 * t43 + t57 * t9 + t10 * t77 + t35 * t34 - g(1) * t100 - g(2) * t98 + (-t110 * t206 - t195 * t4) * qJD(6) -(t180 * t11 + t183 * t13) * t110 + t206 * t81 - t213 * t195 + t4 * t114 - t38 * t205 + t57 * t8 - t10 * t78 + t35 * t33 + g(1) * t99 - g(2) * t97 + (t110 * t207 - t195 * t3) * qJD(6); 0, 0, 0, qJDD(1), -t188, t220 + t196, 0, 0, 0, 0, 0, t181 * t246 + t163, t184 * t246 - t235, t114 * t282 - t115 * t116 + t123 * t84 + t195 * t83, t190 - t239, t195 * t269 - t115 * t91 - t123 * t72 + (-qJD(1) * t178 + t114 * t177) * t282, t195 * t268 - t115 * t93 - t123 * t73 + (qJD(1) * t177 + t114 * t178) * t282 (qJD(1) * t93 + t114 * t91 - t195 * t72) * t178 + (qJD(1) * t91 - t114 * t93 + t195 * t73) * t177 (-qJD(1) * t26 - t211) * t178 + (-qJD(1) * t27 - t212) * t177 - t192, 0, 0, 0, 0, 0, -t115 * t43 - t123 * t9 + t114 * t199 + t125 * t110 * qJD(1) - (-t110 * t118 - t270) * t195, qJD(1) * t199 + t115 * t205 - t123 * t8 - t125 * t262 - t195 * t286; 0, 0, 0, 0, 0, 0, t184 * t188 * t181, -t247 * t188, t233, -t234, qJDD(3), t184 * t196 + t127 + t274, g(3) * t184 + (-t134 - t196) * t181 (t63 - t74) * t116 + (-t62 + t75) * t282 + (t264 * t83 + t265 * t84) * pkin(3), t62 * t74 - t63 * t75 + (t264 * t29 + t265 * t28 + t274 + (-t283 - t239) * t184) * pkin(3), t146 * t269 - t26 * t116 + t151 * t72 - t74 * t91 + (t177 * t250 - t30) * t282 - t285 * t178, t146 * t268 + t27 * t116 + t151 * t73 - t74 * t93 + (t178 * t250 + t31) * t282 + t285 * t177, -t275 + t30 * t93 + t31 * t91 - t283 * t160 + (-qJD(5) * t91 - t146 * t72 - t26 * t282 + t7) * t178 + (qJD(5) * t93 + t146 * t73 - t27 * t282 - t6) * t177, t25 * t151 - t27 * t31 - t26 * t30 - t48 * t74 - g(3) * (-t166 - t210) + (-t6 * t177 + t7 * t178) * t146 + t208 * qJD(5) - t283 * (pkin(3) * t184 + pkin(4) * t162 + qJ(5) * t160) t126 * t8 - t205 * t267, -t125 * t8 - t126 * t9 + t205 * t266 - t267 * t43, t271 - t279, t219 + t272, -t110 * t116 -(-t183 * t120 - t180 * t121) * t81 + t132 * t9 + t10 * t125 - t3 * t116 - t39 * t43 + t266 * t35 + (t180 * t201 - t183 * t200) * t110 - t194 * t161 (-t180 * t120 + t183 * t121) * t81 + t132 * t8 + t10 * t126 + t4 * t116 + t39 * t205 + t267 * t35 + (t180 * t200 + t183 * t201) * t110 + t194 * t159; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t116 ^ 2 - t111, t62 * t116 + t282 * t63 + t204 - t216, -t111 * t177 - t116 * t91 - t268, -t178 * t111 - t116 * t93 + t269, -t177 * t72 - t178 * t73 + (t177 * t93 - t178 * t91) * t282, -t48 * t116 + t208 * t282 - t214 - t216, 0, 0, 0, 0, 0, t219 - t272, t271 + t279; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t282 * t93 + t72, -t282 * t91 + t73, -t91 ^ 2 - t93 ^ 2, t26 * t93 + t27 * t91 + t285, 0, 0, 0, 0, 0, t9 - t287, t8 - t288; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t205 * t43, t205 ^ 2 - t43 ^ 2, t8 + t288, -t9 - t287, -t81, -g(1) * t97 - g(2) * t99 + t159 * t275 + t205 * t35 + t281 * t4 + t231, g(1) * t98 - g(2) * t100 + t161 * t275 + t281 * t3 + t35 * t43 - t213;];
tau_reg  = t1;
