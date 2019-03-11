% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RRPPRP2
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
% 
% Output:
% tau_reg [6x25]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRPPRP2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP2_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP2_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRP2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP2_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:32:00
% EndTime: 2019-03-09 08:32:08
% DurationCPUTime: 3.53s
% Computational Cost: add. (4213->413), mult. (9558->500), div. (0->0), fcn. (6778->10), ass. (0->216)
t150 = sin(pkin(9));
t151 = cos(pkin(9));
t155 = sin(qJ(2));
t158 = cos(qJ(2));
t111 = t150 * t158 + t151 * t155;
t101 = t111 * qJD(1);
t280 = qJD(5) + t101;
t154 = sin(qJ(5));
t157 = cos(qJ(5));
t225 = qJD(1) * t155;
t236 = t151 * t158;
t98 = -qJD(1) * t236 + t150 * t225;
t79 = qJD(2) * t154 - t157 * t98;
t201 = t280 * t79;
t100 = t111 * qJD(2);
t215 = t158 * qJDD(1);
t124 = t151 * t215;
t216 = t155 * qJDD(1);
t186 = -t150 * t216 + t124;
t164 = qJD(1) * t100 - t186;
t221 = qJD(5) * t157;
t222 = qJD(5) * t154;
t30 = qJD(2) * t222 - t157 * qJDD(2) - t154 * t164 - t98 * t221;
t281 = t30 - t201;
t144 = qJ(2) + pkin(9);
t140 = cos(t144);
t132 = g(3) * t140;
t139 = sin(t144);
t156 = sin(qJ(1));
t159 = cos(qJ(1));
t272 = -g(1) * t159 - g(2) * t156;
t174 = t139 * t272 + t132;
t217 = qJD(1) * qJD(2);
t208 = t155 * t217;
t214 = pkin(2) * t208 + qJDD(3);
t207 = t158 * t217;
t72 = qJDD(1) * t111 - t150 * t208 + t151 * t207;
t188 = -t72 * qJ(4) + t214;
t219 = t101 * qJD(4);
t238 = t150 * t155;
t264 = pkin(3) + pkin(8);
t142 = t158 * pkin(2);
t136 = t142 + pkin(1);
t279 = t136 * qJDD(1);
t11 = -t279 + t188 - t219 + (t238 * qJDD(1) + t217 * t111 - t124) * t264;
t153 = -qJ(3) - pkin(7);
t116 = t153 * t155;
t204 = qJD(2) * t153;
t96 = -qJD(3) * t155 + t158 * t204;
t67 = qJDD(2) * pkin(2) + qJD(1) * t96 + qJDD(1) * t116;
t117 = t153 * t158;
t95 = qJD(3) * t158 + t155 * t204;
t73 = qJD(1) * t95 - qJDD(1) * t117;
t26 = -t150 * t73 + t151 * t67;
t191 = qJDD(4) - t26;
t18 = pkin(4) * t72 - t264 * qJDD(2) + t191;
t114 = qJD(1) * t117;
t104 = t150 * t114;
t113 = qJD(1) * t116;
t109 = qJD(2) * pkin(2) + t113;
t70 = t109 * t151 + t104;
t187 = qJD(4) - t70;
t257 = pkin(4) * t101;
t38 = -t264 * qJD(2) + t187 + t257;
t213 = -t157 * t11 - t154 * t18 - t38 * t221;
t115 = -qJD(1) * t136 + qJD(3);
t171 = -qJ(4) * t101 + t115;
t33 = t264 * t98 + t171;
t179 = -t33 * t222 - t213;
t197 = qJDD(2) * t154 - t157 * t164;
t81 = qJD(2) * t157 + t154 * t98;
t31 = qJD(5) * t81 + t197;
t2 = -qJ(6) * t31 - qJD(6) * t79 + t179;
t14 = -t154 * t33 + t157 * t38;
t7 = -qJ(6) * t81 + t14;
t6 = pkin(5) * t280 + t7;
t268 = t280 * t6 - t2;
t15 = t154 * t38 + t157 * t33;
t205 = -t11 * t154 + t157 * t18;
t166 = -qJD(5) * t15 + t205;
t68 = qJDD(5) + t72;
t1 = pkin(5) * t68 + qJ(6) * t30 - qJD(6) * t81 + t166;
t8 = -qJ(6) * t79 + t15;
t269 = t280 * t8 + t1;
t283 = -t154 * t268 + t157 * t269 + t174;
t276 = t157 * t280;
t282 = t81 * t276;
t200 = t154 * t280;
t61 = t157 * t68;
t181 = -t200 * t280 + t61;
t97 = t101 ^ 2;
t277 = -t98 ^ 2 - t97;
t273 = g(1) * t156 - g(2) * t159;
t275 = t273 * t140;
t75 = t113 * t151 + t104;
t229 = -qJD(4) + t75;
t129 = t139 * qJ(4);
t274 = -t140 * pkin(3) - t129;
t271 = qJD(2) * t101;
t231 = t157 * t159;
t234 = t154 * t156;
t91 = t139 * t231 - t234;
t232 = t156 * t157;
t233 = t154 * t159;
t93 = t139 * t232 + t233;
t270 = -g(1) * t91 - g(2) * t93 + t157 * t132;
t134 = -pkin(2) * t151 - pkin(3);
t127 = -pkin(8) + t134;
t262 = pkin(4) * t98;
t237 = t151 * t114;
t71 = t150 * t109 - t237;
t62 = -qJD(2) * qJ(4) - t71;
t43 = -t62 - t262;
t267 = t127 * t68 + t280 * t43;
t266 = t81 ^ 2;
t263 = -t7 + t6;
t220 = qJD(6) * t157;
t230 = qJ(6) - t127;
t203 = pkin(2) * t225 + qJ(4) * t98;
t39 = t264 * t101 + t203;
t74 = t113 * t150 - t237;
t48 = t74 - t262;
t46 = t157 * t48;
t259 = t230 * t222 - t220 + pkin(5) * t98 - t46 - (-qJ(6) * t101 - t39) * t154;
t258 = pkin(2) * t155;
t256 = pkin(5) * t154;
t252 = g(3) * t158;
t251 = t79 * t98;
t250 = t81 * t98;
t249 = t154 * t48 + t157 * t39;
t110 = -t236 + t238;
t195 = -qJ(4) * t111 - t136;
t47 = t264 * t110 + t195;
t76 = -t151 * t116 - t117 * t150;
t57 = pkin(4) * t111 + t76;
t248 = t154 * t57 + t157 * t47;
t27 = t150 * t67 + t151 * t73;
t108 = t230 * t157;
t243 = qJ(6) * t157;
t247 = -qJD(5) * t108 - qJD(6) * t154 - t101 * t243 - t249;
t245 = t154 * t68;
t244 = t157 * t30;
t242 = qJDD(2) * pkin(3);
t241 = t100 * t154;
t240 = t100 * t157;
t239 = t140 * t159;
t235 = t153 * t159;
t228 = t257 - t229;
t135 = pkin(5) * t157 + pkin(4);
t227 = t135 - t153;
t148 = t155 ^ 2;
t226 = -t158 ^ 2 + t148;
t224 = qJD(2) * t155;
t223 = qJD(5) * t127;
t211 = qJDD(2) * qJ(4) + t27;
t138 = pkin(2) * t224;
t210 = t142 - t274;
t130 = pkin(2) * t150 + qJ(4);
t55 = t150 * t95 - t151 * t96;
t202 = -qJ(6) * t110 - t47;
t146 = qJD(2) * qJD(4);
t23 = -t146 - t211;
t123 = t159 * t136;
t196 = g(2) * (pkin(3) * t239 + t159 * t129 + t123);
t194 = -pkin(3) * t139 - t258;
t56 = t150 * t96 + t151 * t95;
t77 = t116 * t150 - t117 * t151;
t152 = -qJ(6) - pkin(8);
t185 = t139 * t256 - t140 * t152;
t183 = -t151 * t208 + t124;
t182 = -t136 + t274;
t180 = -0.2e1 * pkin(1) * t217 - pkin(7) * qJDD(2);
t103 = qJD(2) * t236 - t150 * t224;
t176 = -qJ(4) * t103 - qJD(4) * t111 + t138;
t28 = t264 * t100 + t176;
t40 = pkin(4) * t103 + t55;
t178 = t154 * t40 + t157 * t28 + t57 * t221 - t47 * t222;
t177 = t110 * t222 - t240;
t41 = -pkin(4) * t100 + t56;
t175 = -t276 * t280 - t245;
t173 = t214 - t279;
t172 = -g(3) * t139 + t272 * t140;
t160 = qJD(2) ^ 2;
t169 = 0.2e1 * qJDD(1) * pkin(1) - pkin(7) * t160 + t273;
t161 = qJD(1) ^ 2;
t168 = pkin(1) * t161 - pkin(7) * qJDD(1) - t272;
t19 = -pkin(4) * t164 - t23;
t167 = t19 + t172;
t165 = t55 * t101 - t56 * t98 + t76 * t72 + t272;
t51 = pkin(3) * t98 + t171;
t163 = t101 * t51 + t174 + t191;
t5 = t31 * pkin(5) + qJDD(6) + t19;
t162 = (pkin(3) * t238 - t136) * qJDD(1) - (-t150 * t207 + t183) * pkin(3) + t188;
t120 = qJ(4) * t239;
t118 = t156 * t140 * qJ(4);
t107 = t230 * t154;
t94 = -t139 * t234 + t231;
t92 = t139 * t233 + t232;
t89 = qJD(2) * t98;
t78 = t79 ^ 2;
t69 = pkin(3) * t110 + t195;
t59 = -qJD(2) * pkin(3) + t187;
t58 = -pkin(4) * t110 + t77;
t54 = pkin(3) * t101 + t203;
t53 = t157 * t57;
t42 = pkin(3) * t100 + t176;
t37 = t157 * t40;
t29 = t157 * t31;
t24 = t191 - t242;
t22 = pkin(5) * t79 + qJD(6) + t43;
t21 = t162 - t219;
t20 = t110 * t243 + t248;
t13 = pkin(5) * t111 + t154 * t202 + t53;
t4 = -qJ(6) * t177 + t110 * t220 + t178;
t3 = pkin(5) * t103 + t37 + t202 * t221 + (-qJ(6) * t100 - qJD(5) * t57 - qJD(6) * t110 - t28) * t154;
t9 = [qJDD(1), t273, -t272, qJDD(1) * t148 + 0.2e1 * t155 * t207, 0.2e1 * t155 * t215 - 0.2e1 * t226 * t217, qJDD(2) * t155 + t158 * t160, qJDD(2) * t158 - t155 * t160, 0, t155 * t180 + t169 * t158, -t155 * t169 + t180 * t158, t77 * (t186 - t271) - t27 * t110 - t71 * t100 - t26 * t111 - t70 * t103 + t165, t27 * t77 + t71 * t56 - t26 * t76 - t70 * t55 - t173 * t136 + t115 * t138 - g(1) * (-t136 * t156 - t235) - g(2) * (-t153 * t156 + t123) t62 * t100 + t59 * t103 + t23 * t110 + t24 * t111 - t164 * t77 + t165, -t42 * t98 + t69 * t186 - t21 * t110 - t51 * t100 + t76 * qJDD(2) - t275 + (-t101 * t69 + t55) * qJD(2), qJD(2) * t56 + qJDD(2) * t77 - t101 * t42 - t103 * t51 - t111 * t21 + t139 * t273 - t69 * t72, t21 * t69 + t51 * t42 - t23 * t77 - t62 * t56 + t24 * t76 + t59 * t55 + g(1) * t235 - t196 + (-g(1) * t182 + g(2) * t153) * t156, t81 * t241 + (-t154 * t30 + t81 * t221) * t110 (-t154 * t79 + t157 * t81) * t100 + (-t154 * t31 - t244 + (-t154 * t81 - t157 * t79) * qJD(5)) * t110, t280 * t241 + t103 * t81 - t111 * t30 + (t221 * t280 + t245) * t110, t280 * t240 - t103 * t79 - t111 * t31 + (-t222 * t280 + t61) * t110, t103 * t280 + t111 * t68 (-t154 * t28 + t37) * t280 + (-t154 * t47 + t53) * t68 + t205 * t111 + t14 * t103 + t41 * t79 + t58 * t31 - g(1) * t94 - g(2) * t92 + (-t43 * t100 - t19 * t110) * t157 + (t43 * t154 * t110 - t15 * t111 - t248 * t280) * qJD(5), -t178 * t280 - t248 * t68 - t179 * t111 - t15 * t103 + t41 * t81 - t58 * t30 + t43 * t241 + g(1) * t93 - g(2) * t91 + (t19 * t154 + t221 * t43) * t110, t13 * t30 - t20 * t31 - t3 * t81 - t4 * t79 + t275 + (-t154 * t6 + t157 * t8) * t100 + (-t1 * t154 + t157 * t2 + (-t154 * t8 - t157 * t6) * qJD(5)) * t110, t2 * t20 + t8 * t4 + t1 * t13 + t6 * t3 + t5 * (-t110 * t135 + t77) + t22 * (pkin(5) * t177 + t41) - t196 + (-g(1) * t227 - g(2) * t185) * t159 + (-g(1) * (t182 - t185) - g(2) * t227) * t156; 0, 0, 0, -t155 * t161 * t158, t226 * t161, t216, t215, qJDD(2), t155 * t168 - t252, g(3) * t155 + t168 * t158 (t75 - t70) * t98 + (t71 - t74) * t101 + (-t151 * t72 + ((-t207 - t216) * t150 + t183) * t150) * pkin(2), t70 * t74 - t71 * t75 + (-t252 + t150 * t27 + t151 * t26 + (-qJD(1) * t115 - t272) * t155) * pkin(2), -t130 * t164 + t134 * t72 + (-t62 - t74) * t101 + (t59 + t229) * t98, -qJD(2) * t74 + t54 * t98 + (-pkin(3) + t134) * qJDD(2) + t163, -qJD(2) * t75 + qJDD(2) * t130 + t101 * t54 - t51 * t98 + 0.2e1 * t146 + t172 + t211, -t23 * t130 + t24 * t134 - t51 * t54 - t59 * t74 - g(1) * (t159 * t194 + t120) - g(2) * (t156 * t194 + t118) - g(3) * t210 + t229 * t62, -t200 * t81 - t244, -t29 - t282 + (t30 + t201) * t154, t181 + t250, t175 - t251, t280 * t98, t130 * t31 + t14 * t98 - t46 * t280 + t228 * t79 + t267 * t157 + ((t39 - t223) * t280 + t167) * t154, -t130 * t30 + t249 * t280 - t15 * t98 + t228 * t81 - t267 * t154 + (-t223 * t280 + t167) * t157, t107 * t31 - t108 * t30 - t247 * t79 - t259 * t81 - t283, -t2 * t107 - t1 * t108 + t5 * (t130 + t256) - g(1) * t120 - g(2) * t118 - g(3) * (t185 + t210) + t247 * t8 + t259 * t6 + (pkin(5) * t276 + t228) * t22 + t272 * (t140 * t256 - t258 + (-pkin(3) + t152) * t139); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t277, t101 * t70 + t71 * t98 + t173 - t273, t277, t186 - 0.2e1 * t271, -t72 + t89, -t62 * t98 + (-qJD(4) - t59) * t101 + t162 - t273, 0, 0, 0, 0, 0, t175 + t251, t250 - t181, -t281 * t154 + t282 - t29, -t154 * t269 - t157 * t268 + t22 * t98 - t273; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t89 + t72, -t101 * t98 + qJDD(2), -t97 - t160, qJD(2) * t62 + t163 - t242, 0, 0, 0, 0, 0, -qJD(2) * t79 + t181, -qJD(2) * t81 + t175, t281 * t157 + (t280 * t81 - t31) * t154, -qJD(2) * t22 + t283; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81 * t79, -t78 + t266, -t281, -t197 + (-qJD(5) + t280) * t81, t68, t15 * t280 - t43 * t81 + t166 + t270, g(1) * t92 - g(2) * t94 + t14 * t280 + t43 * t79 + (qJD(5) * t33 - t132) * t154 + t213, pkin(5) * t30 - t263 * t79, t263 * t8 + (-t22 * t81 + t1 + t270) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t78 - t266, t6 * t81 + t8 * t79 + t172 + t5;];
tau_reg  = t9;
