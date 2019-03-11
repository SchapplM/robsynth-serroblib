% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6RPPPRR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3,theta4]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPPPRR3_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR3_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR3_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPPRR3_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR3_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR3_invdynJ_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:34:11
% EndTime: 2019-03-09 01:34:13
% DurationCPUTime: 2.52s
% Computational Cost: add. (5649->404), mult. (10434->506), div. (0->0), fcn. (7107->12), ass. (0->212)
t153 = cos(pkin(10));
t269 = cos(qJ(5));
t215 = t269 * t153;
t111 = qJD(1) * t215;
t151 = sin(pkin(10));
t158 = sin(qJ(5));
t242 = t158 * t151;
t213 = qJD(1) * t242;
t90 = -t111 + t213;
t241 = qJD(6) - t90;
t157 = sin(qJ(6));
t159 = cos(qJ(6));
t132 = t153 * qJD(3);
t160 = -pkin(1) - pkin(2);
t107 = t160 * qJD(1) + qJD(2);
t152 = sin(pkin(9));
t154 = cos(pkin(9));
t234 = qJD(1) * t154;
t84 = qJ(2) * t234 + t152 * t107;
t72 = -qJD(1) * qJ(4) + t84;
t45 = t132 + (pkin(7) * qJD(1) - t72) * t151;
t235 = qJD(1) * t153;
t58 = t151 * qJD(3) + t153 * t72;
t46 = -pkin(7) * t235 + t58;
t20 = t158 * t45 + t269 * t46;
t18 = qJD(5) * pkin(8) + t20;
t236 = qJD(1) * t152;
t83 = -qJ(2) * t236 + t154 * t107;
t71 = qJD(1) * pkin(3) + qJD(4) - t83;
t63 = pkin(4) * t235 + t71;
t96 = t269 * t151 + t158 * t153;
t91 = t96 * qJD(1);
t26 = -t90 * pkin(5) + t91 * pkin(8) + t63;
t9 = -t157 * t18 + t159 * t26;
t284 = t241 * t9;
t176 = t215 - t242;
t93 = t96 * qJD(5);
t256 = t152 * t93 + t176 * t234;
t10 = t157 * t26 + t159 * t18;
t283 = t241 * t10;
t149 = pkin(10) + qJ(5);
t133 = sin(t149);
t134 = cos(t149);
t268 = sin(qJ(1));
t270 = cos(qJ(1));
t95 = -t268 * t152 - t270 * t154;
t273 = g(1) * t95;
t97 = t270 * t152 - t268 * t154;
t200 = g(2) * t97 + t273;
t169 = g(3) * t134 - t200 * t133;
t130 = t153 * qJDD(3);
t227 = qJD(1) * qJD(2);
t117 = t154 * t227;
t106 = t160 * qJDD(1) + qJDD(2);
t223 = t154 * qJDD(1);
t240 = qJ(2) * t223 + t152 * t106;
t65 = t117 + t240;
t56 = -qJDD(1) * qJ(4) - qJD(1) * qJD(4) + t65;
t34 = t130 + (pkin(7) * qJDD(1) - t56) * t151;
t224 = t153 * qJDD(1);
t38 = t151 * qJDD(3) + t153 * t56;
t35 = -pkin(7) * t224 + t38;
t211 = t158 * t35 - t269 * t34;
t8 = -t20 * qJD(5) - t211;
t6 = -qJDD(5) * pkin(5) - t8;
t167 = -t6 + t169;
t282 = -pkin(8) * qJD(6) * t241 + t167;
t207 = t159 * t241;
t210 = qJDD(1) * t269;
t226 = t151 * qJDD(1);
t187 = -t153 * t210 + t158 * t226;
t53 = qJD(1) * t93 + t187;
t50 = -qJDD(6) + t53;
t281 = -t157 * t50 + t207 * t241;
t272 = g(2) * t95;
t201 = g(1) * t97 - t272;
t279 = t201 * t133;
t212 = qJD(5) * t269;
t233 = qJD(5) * t158;
t278 = -t151 * t233 + t153 * t212;
t70 = t157 * qJD(5) - t159 * t91;
t244 = qJD(6) * t70;
t217 = -qJD(5) * t111 - t151 * t210 - t158 * t224;
t52 = qJD(5) * t213 + t217;
t25 = -t159 * qJDD(5) + t157 * t52 + t244;
t277 = qJDD(1) * pkin(3) + qJDD(4);
t125 = pkin(4) * t224;
t115 = t152 * t227;
t225 = t152 * qJDD(1);
t206 = -qJ(2) * t225 + t154 * t106;
t64 = -t115 + t206;
t59 = -t64 + t277;
t48 = t125 + t59;
t14 = -t53 * pkin(5) - t52 * pkin(8) + t48;
t13 = t159 * t14;
t222 = -t158 * t34 - t45 * t212 - t269 * t35;
t7 = -t46 * t233 - t222;
t5 = qJDD(5) * pkin(8) + t7;
t2 = -qJD(6) * t10 - t157 * t5 + t13;
t276 = t2 + t283;
t267 = g(3) * t133;
t170 = t200 * t134 + t267;
t194 = -t241 * t278 + t50 * t96;
t231 = qJD(6) * t159;
t218 = t96 * t231;
t275 = -t157 * t194 + t218 * t241;
t274 = t91 ^ 2;
t102 = t154 * qJ(2) + t152 * t160;
t98 = -qJ(4) + t102;
t271 = pkin(7) - t98;
t229 = t159 * qJD(5);
t68 = -t157 * t91 - t229;
t265 = t68 * t90;
t264 = t70 * t68;
t263 = t70 * t91;
t262 = t90 * t91;
t261 = t91 * t68;
t260 = -t157 * t25 - t68 * t231;
t82 = t176 * t152;
t60 = -t159 * t154 - t157 * t82;
t259 = qJD(6) * t60 - t157 * t236 - t256 * t159;
t183 = t157 * t154 - t159 * t82;
t258 = qJD(6) * t183 + t256 * t157 - t159 * t236;
t257 = -t152 * t278 + t154 * t91;
t255 = t134 * t95;
t254 = t134 * t97;
t252 = t157 * t68;
t251 = t157 * t70;
t250 = t158 * t46;
t249 = t159 * t68;
t248 = t159 * t70;
t232 = qJD(6) * t157;
t24 = -qJD(6) * t229 - t157 * qJDD(5) - t159 * t52 - t91 * t232;
t247 = t24 * t157;
t246 = t25 * t159;
t245 = pkin(1) * qJDD(1);
t161 = qJD(1) ^ 2;
t243 = t154 * t161;
t239 = t270 * pkin(1) + t268 * qJ(2);
t238 = g(1) * t268 - g(2) * t270;
t147 = t151 ^ 2;
t148 = t153 ^ 2;
t237 = t147 + t148;
t230 = t152 * qJD(2);
t228 = qJ(2) * qJDD(1);
t220 = 0.2e1 * t227;
t219 = t96 * t232;
t216 = t270 * pkin(2) + t239;
t208 = t157 * t241;
t101 = -t152 * qJ(2) + t154 * t160;
t205 = qJDD(1) * t237;
t203 = qJDD(2) - t245;
t99 = pkin(3) - t101;
t199 = -t268 * pkin(1) + t270 * qJ(2);
t1 = qJD(6) * t9 + t157 * t14 + t159 * t5;
t198 = t1 - t284;
t197 = t134 * pkin(5) + t133 * pkin(8);
t196 = t176 * t24 + t93 * t70;
t195 = -t176 * t25 + t93 * t68;
t193 = -t176 * t52 - t91 * t93;
t192 = -t278 * t90 - t96 * t53;
t141 = t153 * pkin(4);
t124 = t141 + pkin(3);
t156 = -pkin(7) - qJ(4);
t191 = -t95 * t124 - t97 * t156 + t216;
t190 = -qJD(6) * t18 + t272;
t189 = t10 * t159 - t157 * t9;
t188 = -t10 * t157 - t159 * t9;
t37 = -t151 * t56 + t130;
t186 = -t37 * t151 + t38 * t153;
t185 = (-t151 * t72 + t132) * t151 - t58 * t153;
t184 = t83 * t152 - t84 * t154;
t77 = t271 * t151;
t78 = t271 * t153;
t33 = t158 * t77 - t269 * t78;
t86 = t141 + t99;
t36 = pkin(5) * t176 + t96 * pkin(8) + t86;
t16 = t157 * t36 + t159 * t33;
t15 = -t157 * t33 + t159 * t36;
t182 = -qJD(5) * t278 - t96 * qJDD(5);
t181 = t93 * qJD(5) - qJDD(5) * t176;
t180 = -t159 * t50 + (t157 * t90 - t232) * t241;
t179 = t152 * t161 + t223;
t19 = t269 * t45 - t250;
t178 = t158 * t78 + t269 * t77;
t17 = -qJD(5) * pkin(5) - t19;
t177 = pkin(8) * t50 + t241 * t17;
t175 = t17 * qJD(6) * t96 + g(1) * t254;
t174 = -t17 * t278 - t6 * t96 - t273;
t173 = g(1) * t270 + g(2) * t268;
t172 = -t268 * pkin(2) + t199;
t171 = -t201 - t206;
t168 = -g(2) * t254 - qJD(6) * t26 - t267 - t5;
t166 = t97 * t124 - t95 * t156 + t172;
t165 = t194 * t159 + t219 * t241;
t164 = t115 + t171 + t277;
t163 = qJD(6) * t188 + t1 * t159 - t2 * t157;
t162 = qJDD(1) * t99 + t115 - t201 + t59;
t113 = t154 * qJD(2) - qJD(4);
t89 = t90 ^ 2;
t81 = t96 * t152;
t51 = -t91 * pkin(5) - t90 * pkin(8);
t43 = t97 * t157 - t159 * t255;
t42 = t157 * t255 + t97 * t159;
t40 = -t93 * pkin(5) + pkin(8) * t278 + t230;
t22 = qJD(5) * t33 + t113 * t96;
t21 = qJD(5) * t178 + t113 * t176;
t12 = t157 * t51 + t159 * t19;
t11 = -t157 * t19 + t159 * t51;
t4 = -qJD(6) * t16 - t157 * t21 + t159 * t40;
t3 = qJD(6) * t15 + t157 * t40 + t159 * t21;
t23 = [0, 0, 0, 0, 0, qJDD(1), t238, t173, 0, 0, 0, 0, 0, qJDD(1), 0, 0, -qJDD(2) + t238 + 0.2e1 * t245, 0, -t173 + t220 + 0.2e1 * t228, -t203 * pkin(1) - g(1) * t199 - g(2) * t239 + (t220 + t228) * qJ(2), 0, 0, 0, 0, 0, qJDD(1), -t101 * qJDD(1) + 0.2e1 * t115 + t171, t102 * qJDD(1) + 0.2e1 * t117 + t200 + t240, 0, -g(1) * t172 - g(2) * t216 - t184 * qJD(2) + t64 * t101 + t65 * t102, t147 * qJDD(1), 0.2e1 * t151 * t224, 0, t148 * qJDD(1), 0, 0, t162 * t153, -t162 * t151, -t237 * t113 * qJD(1) - t98 * t205 - t186 - t200, t59 * t99 + t71 * t230 - g(1) * (t97 * pkin(3) + t95 * qJ(4) + t172) - g(2) * (-t95 * pkin(3) + t97 * qJ(4) + t216) + t186 * t98 - t185 * t113, t278 * t91 - t52 * t96, t192 + t193, t182, -t176 * t53 + t90 * t93, t181, 0, -t22 * qJD(5) + qJDD(5) * t178 - t134 * t201 + t176 * t48 - t230 * t90 - t86 * t53 - t63 * t93, -t21 * qJD(5) - t33 * qJDD(5) - t230 * t91 - t278 * t63 - t48 * t96 + t86 * t52 + t279, -t176 * t7 - t178 * t52 + t19 * t278 + t20 * t93 + t21 * t90 - t22 * t91 + t33 * t53 + t8 * t96 - t200, -g(1) * t166 - g(2) * t191 + t178 * t8 - t19 * t22 + t20 * t21 + t230 * t63 + t7 * t33 + t48 * t86, t70 * t219 + (t24 * t96 - t278 * t70) * t159 -(-t249 - t251) * t278 + (-t247 + t246 + (t248 - t252) * qJD(6)) * t96, t165 - t196, -t68 * t218 + (-t25 * t96 - t278 * t68) * t157, t195 + t275, -t176 * t50 - t241 * t93, -g(2) * t43 - t15 * t50 + t157 * t174 - t159 * t175 + t176 * t2 - t178 * t25 + t22 * t68 + t241 * t4 - t9 * t93, -g(2) * t42 - t1 * t176 + t10 * t93 + t157 * t175 + t159 * t174 + t16 * t50 + t178 * t24 + t22 * t70 - t241 * t3, t15 * t24 - t16 * t25 - t3 * t68 - t4 * t70 - t188 * t278 - t279 + (qJD(6) * t189 + t1 * t157 + t2 * t159) * t96, t1 * t16 + t10 * t3 + t2 * t15 + t9 * t4 - t6 * t178 + t17 * t22 - g(1) * (t197 * t97 + t166) - g(2) * (-t197 * t95 + t191); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(1), 0, -t161, -t161 * qJ(2) + t203 - t238, 0, 0, 0, 0, 0, 0, -t179, t225 - t243, 0, qJD(1) * t184 + t65 * t152 + t64 * t154 - t238, 0, 0, 0, 0, 0, 0, -t179 * t153, t179 * t151, -t152 * t205 + t237 * t243, -t59 * t154 + t186 * t152 + (-t152 * t71 + t154 * t185) * qJD(1) - t238, 0, 0, 0, 0, 0, 0, t257 * qJD(5) - t81 * qJDD(5) + t154 * t53 + t90 * t236, t256 * qJD(5) - t82 * qJDD(5) - t154 * t52 + t91 * t236, -t256 * t90 + t257 * t91 + t81 * t52 + t82 * t53, -t48 * t154 + t257 * t19 - t256 * t20 - t63 * t236 + t7 * t82 - t8 * t81 - t238, 0, 0, 0, 0, 0, 0, t241 * t258 + t81 * t25 - t257 * t68 - t60 * t50, -t183 * t50 - t81 * t24 - t241 * t259 - t257 * t70, t183 * t25 + t60 * t24 - t258 * t70 - t259 * t68, -t1 * t183 + t259 * t10 - t257 * t17 + t2 * t60 + t258 * t9 + t6 * t81 - t238; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3) + g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t38 * t151 + t37 * t153 + g(3), 0, 0, 0, 0, 0, 0, -t181, t182, -t192 + t193, t176 * t8 - t19 * t93 + t20 * t278 + t7 * t96 + g(3), 0, 0, 0, 0, 0, 0, t195 - t275, t165 + t196 -(t249 - t251) * t278 + (-t247 - t246 + (t248 + t252) * qJD(6)) * t96, t163 * t96 + t17 * t93 - t176 * t6 + t189 * t278 + g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t224, -t226, -t237 * t161, -qJD(1) * t185 + t164, 0, 0, 0, 0, 0, 0, -0.2e1 * t91 * qJD(5) - t187 (t90 + t213) * qJD(5) + t217, -t89 - t274, -t19 * t91 - t20 * t90 + t125 + t164, 0, 0, 0, 0, 0, 0, t180 + t261, t263 - t281 (t24 + t265) * t159 + t70 * t208 + t260, t198 * t157 + t159 * t276 + t17 * t91 - t201; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t262, -t89 + t274 (-t90 + t213) * qJD(5) + t217, -t262, t187, qJDD(5), t63 * t91 + t169 - t211, -t63 * t90 + (t19 + t250) * qJD(5) + t222 - t170, 0, 0, t207 * t70 - t247 (-t24 + t265) * t159 - t241 * t251 + t260, t263 + t281, t208 * t68 - t246, t180 - t261, t241 * t91, -pkin(5) * t25 - t11 * t241 + t177 * t157 + t282 * t159 - t20 * t68 + t9 * t91, pkin(5) * t24 - t10 * t91 + t12 * t241 - t282 * t157 + t177 * t159 - t20 * t70, t11 * t70 + t12 * t68 + ((-t25 + t244) * pkin(8) + t198) * t159 + ((qJD(6) * t68 - t24) * pkin(8) - t276) * t157 + t170, -t10 * t12 - t9 * t11 - t17 * t20 + t167 * pkin(5) + (t163 + t170) * pkin(8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t264, -t68 ^ 2 + t70 ^ 2, t241 * t68 - t24, -t264, t241 * t70 - t25, -t50, -g(1) * t42 + t157 * t168 + t159 * t190 - t17 * t70 + t13 + t283, g(1) * t43 + t17 * t68 + t284 + (-t14 - t190) * t157 + t168 * t159, 0, 0;];
tau_reg  = t23;
