% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRPRP9
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% tau_reg [5x25]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPRP9_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP9_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP9_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP9_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP9_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP9_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:07:31
% EndTime: 2019-12-31 20:07:40
% DurationCPUTime: 3.49s
% Computational Cost: add. (3707->412), mult. (8333->526), div. (0->0), fcn. (5789->10), ass. (0->194)
t163 = cos(qJ(2));
t148 = t163 * qJDD(1);
t161 = sin(qJ(2));
t222 = qJD(1) * qJD(2);
t183 = t161 * t222 - t148;
t108 = qJDD(4) + t183;
t154 = pkin(8) + qJ(4);
t146 = sin(t154);
t153 = g(3) * t163;
t162 = sin(qJ(1));
t164 = cos(qJ(1));
t203 = g(1) * t164 + g(2) * t162;
t192 = t203 * t161;
t175 = t192 - t153;
t157 = sin(pkin(8));
t254 = pkin(7) + qJ(3);
t120 = t254 * t157;
t158 = cos(pkin(8));
t121 = t254 * t158;
t160 = sin(qJ(4));
t264 = cos(qJ(4));
t69 = -t160 * t120 + t264 * t121;
t276 = -t69 * t108 - t175 * t146;
t229 = qJD(1) * t161;
t213 = t157 * t229;
t226 = qJD(2) * t158;
t105 = -t213 + t226;
t216 = t158 * t229;
t227 = qJD(2) * t157;
t106 = t216 + t227;
t55 = -t264 * t105 + t106 * t160;
t275 = t55 ^ 2;
t188 = -t160 * t105 - t264 * t106;
t265 = t188 ^ 2;
t217 = t264 * t158;
t186 = -t160 * t157 + t217;
t228 = qJD(1) * t163;
t212 = qJD(4) * t264;
t223 = qJD(4) * t160;
t269 = -t157 * t223 + t158 * t212;
t249 = -t186 * t228 + t269;
t110 = t264 * t157 + t160 * t158;
t178 = t163 * t110;
t97 = t110 * qJD(4);
t248 = -qJD(1) * t178 + t97;
t137 = -qJD(4) + t228;
t274 = t137 * t55;
t273 = t137 * t188;
t210 = t163 * t222;
t220 = t161 * qJDD(1);
t184 = t210 + t220;
t187 = -t264 * t120 - t160 * t121;
t240 = t158 * t163;
t198 = pkin(3) * t161 - pkin(7) * t240;
t200 = pkin(2) * t161 - qJ(3) * t163;
t112 = t200 * qJD(1);
t74 = pkin(6) * t213 + t158 * t112;
t47 = t198 * qJD(1) + t74;
t241 = t158 * t161;
t242 = t157 * t163;
t185 = -pkin(6) * t241 - pkin(7) * t242;
t98 = t157 * t112;
t60 = t185 * qJD(1) + t98;
t272 = t186 * qJD(3) + t187 * qJD(4) - t160 * t47 - t264 * t60;
t271 = t110 * qJD(3) + t69 * qJD(4) - t160 * t60 + t264 * t47;
t268 = pkin(6) * t210 + qJDD(3);
t221 = qJDD(2) * t157;
t171 = t184 * t158 + t221;
t233 = t184 * t157;
t199 = qJDD(2) * t158 - t233;
t19 = -t188 * qJD(4) + t160 * t171 - t264 * t199;
t201 = pkin(2) * t163 + qJ(3) * t161;
t117 = -pkin(1) - t201;
t104 = t158 * t117;
t61 = -pkin(7) * t241 + t104 + (-pkin(6) * t157 - pkin(3)) * t163;
t243 = t157 * t161;
t78 = pkin(6) * t240 + t157 * t117;
t67 = -pkin(7) * t243 + t78;
t193 = t160 * t61 + t264 * t67;
t225 = qJD(2) * t161;
t218 = pkin(6) * t225;
t94 = t200 * qJD(2) - qJD(3) * t161;
t65 = t157 * t218 + t158 * t94;
t39 = t198 * qJD(2) + t65;
t81 = t157 * t94;
t48 = t185 * qJD(2) + t81;
t267 = -t193 * qJD(4) - t160 * t48 + t264 * t39;
t266 = -0.2e1 * pkin(1);
t263 = pkin(3) * t157;
t262 = pkin(4) * t108;
t261 = pkin(6) * t105;
t260 = g(1) * t162;
t257 = g(2) * t164;
t256 = g(3) * t161;
t255 = t188 * t55;
t144 = pkin(6) * t228;
t219 = pkin(3) * t228;
t101 = t157 * t219 + t144;
t253 = t248 * pkin(4) - t249 * qJ(5) - qJD(5) * t110 - t101;
t252 = -qJ(5) * t229 + t272;
t251 = pkin(4) * t229 + t271;
t52 = t94 * qJD(1) + t117 * qJDD(1);
t83 = -t183 * pkin(6) + qJDD(2) * qJ(3) + qJD(2) * qJD(3);
t27 = t157 * t52 + t158 * t83;
t100 = t117 * qJD(1);
t122 = qJD(2) * qJ(3) + t144;
t62 = t158 * t100 - t122 * t157;
t31 = -pkin(7) * t106 - t219 + t62;
t63 = t157 * t100 + t158 * t122;
t35 = pkin(7) * t105 + t63;
t11 = t160 * t31 + t264 * t35;
t247 = t11 * t137;
t245 = qJ(5) * t108;
t244 = qJDD(2) * pkin(2);
t239 = t161 * t164;
t147 = cos(t154);
t238 = t162 * t147;
t237 = t162 * t163;
t236 = t163 * t164;
t235 = t164 * t146;
t10 = -t160 * t35 + t264 * t31;
t234 = qJD(5) - t10;
t142 = pkin(6) * t220;
t232 = -t142 - t153;
t224 = qJD(2) * t163;
t215 = t157 * t224;
t102 = pkin(3) * t215 + pkin(6) * t224;
t113 = pkin(3) * t243 + t161 * pkin(6);
t231 = t164 * pkin(1) + t162 * pkin(6);
t155 = t161 ^ 2;
t230 = -t163 ^ 2 + t155;
t141 = t158 * pkin(3) + pkin(2);
t26 = -t157 * t83 + t158 * t52;
t14 = t183 * pkin(3) - t171 * pkin(7) + t26;
t20 = t199 * pkin(7) + t27;
t209 = -t264 * t14 + t160 * t20 + t35 * t212 + t31 * t223;
t208 = -qJD(2) * pkin(2) + qJD(3);
t84 = t146 * t237 + t147 * t164;
t86 = t163 * t235 - t238;
t206 = -g(1) * t84 + g(2) * t86;
t85 = t147 * t237 - t235;
t87 = t146 * t162 + t147 * t236;
t205 = g(1) * t85 - g(2) * t87;
t139 = t161 * t260;
t204 = -g(2) * t239 + t139;
t202 = -t257 + t260;
t116 = pkin(6) * t229 + t208;
t197 = pkin(4) * t147 + qJ(5) * t146 + t141;
t92 = t142 - t244 + t268;
t195 = -t160 * t67 + t264 * t61;
t191 = -pkin(6) * qJDD(2) + t222 * t266;
t190 = t160 * t14 + t264 * t20 + t31 * t212 - t35 * t223;
t189 = t160 * t39 + t61 * t212 - t67 * t223 + t264 * t48;
t18 = -t105 * t212 + t106 * t223 - t160 * t199 - t264 * t171;
t182 = t158 * t220 + t221;
t166 = qJD(1) ^ 2;
t181 = pkin(1) * t166 + t203;
t73 = -pkin(3) * t105 + t116;
t165 = qJD(2) ^ 2;
t180 = pkin(6) * t165 + qJDD(1) * t266 + t257;
t179 = g(2) * t161 * t238 + t187 * t108 + (g(1) * t239 - t153) * t147;
t177 = t18 - t274;
t176 = g(1) * t86 + g(2) * t84 + t146 * t256 - t209;
t174 = -t203 * t163 - t256;
t173 = -t192 - t244;
t172 = -t92 + t175;
t15 = pkin(4) * t55 + qJ(5) * t188 + t73;
t170 = -t15 * t188 + qJDD(5) - t176;
t40 = -t199 * pkin(3) + t92;
t169 = -g(1) * t87 - g(2) * t85 - t147 * t256 + t190;
t3 = t19 * pkin(4) + t18 * qJ(5) + qJD(5) * t188 + t40;
t167 = t19 + t273;
t151 = t164 * pkin(6);
t91 = t186 * t161;
t90 = t110 * t161;
t77 = -pkin(6) * t242 + t104;
t75 = -pkin(6) * t216 + t98;
t66 = -t158 * t218 + t81;
t53 = -pkin(4) * t186 - qJ(5) * t110 - t141;
t43 = qJD(2) * t178 + t269 * t161;
t42 = t160 * t215 + t161 * t97 - t217 * t224;
t32 = pkin(4) * t90 - qJ(5) * t91 + t113;
t25 = -pkin(4) * t188 + qJ(5) * t55;
t24 = t163 * pkin(4) - t195;
t23 = -qJ(5) * t163 + t193;
t9 = -t137 * qJ(5) + t11;
t8 = t137 * pkin(4) + t234;
t7 = pkin(4) * t43 + qJ(5) * t42 - qJD(5) * t91 + t102;
t6 = -t18 - t274;
t5 = -pkin(4) * t225 - t267;
t4 = qJ(5) * t225 - qJD(5) * t163 + t189;
t2 = qJDD(5) + t209 - t262;
t1 = -qJD(5) * t137 + t190 + t245;
t12 = [qJDD(1), t202, t203, qJDD(1) * t155 + 0.2e1 * t161 * t210, 0.2e1 * t161 * t148 - 0.2e1 * t230 * t222, qJDD(2) * t161 + t163 * t165, qJDD(2) * t163 - t161 * t165, 0, t191 * t161 + (-t180 + t260) * t163, t180 * t161 + t191 * t163 - t139, -t203 * t157 + (-pkin(6) * t199 + t92 * t157 + (qJD(1) * t77 + t62) * qJD(2)) * t161 + (-t65 * qJD(1) - t77 * qJDD(1) - t26 + t202 * t158 + (t116 * t157 - t261) * qJD(2)) * t163, -t203 * t158 + (t92 * t158 + (-qJD(1) * t78 - t63) * qJD(2) + t182 * pkin(6)) * t161 + (t66 * qJD(1) + t78 * qJDD(1) + t27 - t202 * t157 + (t116 * t158 + (t106 + t216) * pkin(6)) * qJD(2)) * t163, t66 * t105 - t78 * t233 - t65 * t106 + (-qJDD(2) * t77 - t161 * t27 - t63 * t224) * t157 + (t78 * qJDD(2) - t26 * t161 - t77 * t184 - t62 * t224) * t158 + t204, t27 * t78 + t63 * t66 + t26 * t77 + t62 * t65 - g(1) * t151 - g(2) * (t201 * t164 + t231) - t117 * t260 + (t116 * t224 + t92 * t161) * pkin(6), -t18 * t91 + t188 * t42, t18 * t90 + t188 * t43 - t19 * t91 + t42 * t55, t108 * t91 + t137 * t42 + t163 * t18 - t188 * t225, -t108 * t90 + t137 * t43 + t163 * t19 - t55 * t225, -t108 * t163 - t137 * t225, t10 * t225 + t102 * t55 + t195 * t108 + t113 * t19 - t267 * t137 + t209 * t163 + t40 * t90 + t73 * t43 + t205, -t102 * t188 - t193 * t108 - t11 * t225 - t113 * t18 + t189 * t137 + t190 * t163 + t40 * t91 - t73 * t42 + t206, -t108 * t24 + t137 * t5 + t15 * t43 + t163 * t2 + t19 * t32 - t8 * t225 + t3 * t90 + t55 * t7 + t205, -t1 * t90 - t18 * t24 - t188 * t5 - t19 * t23 + t2 * t91 - t4 * t55 - t42 * t8 - t43 * t9 + t204, -t1 * t163 + t108 * t23 - t137 * t4 + t15 * t42 + t18 * t32 + t188 * t7 + t225 * t9 - t3 * t91 - t206, t1 * t23 + t9 * t4 + t3 * t32 + t15 * t7 + t2 * t24 + t8 * t5 - g(1) * (-pkin(4) * t85 - qJ(5) * t84 + t164 * t263 + t151) - g(2) * (pkin(4) * t87 + qJ(5) * t86 + t141 * t236 + t239 * t254 + t231) + (-g(1) * (-t141 * t163 - t161 * t254 - pkin(1)) - g(2) * t263) * t162; 0, 0, 0, -t161 * t166 * t163, t230 * t166, t220, t148, qJDD(2), t181 * t161 + t232, t256 + (-pkin(6) * qJDD(1) + t181) * t163, t157 * qJ(3) * t148 - pkin(2) * t233 + (t172 + t244) * t158 + ((-qJ(3) * t227 - t62) * t161 + (t261 + t74 + (qJD(3) - t116) * t157) * t163) * qJD(1), -t200 * t158 * qJDD(1) + (t173 + t92 + t153) * t157 + ((-qJ(3) * t226 + t63) * t161 + (-pkin(6) * t106 - t75 + (-t116 + t208) * t158) * t163) * qJD(1), -t75 * t105 + t74 * t106 + (t199 * qJ(3) + qJD(3) * t105 + t62 * t228 + t27) * t158 + (t171 * qJ(3) + qJD(3) * t106 + t63 * t228 - t26) * t157 + t174, -t116 * t144 - t62 * t74 - t63 * t75 + (-t157 * t62 + t158 * t63) * qJD(3) + t172 * pkin(2) + (-t26 * t157 + t27 * t158 + t174) * qJ(3), -t110 * t18 - t188 * t249, -t110 * t19 - t18 * t186 + t188 * t248 - t249 * t55, t108 * t110 - t249 * t137 + t188 * t229, t108 * t186 + t248 * t137 + t55 * t229, t137 * t229, -t10 * t229 - t101 * t55 + t271 * t137 - t141 * t19 - t186 * t40 + t248 * t73 + t179, t101 * t188 + t11 * t229 + t40 * t110 + t272 * t137 + t141 * t18 + t249 * t73 + t276, t251 * t137 + t248 * t15 - t186 * t3 + t19 * t53 + t8 * t229 + t253 * t55 + t179, t1 * t186 + t110 * t2 + t18 * t187 - t188 * t251 - t19 * t69 - t248 * t9 + t249 * t8 - t252 * t55 + t174, -t110 * t3 - t252 * t137 - t249 * t15 + t18 * t53 + t188 * t253 - t9 * t229 - t276, t1 * t69 - t2 * t187 + t3 * t53 + t252 * t9 + t251 * t8 + t253 * t15 + (-g(3) * t197 - t203 * t254) * t163 + (-g(3) * t254 + t203 * t197) * t161; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t106 * t228 - t199, (-t105 + t226) * t228 + t182, -t105 ^ 2 - t106 ^ 2, -t105 * t63 + t106 * t62 + t173 - t232 + t268, 0, 0, 0, 0, 0, t167, -t177, t167, -t265 - t275, t177, t188 * t8 + t55 * t9 - t175 + t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t255, t265 - t275, t6, -t19 + t273, t108, t188 * t73 + t176 - t247, -t10 * t137 + t55 * t73 - t169, -t25 * t55 - t170 - t247 + 0.2e1 * t262, pkin(4) * t18 - qJ(5) * t19 - (-t11 + t9) * t188 + (t8 - t234) * t55, 0.2e1 * t245 - t15 * t55 - t25 * t188 + (-0.2e1 * qJD(5) + t10) * t137 + t169, t1 * qJ(5) - t2 * pkin(4) - t15 * t25 - t8 * t11 - g(1) * (-pkin(4) * t86 + qJ(5) * t87) - g(2) * (-pkin(4) * t84 + qJ(5) * t85) + t234 * t9 - (-pkin(4) * t146 + qJ(5) * t147) * t256; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t108 - t255, t6, -t137 ^ 2 - t265, t137 * t9 + t170 - t262;];
tau_reg = t12;
