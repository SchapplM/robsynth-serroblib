% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% 
% Output:
% tauc_reg [5x31]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRRRR8_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR8_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR8_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR8_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:25:47
% EndTime: 2019-12-31 22:25:58
% DurationCPUTime: 3.66s
% Computational Cost: add. (4351->330), mult. (10590->473), div. (0->0), fcn. (7896->8), ass. (0->195)
t168 = qJD(2) + qJD(3);
t173 = sin(qJ(3));
t177 = cos(qJ(2));
t265 = cos(qJ(3));
t215 = qJD(1) * t265;
t174 = sin(qJ(2));
t233 = qJD(1) * t174;
t284 = -t173 * t233 + t177 * t215;
t90 = t284 * t168;
t268 = qJD(4) + qJD(5);
t285 = t284 - t268;
t214 = t265 * qJD(3);
t267 = pkin(6) + pkin(7);
t153 = t267 * t177;
t143 = qJD(1) * t153;
t129 = t173 * t143;
t151 = t267 * t174;
t141 = qJD(1) * t151;
t98 = -t265 * t141 - t129;
t280 = pkin(2) * t214 - t98;
t172 = sin(qJ(4));
t176 = cos(qJ(4));
t239 = t173 * t177;
t128 = -qJD(1) * t239 - t174 * t215;
t92 = -t128 * pkin(3) - pkin(8) * t284;
t78 = pkin(2) * t233 + t92;
t283 = -t280 * t172 - t176 * t78;
t107 = -t172 * t128 - t176 * t168;
t171 = sin(qJ(5));
t175 = cos(qJ(5));
t193 = t176 * t128 - t172 * t168;
t194 = t171 * t107 + t175 * t193;
t56 = t175 * t107 - t171 * t193;
t282 = t194 * t56;
t136 = t171 * t172 - t175 * t176;
t253 = t285 * t136;
t138 = t171 * t176 + t175 * t172;
t281 = t285 * t138;
t279 = t194 ^ 2 - t56 ^ 2;
t229 = qJD(5) * t171;
t166 = -t177 * pkin(2) - pkin(1);
t149 = t166 * qJD(1);
t76 = -pkin(3) * t284 + t128 * pkin(8) + t149;
t130 = t265 * t143;
t260 = qJD(2) * pkin(2);
t131 = -t141 + t260;
t96 = t173 * t131 + t130;
t80 = t168 * pkin(8) + t96;
t36 = t172 * t76 + t176 * t80;
t26 = -t107 * pkin(9) + t36;
t24 = t26 * t229;
t95 = t265 * t131 - t129;
t79 = -t168 * pkin(3) - t95;
t53 = t107 * pkin(4) + t79;
t278 = t53 * t56 + t24;
t123 = qJD(4) - t284;
t118 = qJD(5) + t123;
t228 = qJD(5) * t175;
t230 = qJD(4) * t176;
t231 = qJD(4) * t172;
t51 = t128 * t231 + t168 * t230 + t176 * t90;
t52 = -t193 * qJD(4) + t172 * t90;
t13 = -t107 * t228 - t171 * t52 + t175 * t51 + t193 * t229;
t277 = t56 * t118 + t13;
t227 = qJD(1) * qJD(2);
t213 = t174 * t227;
t139 = t265 * t174 + t239;
t103 = t168 * t139;
t91 = t103 * qJD(1);
t40 = pkin(2) * t213 + t91 * pkin(3) - t90 * pkin(8);
t38 = t176 * t40;
t221 = qJD(2) * t267;
t204 = qJD(1) * t221;
t132 = t174 * t204;
t133 = t177 * t204;
t232 = qJD(3) * t173;
t43 = t131 * t214 - t265 * t132 - t173 * t133 - t143 * t232;
t182 = -t36 * qJD(4) - t172 * t43 + t38;
t4 = t91 * pkin(4) - t51 * pkin(9) + t182;
t191 = t172 * t40 + t176 * t43 + t76 * t230 - t80 * t231;
t5 = -t52 * pkin(9) + t191;
t222 = -t171 * t5 + t175 * t4;
t35 = -t172 * t80 + t176 * t76;
t25 = pkin(9) * t193 + t35;
t23 = t123 * pkin(4) + t25;
t257 = t175 * t26;
t9 = t171 * t23 + t257;
t276 = -t9 * qJD(5) + t53 * t194 + t222;
t181 = t194 * qJD(5) - t171 * t51 - t175 * t52;
t275 = -t118 * t194 + t181;
t274 = -0.2e1 * t227;
t44 = t131 * t232 - t173 * t132 + t265 * t133 + t143 * t214;
t273 = t44 * t176 - t79 * t231;
t97 = -t173 * t141 + t130;
t203 = pkin(2) * t232 - t97;
t85 = t138 * t139;
t244 = t284 * t172;
t272 = (t231 - t244) * pkin(4);
t190 = -t173 * t174 + t265 * t177;
t102 = t168 * t190;
t238 = t176 * t102;
t188 = -t139 * t231 + t238;
t271 = -t265 * t151 - t173 * t153;
t270 = t172 * t78 - t280 * t176;
t269 = qJD(1) * t139;
t266 = -pkin(8) - pkin(9);
t264 = t176 * pkin(4);
t163 = t173 * pkin(2) + pkin(8);
t263 = -pkin(9) - t163;
t262 = t172 * t92 + t176 * t95;
t258 = t172 * t91;
t256 = t176 * t91;
t254 = t51 * t172;
t111 = -t173 * t151 + t265 * t153;
t106 = t176 * t111;
t94 = -pkin(3) * t190 - t139 * pkin(8) + t166;
t251 = t172 * t94 + t106;
t250 = t272 + t203;
t249 = t102 * t172;
t248 = t107 * t123;
t247 = t193 * t123;
t246 = t118 * t128;
t245 = t123 * t128;
t243 = t284 * t176;
t242 = t128 * t284;
t241 = t139 * t172;
t240 = t139 * t176;
t179 = qJD(1) ^ 2;
t237 = t177 * t179;
t178 = qJD(2) ^ 2;
t236 = t178 * t174;
t235 = t178 * t177;
t234 = t174 ^ 2 - t177 ^ 2;
t226 = pkin(9) * t244;
t225 = t174 * t260;
t220 = qJD(4) * t266;
t218 = t139 * t230;
t212 = qJD(5) * t23 + t5;
t210 = -t172 * t95 + t176 * t92;
t209 = qJD(4) * t263;
t208 = pkin(1) * t274;
t207 = t123 * t176;
t164 = -t265 * pkin(2) - pkin(3);
t205 = -t36 * t128 + t44 * t172 + t79 * t230;
t202 = -t96 + t272;
t201 = -t128 * pkin(4) - pkin(9) * t243;
t134 = t263 * t172;
t200 = -qJD(5) * t134 - t172 * t209 - t226 + t270;
t167 = t176 * pkin(9);
t135 = t176 * t163 + t167;
t199 = qJD(5) * t135 - t176 * t209 + t201 - t283;
t150 = t266 * t172;
t198 = -qJD(5) * t150 - t172 * t220 - t226 + t262;
t152 = t176 * pkin(8) + t167;
t197 = qJD(5) * t152 - t176 * t220 + t201 + t210;
t196 = -t163 * t91 - t284 * t79;
t192 = t35 * t128 - t273;
t189 = t218 + t249;
t50 = t103 * pkin(3) - t102 * pkin(8) + t225;
t142 = t174 * t221;
t144 = t177 * t221;
t60 = t271 * qJD(3) - t265 * t142 - t173 * t144;
t187 = -t111 * t231 + t172 * t50 + t176 * t60 + t94 * t230;
t20 = t52 * pkin(4) + t44;
t8 = -t171 * t26 + t175 * t23;
t186 = t8 * t128 + t20 * t136 - t281 * t53;
t185 = -t9 * t128 + t20 * t138 + t253 * t53;
t184 = t149 * t128 - t44;
t180 = -t149 * t284 - t43;
t61 = t111 * qJD(3) - t173 * t142 + t265 * t144;
t165 = -pkin(3) - t264;
t148 = t164 - t264;
t88 = t176 * t94;
t86 = t136 * t139;
t77 = pkin(4) * t241 - t271;
t67 = t128 ^ 2 - t284 ^ 2;
t66 = t91 * t190;
t63 = (-t128 - t269) * t168;
t46 = t176 * t50;
t39 = -pkin(9) * t241 + t251;
t31 = -pkin(4) * t190 - pkin(9) * t240 - t172 * t111 + t88;
t29 = t189 * pkin(4) + t61;
t22 = t123 * t207 - t128 * t193 + t258;
t21 = -t123 ^ 2 * t172 - t107 * t128 + t256;
t19 = -t193 * t207 + t254;
t18 = -t229 * t241 + (t268 * t240 + t249) * t175 + t188 * t171;
t17 = -t136 * t102 - t268 * t85;
t12 = t118 * t281 - t56 * t128 - t136 * t91;
t11 = t253 * t118 - t128 * t194 + t138 * t91;
t10 = -t189 * pkin(9) + t187;
t7 = -pkin(9) * t238 + t103 * pkin(4) - t172 * t60 + t46 + (-t106 + (pkin(9) * t139 - t94) * t172) * qJD(4);
t6 = (t51 - t248) * t176 + (-t52 + t247) * t172;
t2 = t13 * t138 - t194 * t253;
t1 = -t13 * t136 + t138 * t181 - t194 * t281 - t253 * t56;
t3 = [0, 0, 0, 0.2e1 * t177 * t213, t234 * t274, t235, -t236, 0, -pkin(6) * t235 + t174 * t208, pkin(6) * t236 + t177 * t208, -t128 * t102 + t90 * t139, t102 * t284 + t128 * t103 - t139 * t91 + t190 * t90, t102 * t168, -t103 * t168, 0, t149 * t103 + t166 * t91 - t61 * t168 + (-qJD(1) * t190 - t284) * t225, t149 * t102 + t166 * t90 - t60 * t168 + (-t128 + t269) * t225, -t188 * t193 + t51 * t240, (-t107 * t176 + t172 * t193) * t102 + (-t254 - t176 * t52 + (t107 * t172 + t176 * t193) * qJD(4)) * t139, -t103 * t193 + t188 * t123 - t190 * t51 + t91 * t240, -t107 * t103 - t189 * t123 + t190 * t52 - t91 * t241, t123 * t103 - t66, t61 * t107 - t271 * t52 + t79 * t218 + (-t111 * t230 + t46) * t123 + t88 * t91 - (-t80 * t230 + t38) * t190 + t35 * t103 + (t44 * t139 + t79 * t102 + (-qJD(4) * t94 - t60) * t123 - t111 * t91 - (-qJD(4) * t76 - t43) * t190) * t172, -t36 * t103 - t187 * t123 + t273 * t139 + t190 * t191 - t193 * t61 + t79 * t238 - t251 * t91 - t271 * t51, -t13 * t86 - t17 * t194, -t13 * t85 - t17 * t56 + t18 * t194 - t181 * t86, -t103 * t194 + t17 * t118 - t13 * t190 - t86 * t91, -t56 * t103 - t18 * t118 - t181 * t190 - t85 * t91, t118 * t103 - t66, (-t171 * t10 + t175 * t7) * t118 + (-t171 * t39 + t175 * t31) * t91 - t222 * t190 + t8 * t103 + t29 * t56 - t77 * t181 + t20 * t85 + t53 * t18 + ((-t171 * t31 - t175 * t39) * t118 + t9 * t190) * qJD(5), -t9 * t103 + t77 * t13 - t24 * t190 + t53 * t17 - t20 * t86 - t29 * t194 + (-(-qJD(5) * t39 + t7) * t118 - t31 * t91 + t4 * t190) * t171 + (-(qJD(5) * t31 + t10) * t118 - t39 * t91 + t212 * t190) * t175; 0, 0, 0, -t174 * t237, t234 * t179, 0, 0, 0, t179 * pkin(1) * t174, pkin(1) * t237, t242, t67, 0, t63, 0, t97 * t168 + (-t168 * t232 + t233 * t284) * pkin(2) + t184, t98 * t168 + (t128 * t233 - t168 * t214) * pkin(2) + t180, t19, t6, t22, t21, t245, t164 * t52 + t196 * t172 + t203 * t107 + (-t163 * t230 + t283) * t123 + t192, t164 * t51 + t196 * t176 - t203 * t193 + (t163 * t231 + t270) * t123 + t205, t2, t1, t11, t12, t246, (t175 * t134 - t171 * t135) * t91 - t148 * t181 + t250 * t56 + (t200 * t171 - t199 * t175) * t118 + t186, -(t171 * t134 + t175 * t135) * t91 + t148 * t13 - t250 * t194 + (t199 * t171 + t200 * t175) * t118 + t185; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t242, t67, 0, t63, 0, t96 * t168 + t184, t95 * t168 + t180, t19, t6, t22, t21, t245, -pkin(3) * t52 - t96 * t107 - t79 * t244 - t210 * t123 + (-t123 * t230 - t258) * pkin(8) + t192, -pkin(3) * t51 + t96 * t193 - t79 * t243 + t262 * t123 + (t123 * t231 - t256) * pkin(8) + t205, t2, t1, t11, t12, t246, (t175 * t150 - t171 * t152) * t91 - t165 * t181 + t202 * t56 + (t198 * t171 - t197 * t175) * t118 + t186, -(t171 * t150 + t175 * t152) * t91 + t165 * t13 - t202 * t194 + (t171 * t197 + t175 * t198) * t118 + t185; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t193 * t107, -t107 ^ 2 + t193 ^ 2, t51 + t248, -t52 - t247, t91, t36 * t123 + t193 * t79 + t182, t79 * t107 + t35 * t123 - t191, -t282, t279, t277, t275, t91, -(-t171 * t25 - t257) * t118 + (-t118 * t229 + t175 * t91 + t193 * t56) * pkin(4) + t276, (-t26 * t118 - t4) * t171 + (t25 * t118 - t212) * t175 + (-t118 * t228 - t171 * t91 - t193 * t194) * pkin(4) + t278; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t282, t279, t277, t275, t91, t9 * t118 + t276, t8 * t118 - t171 * t4 - t175 * t212 + t278;];
tauc_reg = t3;
