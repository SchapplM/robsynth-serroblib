% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5,theta4]';
% 
% Output:
% tauc_reg [5x26]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRRPR10_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR10_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR10_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR10_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:29:51
% EndTime: 2019-12-31 21:30:03
% DurationCPUTime: 4.03s
% Computational Cost: add. (5202->358), mult. (14212->525), div. (0->0), fcn. (11121->10), ass. (0->190)
t174 = sin(pkin(5));
t182 = cos(qJ(2));
t243 = qJD(1) * t182;
t226 = t174 * t243;
t283 = qJD(3) - t226;
t178 = sin(qJ(3));
t181 = cos(qJ(3));
t208 = t178 * t226;
t267 = -qJ(4) - pkin(8);
t218 = qJD(3) * t267;
t179 = sin(qJ(2));
t245 = qJD(1) * t174;
t227 = t179 * t245;
t176 = cos(pkin(5));
t244 = qJD(1) * t176;
t233 = pkin(1) * t244;
t130 = -pkin(7) * t227 + t182 * t233;
t195 = (pkin(2) * t179 - pkin(8) * t182) * t174;
t131 = qJD(1) * t195;
t247 = t181 * t130 + t178 * t131;
t282 = qJ(4) * t208 + t181 * qJD(4) + t178 * t218 - t247;
t212 = -t178 * t130 + t181 * t131;
t250 = t181 * t182;
t281 = -t178 * qJD(4) + t181 * t218 - (pkin(3) * t179 - qJ(4) * t250) * t245 - t212;
t164 = qJD(2) + t244;
t209 = t178 * t227;
t113 = -t181 * t164 + t209;
t115 = t164 * t178 + t181 * t227;
t173 = sin(pkin(10));
t175 = cos(pkin(10));
t214 = -t175 * t113 - t115 * t173;
t274 = qJD(5) - t214;
t177 = sin(qJ(5));
t180 = cos(qJ(5));
t197 = -t113 * t173 + t175 * t115;
t59 = t177 * t197 - t180 * t283;
t280 = t274 * t59;
t143 = t173 * t181 + t175 * t178;
t259 = t283 * t143;
t142 = t173 * t178 - t175 * t181;
t138 = t142 * qJD(3);
t97 = t142 * t226;
t279 = t138 - t97;
t198 = -t177 * t283 - t180 * t197;
t278 = t198 * t274;
t216 = t180 * t274;
t234 = qJD(1) * qJD(2);
t219 = t174 * t234;
t206 = t182 * t219;
t238 = qJD(3) * t181;
t89 = -qJD(3) * t209 + t164 * t238 + t181 * t206;
t240 = qJD(2) * t182;
t223 = t178 * t240;
t239 = qJD(3) * t178;
t90 = (t179 * t238 + t223) * t245 + t164 * t239;
t55 = t173 * t89 + t175 * t90;
t261 = t177 * t55;
t277 = -t216 * t274 - t261;
t266 = -t282 * t173 + t281 * t175;
t264 = t281 * t173 + t282 * t175;
t251 = t174 * t182;
t127 = pkin(7) * t251 + (pkin(1) * t179 + pkin(8)) * t176;
t128 = (-pkin(2) * t182 - pkin(8) * t179 - pkin(1)) * t174;
t248 = t181 * t127 + t178 * t128;
t161 = t179 * t233;
t133 = pkin(7) * t226 + t161;
t275 = -t133 + (-t208 + t239) * pkin(3);
t168 = pkin(3) * t173 + pkin(9);
t205 = t179 * t219;
t132 = qJD(2) * t195;
t122 = qJD(1) * t132;
t252 = t174 * t179;
t165 = pkin(7) * t252;
t271 = pkin(1) * t182;
t134 = (t176 * t271 - t165) * qJD(2);
t123 = qJD(1) * t134;
t102 = pkin(8) * t164 + t133;
t109 = qJD(1) * t128;
t71 = t102 * t181 + t109 * t178;
t185 = -qJD(3) * t71 + t181 * t122 - t178 * t123;
t21 = pkin(3) * t205 - qJ(4) * t89 - qJD(4) * t115 + t185;
t194 = -t102 * t239 + t109 * t238 + t178 * t122 + t181 * t123;
t25 = -qJ(4) * t90 - qJD(4) * t113 + t194;
t5 = -t173 * t25 + t175 * t21;
t3 = -pkin(4) * t205 - t5;
t273 = (pkin(3) * t115 + pkin(4) * t197 - pkin(9) * t214 + qJD(5) * t168) * t274 + t3;
t158 = t267 * t181;
t221 = t267 * t178;
t100 = -t175 * t158 + t173 * t221;
t228 = -pkin(3) * t181 - pkin(2);
t91 = pkin(4) * t142 - pkin(9) * t143 + t228;
t272 = (-t259 * pkin(4) - t279 * pkin(9) + qJD(5) * t100 - t275) * t274 - t91 * t55;
t56 = -t173 * t90 + t175 * t89;
t27 = -qJD(5) * t198 + t177 * t56 - t180 * t205;
t124 = pkin(7) * t206 + qJD(2) * t161;
t68 = pkin(3) * t90 + t124;
t14 = pkin(4) * t55 - pkin(9) * t56 + t68;
t70 = -t102 * t178 + t181 * t109;
t53 = -qJ(4) * t115 + t70;
t47 = pkin(3) * t283 + t53;
t54 = -qJ(4) * t113 + t71;
t49 = t175 * t54;
t20 = t173 * t47 + t49;
t16 = pkin(9) * t283 + t20;
t101 = -pkin(2) * t164 - t130;
t76 = pkin(3) * t113 + qJD(4) + t101;
t28 = -pkin(4) * t214 - pkin(9) * t197 + t76;
t201 = t16 * t177 - t180 * t28;
t6 = t173 * t21 + t175 * t25;
t4 = pkin(9) * t205 + t6;
t1 = -qJD(5) * t201 + t177 * t14 + t180 * t4;
t270 = t59 * t197;
t269 = t198 * t197;
t140 = t176 * t178 + t181 * t252;
t184 = -t248 * qJD(3) + t181 * t132 - t178 * t134;
t242 = qJD(2) * t179;
t225 = t174 * t242;
t139 = -t176 * t181 + t178 * t252;
t224 = t174 * t240;
t94 = -qJD(3) * t139 + t181 * t224;
t32 = pkin(3) * t225 - qJ(4) * t94 - qJD(4) * t140 + t184;
t193 = -t127 * t239 + t128 * t238 + t178 * t132 + t181 * t134;
t93 = qJD(3) * t140 + t174 * t223;
t36 = -qJ(4) * t93 - qJD(4) * t139 + t193;
t12 = t173 * t32 + t175 * t36;
t213 = -t127 * t178 + t181 * t128;
t58 = -pkin(3) * t251 - qJ(4) * t140 + t213;
t66 = -qJ(4) * t139 + t248;
t35 = t173 * t58 + t175 * t66;
t265 = pkin(4) * t227 - t266;
t263 = t100 * t55;
t262 = t173 * t54;
t236 = qJD(5) * t180;
t237 = qJD(5) * t177;
t26 = t177 * t205 + t180 * t56 - t197 * t237 + t236 * t283;
t260 = t26 * t177;
t258 = t113 * t283;
t257 = t115 * t283;
t256 = t143 * t180;
t255 = t283 * t178;
t254 = t283 * t181;
t170 = t174 ^ 2;
t253 = t170 * qJD(1) ^ 2;
t135 = t176 * pkin(1) * t242 + pkin(7) * t224;
t246 = t179 ^ 2 - t182 ^ 2;
t241 = qJD(2) * t181;
t235 = qJD(2) - t164;
t231 = t179 * t253;
t230 = t177 * t251;
t220 = t170 * t234;
t211 = t164 + t244;
t210 = 0.2e1 * t220;
t207 = pkin(3) * t93 + t135;
t203 = -0.2e1 * pkin(1) * t220;
t8 = t16 * t180 + t177 * t28;
t11 = -t173 * t36 + t175 * t32;
t19 = t175 * t47 - t262;
t34 = -t173 * t66 + t175 * t58;
t31 = -pkin(9) * t251 + t35;
t126 = t165 + (-pkin(2) - t271) * t176;
t188 = t139 * pkin(3) + t126;
t86 = t175 * t139 + t140 * t173;
t87 = -t139 * t173 + t140 * t175;
t44 = pkin(4) * t86 - pkin(9) * t87 + t188;
t200 = t177 * t44 + t180 * t31;
t199 = -t177 * t31 + t180 * t44;
t196 = t180 * t55 + (t177 * t214 - t237) * t274;
t73 = t177 * t87 + t180 * t251;
t84 = -t177 * t97 - t180 * t227;
t192 = -t138 * t177 + t143 * t236 - t84;
t85 = t177 * t227 - t180 * t97;
t191 = -t138 * t180 - t143 * t237 - t85;
t15 = -pkin(4) * t283 - t19;
t24 = t175 * t53 - t262;
t187 = -t168 * t55 + (t15 + t24) * t274;
t2 = -t8 * qJD(5) + t180 * t14 - t177 * t4;
t186 = -t263 + t3 * t143 + (pkin(9) * t227 - qJD(5) * t91 - t264) * t274;
t169 = -pkin(3) * t175 - pkin(4);
t99 = -t158 * t173 - t175 * t221;
t74 = t180 * t87 - t230;
t65 = -t173 * t93 + t175 * t94;
t64 = t173 * t94 + t175 * t93;
t38 = -qJD(5) * t230 + t177 * t65 - t180 * t225 + t87 * t236;
t37 = -qJD(5) * t73 + t177 * t225 + t180 * t65;
t30 = pkin(4) * t251 - t34;
t23 = t173 * t53 + t49;
t17 = pkin(4) * t64 - pkin(9) * t65 + t207;
t10 = pkin(9) * t225 + t12;
t9 = -pkin(4) * t225 - t11;
t7 = [0, 0, 0, t179 * t182 * t210, -t246 * t210, t211 * t224, -t211 * t225, 0, -t124 * t176 - t135 * t164 + t179 * t203, -t123 * t176 - t134 * t164 + t182 * t203, t115 * t94 + t140 * t89, -t113 * t94 - t115 * t93 - t139 * t89 - t140 * t90, t94 * t283 + (-t182 * t89 + (qJD(1) * t140 + t115) * t242) * t174, -t93 * t283 + (t182 * t90 + (-qJD(1) * t139 - t113) * t242) * t174, (-t170 * t243 + t174 * t283) * t242, t184 * t283 + t135 * t113 + t126 * t90 + t124 * t139 + t101 * t93 + (-t185 * t182 + (t213 * qJD(1) + t70) * t242) * t174, -t193 * t283 + t135 * t115 + t126 * t89 + t124 * t140 + t101 * t94 + (t194 * t182 + (-t248 * qJD(1) - t71) * t242) * t174, -t11 * t197 + t12 * t214 - t19 * t65 - t20 * t64 - t34 * t56 - t35 * t55 - t5 * t87 - t6 * t86, t19 * t11 + t20 * t12 + t188 * t68 + t207 * t76 + t5 * t34 + t6 * t35, -t198 * t37 + t26 * t74, t198 * t38 - t26 * t73 - t27 * t74 - t37 * t59, -t198 * t64 + t26 * t86 + t274 * t37 + t55 * t74, -t27 * t86 - t274 * t38 - t55 * t73 - t59 * t64, t274 * t64 + t55 * t86, (-qJD(5) * t200 - t177 * t10 + t180 * t17) * t274 + t199 * t55 + t2 * t86 - t201 * t64 + t9 * t59 + t30 * t27 + t3 * t73 + t15 * t38, -(qJD(5) * t199 + t180 * t10 + t177 * t17) * t274 - t200 * t55 - t1 * t86 - t8 * t64 - t9 * t198 + t30 * t26 + t3 * t74 + t15 * t37; 0, 0, 0, -t182 * t231, t246 * t253, t235 * t226, -t235 * t227, 0, pkin(1) * t231 + t133 * t164 - t124, pkin(7) * t205 + t130 * t164 + (-t176 * t234 + t253) * t271, t115 * t254 + t89 * t178, (t89 - t258) * t181 + (-t257 - t90) * t178, t283 * t238 + (-t283 * t250 + (qJD(2) * t178 - t115) * t179) * t245, -t283 * t239 + (t182 * t255 + (t113 + t241) * t179) * t245, -t283 * t227, -pkin(2) * t90 - t124 * t181 - t212 * t283 - t133 * t113 + (-pkin(8) * t254 + t101 * t178) * qJD(3) + (-t70 * t179 + (-pkin(8) * t242 - t101 * t182) * t178) * t245, -pkin(2) * t89 + t124 * t178 + t247 * t283 - t133 * t115 + (pkin(8) * t255 + t101 * t181) * qJD(3) + (-t101 * t250 + (-pkin(8) * t241 + t71) * t179) * t245, -t142 * t6 - t143 * t5 + t279 * t19 - t266 * t197 - t259 * t20 + t264 * t214 + t56 * t99 - t263, t6 * t100 + t266 * t19 + t264 * t20 + t68 * t228 + t275 * t76 - t5 * t99, -t191 * t198 + t26 * t256, t59 * t85 - t198 * t84 - (t177 * t198 - t180 * t59) * t138 + (-t260 - t180 * t27 + (t177 * t59 + t180 * t198) * qJD(5)) * t143, t142 * t26 + t191 * t274 - t198 * t259 + t55 * t256, -t142 * t27 - t143 * t261 - t192 * t274 - t259 * t59, t55 * t142 + t259 * t274, t2 * t142 + t192 * t15 + t186 * t177 - t272 * t180 - t201 * t259 + t265 * t59 + t99 * t27, -t1 * t142 + t191 * t15 + t272 * t177 + t186 * t180 - t198 * t265 - t259 * t8 + t99 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t115 * t113, -t113 ^ 2 + t115 ^ 2, t89 + t258, -t90 + t257, t205, -t101 * t115 + t283 * t71 + t185, t101 * t113 + t283 * t70 - t194, (-t173 * t55 - t175 * t56) * pkin(3) + (t19 - t24) * t214 + (t20 - t23) * t197, t19 * t23 - t20 * t24 + (-t115 * t76 + t173 * t6 + t175 * t5) * pkin(3), -t198 * t216 + t260, (t26 - t280) * t180 + (-t27 + t278) * t177, t269 - t277, t196 + t270, -t274 * t197, t169 * t27 + t187 * t177 - t273 * t180 + t197 * t201 - t23 * t59, t169 * t26 + t273 * t177 + t187 * t180 + t8 * t197 + t198 * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t197 ^ 2 - t214 ^ 2, t19 * t197 - t20 * t214 + t68, 0, 0, 0, 0, 0, t196 - t270, t269 + t277; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t198 * t59, t198 ^ 2 - t59 ^ 2, t26 + t280, -t27 - t278, t55, t15 * t198 + t274 * t8 + t2, t15 * t59 - t201 * t274 - t1;];
tauc_reg = t7;
