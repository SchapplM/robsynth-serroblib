% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RRRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRRRP9_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP9_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP9_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP9_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:06:04
% EndTime: 2019-12-31 22:06:16
% DurationCPUTime: 4.49s
% Computational Cost: add. (5877->439), mult. (14571->574), div. (0->0), fcn. (9658->6), ass. (0->210)
t180 = cos(qJ(2));
t240 = t180 * qJD(1);
t160 = -qJD(3) + t240;
t147 = -qJD(4) + t160;
t177 = sin(qJ(3));
t178 = sin(qJ(2));
t249 = qJD(1) * t178;
t222 = t177 * t249;
t179 = cos(qJ(3));
t246 = qJD(2) * t179;
t123 = -t222 + t246;
t227 = t179 * t249;
t248 = qJD(2) * t177;
t124 = t227 + t248;
t176 = sin(qJ(4));
t292 = cos(qJ(4));
t187 = t176 * t123 + t292 * t124;
t243 = qJD(3) * t178;
t218 = qJD(1) * t243;
t239 = qJD(1) * qJD(2);
t219 = t180 * t239;
t308 = qJD(2) * qJD(3) + t219;
t231 = t308 * t177 + t179 * t218;
t93 = t177 * t218 - t308 * t179;
t32 = qJD(4) * t187 - t176 * t93 + t292 * t231;
t310 = t147 * t187 + t32;
t75 = -t292 * t123 + t124 * t176;
t170 = pkin(6) * t249;
t280 = qJD(2) * pkin(2);
t139 = t170 - t280;
t96 = -pkin(3) * t123 + t139;
t27 = pkin(4) * t75 - qJ(5) * t187 + t96;
t313 = t27 * t75;
t312 = t75 * t96;
t311 = t75 * t187;
t260 = t176 * t179;
t126 = t292 * t177 + t260;
t298 = qJD(3) + qJD(4);
t86 = t298 * t126;
t273 = t126 * t240 - t86;
t221 = t292 * qJD(3);
t205 = t179 * t221;
t220 = t292 * qJD(4);
t228 = t177 * t240;
t229 = t292 * t179;
t261 = t176 * t177;
t272 = -t176 * t228 - t179 * t220 + t229 * t240 + t298 * t261 - t205;
t242 = qJD(3) * t179;
t224 = t178 * t242;
t245 = qJD(2) * t180;
t226 = t177 * t245;
t307 = t224 + t226;
t294 = t187 ^ 2;
t235 = t75 ^ 2 - t294;
t241 = qJD(4) * t176;
t31 = -t123 * t220 + t124 * t241 + t176 * t231 + t292 * t93;
t306 = -t147 * t75 - t31;
t43 = pkin(4) * t187 + qJ(5) * t75;
t305 = -0.2e1 * t239;
t293 = -pkin(8) - pkin(7);
t230 = qJD(3) * t293;
t129 = t177 * t230;
t256 = t179 * t180;
t193 = pkin(3) * t178 - pkin(8) * t256;
t200 = pkin(2) * t178 - pkin(7) * t180;
t128 = t200 * qJD(1);
t94 = pkin(6) * t222 + t179 * t128;
t68 = qJD(1) * t193 + t94;
t112 = t177 * t128;
t257 = t178 * t179;
t258 = t177 * t180;
t80 = t112 + (-pkin(6) * t257 - pkin(8) * t258) * qJD(1);
t141 = t293 * t177;
t142 = t293 * t179;
t92 = t176 * t141 - t292 * t142;
t283 = t92 * qJD(4) - t293 * t205 + t292 * t68 + (t129 - t80) * t176;
t290 = t27 * t187;
t302 = t96 * t187;
t133 = -pkin(2) * t180 - pkin(7) * t178 - pkin(1);
t117 = t133 * qJD(1);
t130 = t200 * qJD(2);
t118 = qJD(1) * t130;
t171 = pkin(6) * t240;
t140 = qJD(2) * pkin(7) + t171;
t164 = t178 * t239;
t208 = pkin(6) * t164;
t244 = qJD(3) * t177;
t49 = t117 * t242 + t177 * t118 - t140 * t244 - t179 * t208;
t82 = t179 * t117 - t140 * t177;
t301 = t160 * t82 + t49;
t83 = t117 * t177 + t140 * t179;
t50 = -qJD(3) * t83 + t179 * t118 + t177 * t208;
t300 = t160 * t83 - t50;
t202 = -t171 + (-t228 + t244) * pkin(3);
t162 = pkin(6) * t256;
t100 = t177 * t133 + t162;
t107 = t126 * t178;
t247 = qJD(2) * t178;
t206 = t292 * t245;
t225 = t177 * t243;
t259 = t177 * t178;
t233 = t176 * t259;
t52 = t177 * t206 - t176 * t225 - qJD(4) * t233 + (t176 * t245 + (t221 + t220) * t178) * t179;
t297 = -t147 * t52 - t180 * t32 + (qJD(1) * t107 + t75) * t247;
t122 = t179 * t133;
t291 = pkin(6) * t177;
t81 = -pkin(8) * t257 + t122 + (-pkin(3) - t291) * t180;
t87 = -pkin(8) * t259 + t100;
t281 = t176 * t81 + t292 * t87;
t251 = t179 * t130 + t247 * t291;
t44 = t193 * qJD(2) + (-t162 + (pkin(8) * t178 - t133) * t177) * qJD(3) + t251;
t63 = t177 * t130 + t133 * t242 + (-t178 * t246 - t180 * t244) * pkin(6);
t48 = -pkin(8) * t307 + t63;
t11 = -qJD(4) * t281 - t176 * t48 + t292 * t44;
t125 = -t229 + t261;
t296 = t273 * t147 + (qJD(2) * t125 - t75) * t249;
t295 = t125 * t31 - t126 * t32 + t187 * t273 + t272 * t75;
t286 = -t273 * pkin(4) + t272 * qJ(5) - qJD(5) * t126 + t202;
t39 = t176 * t68 + t292 * t80;
t36 = qJ(5) * t249 + t39;
t186 = t292 * t141 + t176 * t142;
t56 = t186 * qJD(4) + t292 * t129 + t230 * t260;
t285 = t36 - t56;
t284 = pkin(4) * t249 + t283;
t282 = -t39 + t56;
t66 = pkin(8) * t123 + t83;
t276 = t176 * t66;
t274 = t93 * t177;
t65 = -pkin(8) * t124 + t82;
t23 = t292 * t65 - t276;
t271 = pkin(3) * t220 + qJD(5) - t23;
t270 = qJD(2) * t186;
t269 = qJD(2) * t92;
t268 = t123 * t160;
t267 = t124 * t123;
t266 = t124 * t160;
t265 = t139 * t177;
t264 = t139 * t179;
t263 = t160 * t177;
t262 = t160 * t179;
t182 = qJD(1) ^ 2;
t255 = t180 * t182;
t181 = qJD(2) ^ 2;
t254 = t181 * t178;
t253 = t181 * t180;
t58 = -pkin(3) * t160 + t65;
t20 = t292 * t58 - t276;
t252 = qJD(5) - t20;
t131 = pkin(3) * t259 + t178 * pkin(6);
t174 = t178 ^ 2;
t250 = -t180 ^ 2 + t174;
t237 = pkin(6) * t258;
t172 = pkin(6) * t245;
t234 = t292 * t66;
t232 = t178 * t255;
t97 = t307 * pkin(3) + t172;
t169 = -pkin(3) * t179 - pkin(2);
t223 = t147 * t249;
t26 = pkin(3) * t164 + pkin(8) * t93 + t50;
t33 = -t231 * pkin(8) + t49;
t216 = -t176 * t26 - t58 * t220 + t66 * t241 - t292 * t33;
t215 = t176 * t33 + t66 * t220 + t58 * t241 - t292 * t26;
t214 = pkin(1) * t305;
t211 = -t123 + t246;
t210 = -t124 + t248;
t209 = pkin(4) * t164;
t207 = t186 * t31 - t92 * t32 - t56 * t75;
t204 = t180 * t164;
t22 = t176 * t65 + t234;
t203 = pkin(3) * t241 - t22;
t201 = t231 * t179;
t198 = t107 * t32 + t52 * t75;
t195 = -t177 * t83 - t179 * t82;
t194 = qJD(1) * t174 - t160 * t180;
t132 = t147 * qJD(5);
t155 = qJ(5) * t164;
t1 = t155 - t132 - t216;
t73 = t231 * pkin(3) + pkin(6) * t219;
t192 = -t147 * t20 + t216;
t21 = t176 * t58 + t234;
t191 = -t147 * t21 - t215;
t45 = -t176 * t87 + t292 * t81;
t10 = t176 * t44 + t81 * t220 - t87 * t241 + t292 * t48;
t188 = t125 * t32 - t273 * t75;
t2 = -t209 + t215;
t108 = t178 * t229 - t233;
t51 = t176 * t226 + t178 * t86 - t179 * t206;
t184 = t107 * t31 - t108 * t32 - t187 * t52 + t51 * t75;
t168 = -t292 * pkin(3) - pkin(4);
t163 = pkin(3) * t176 + qJ(5);
t99 = t122 - t237;
t98 = (-t147 - t240) * t247;
t95 = -pkin(6) * t227 + t112;
t72 = pkin(4) * t125 - qJ(5) * t126 + t169;
t64 = -t100 * qJD(3) + t251;
t59 = pkin(4) * t107 - qJ(5) * t108 + t131;
t42 = t180 * pkin(4) - t45;
t41 = -qJ(5) * t180 + t281;
t35 = pkin(3) * t124 + t43;
t18 = -t147 * qJ(5) + t21;
t17 = t147 * pkin(4) + t252;
t16 = t272 * t147 + (qJD(2) * t126 - t187) * t249;
t13 = pkin(4) * t52 + qJ(5) * t51 - qJD(5) * t108 + t97;
t12 = -t108 * t31 - t187 * t51;
t9 = -pkin(4) * t247 - t11;
t8 = qJ(5) * t247 - qJD(5) * t180 + t10;
t7 = -t126 * t31 - t187 * t272;
t6 = t147 * t51 + t180 * t31 + (qJD(1) * t108 + t187) * t247;
t5 = t32 * pkin(4) + t31 * qJ(5) - qJD(5) * t187 + t73;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t204, t250 * t305, t253, -0.2e1 * t204, -t254, 0, -pkin(6) * t253 + t178 * t214, pkin(6) * t254 + t180 * t214, 0, 0, -t93 * t257 + (t179 * t245 - t225) * t124, (t123 * t179 - t124 * t177) * t245 + (-t201 + t274 + (-t123 * t177 - t124 * t179) * qJD(3)) * t178, t160 * t225 + t180 * t93 + (t124 * t178 + t179 * t194) * qJD(2), -t123 * t307 + t231 * t259, t160 * t224 + t231 * t180 + (t123 * t178 - t177 * t194) * qJD(2), (-t160 - t240) * t247, -t64 * t160 - t50 * t180 + (pkin(6) * t231 + t139 * t242) * t178 + ((-pkin(6) * t123 + t265) * t180 + (t82 + (t99 + t237) * qJD(1)) * t178) * qJD(2), t160 * t63 + t180 * t49 + (-pkin(6) * t93 - t139 * t244) * t178 + ((pkin(6) * t124 + t264) * t180 + (-t83 + (-t100 + t162) * qJD(1)) * t178) * qJD(2), -t64 * t124 + t99 * t93 + t63 * t123 - t100 * t231 + t195 * t245 + (-t177 * t49 - t179 * t50 + (t177 * t82 - t179 * t83) * qJD(3)) * t178, t100 * t49 + t50 * t99 + t63 * t83 + t64 * t82 + (t139 + t170) * t172, t12, t184, t6, t198, -t297, t98, t107 * t73 - t11 * t147 + t131 * t32 + t180 * t215 + t52 * t96 + t75 * t97 + (qJD(1) * t45 + t20) * t247, t10 * t147 + t108 * t73 - t131 * t31 - t180 * t216 - t51 * t96 + t187 * t97 + (-qJD(1) * t281 - t21) * t247, -t10 * t75 + t107 * t216 + t108 * t215 - t11 * t187 + t20 * t51 - t21 * t52 - t281 * t32 + t31 * t45, t10 * t21 + t11 * t20 + t131 * t73 - t215 * t45 - t216 * t281 + t96 * t97, t12, t6, -t184, t98, t297, t198, t107 * t5 + t13 * t75 + t147 * t9 + t180 * t2 + t27 * t52 + t32 * t59 + (-qJD(1) * t42 - t17) * t247, -t1 * t107 + t108 * t2 - t17 * t51 - t18 * t52 + t187 * t9 - t31 * t42 - t32 * t41 - t75 * t8, -t1 * t180 - t108 * t5 - t13 * t187 - t147 * t8 + t27 * t51 + t31 * t59 + (qJD(1) * t41 + t18) * t247, t1 * t41 + t13 * t27 + t17 * t9 + t18 * t8 + t2 * t42 + t5 * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t232, t250 * t182, 0, t232, 0, 0, t182 * pkin(1) * t178, pkin(1) * t255, 0, 0, -t124 * t262 - t274, (-t93 - t268) * t179 + (-t231 + t266) * t177, -t160 * t242 + (t160 * t256 + t178 * t210) * qJD(1), t123 * t263 - t201, t160 * t244 + (-t160 * t258 + t178 * t211) * qJD(1), t160 * t249, -pkin(2) * t231 + t94 * t160 + (pkin(7) * t262 + t265) * qJD(3) + ((-pkin(7) * t248 - t82) * t178 + (-pkin(6) * t211 - t265) * t180) * qJD(1), pkin(2) * t93 - t160 * t95 + (-pkin(7) * t263 + t264) * qJD(3) + ((-pkin(7) * t246 + t83) * t178 + (pkin(6) * t210 - t264) * t180) * qJD(1), -t95 * t123 + t94 * t124 + ((qJD(3) * t124 - t231) * pkin(7) + t301) * t179 + ((-qJD(3) * t123 - t93) * pkin(7) + t300) * t177, -t82 * t94 - t83 * t95 + (-t139 - t280) * t171 + (qJD(3) * t195 - t50 * t177 + t49 * t179) * pkin(7), t7, t295, t16, t188, -t296, t223, t125 * t73 + t169 * t32 - t273 * t96 + t202 * t75 + t283 * t147 + (-t20 + t270) * t249, t126 * t73 - t169 * t31 - t272 * t96 + t202 * t187 + t282 * t147 + (t21 - t269) * t249, t125 * t216 + t126 * t215 + t187 * t283 + t272 * t20 + t273 * t21 + t39 * t75 + t207, t169 * t73 - t186 * t215 - t283 * t20 + t202 * t96 + t282 * t21 - t216 * t92, t7, t16, -t295, t223, t296, t188, t125 * t5 + t32 * t72 + t286 * t75 - t273 * t27 + t284 * t147 + (t17 + t270) * t249, -t1 * t125 + t126 * t2 - t272 * t17 + t273 * t18 + t187 * t284 + t36 * t75 + t207, -t126 * t5 + t31 * t72 - t286 * t187 + t272 * t27 + t285 * t147 + (-t18 + t269) * t249, t1 * t92 + t284 * t17 - t285 * t18 - t186 * t2 + t286 * t27 + t5 * t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t267, -t123 ^ 2 + t124 ^ 2, -t93 + t268, t267, -t231 - t266, t164, -t124 * t139 - t300, -t123 * t139 - t301, 0, 0, t311, -t235, t306, -t311, -t310, t164, -t22 * t147 - t302 + (-t124 * t75 + t147 * t241 + t292 * t164) * pkin(3) - t215, -t23 * t147 + t312 + (-t124 * t187 + t147 * t220 - t164 * t176) * pkin(3) + t216, -t20 * t75 + t21 * t187 - t22 * t187 + t23 * t75 + (t292 * t31 - t176 * t32 + (t176 * t187 - t292 * t75) * qJD(4)) * pkin(3), t20 * t22 - t21 * t23 + (-t292 * t215 - t124 * t96 - t176 * t216 + (-t176 * t20 + t292 * t21) * qJD(4)) * pkin(3), t311, t306, t235, t164, t310, -t311, -t290 - t35 * t75 + t203 * t147 + (pkin(4) - t168) * t164 - t215, -t163 * t32 - t168 * t31 + (t18 + t203) * t187 + (t17 - t271) * t75, -t147 * t271 + t163 * t164 + t187 * t35 + t1 - t313, t1 * t163 + t168 * t2 + t17 * t203 + t18 * t271 - t27 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t311, -t235, t306, -t311, -t310, t164, t191 - t302, t192 + t312, 0, 0, t311, t306, t235, t164, t310, -t311, -t43 * t75 + t191 + 0.2e1 * t209 - t290, pkin(4) * t31 - qJ(5) * t32 + (t18 - t21) * t187 + (t17 - t252) * t75, t187 * t43 - 0.2e1 * t132 + 0.2e1 * t155 - t192 - t313, -pkin(4) * t2 + qJ(5) * t1 - t17 * t21 + t18 * t252 - t27 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t164 + t311, t306, -t147 ^ 2 - t294, t147 * t18 + t2 + t290;];
tauc_reg = t3;
