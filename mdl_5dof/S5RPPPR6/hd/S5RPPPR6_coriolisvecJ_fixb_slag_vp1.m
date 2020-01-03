% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPPPR6_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR6_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR6_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR6_coriolisvecJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR6_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPPR6_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPPR6_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:47:37
% EndTime: 2019-12-31 17:47:45
% DurationCPUTime: 5.77s
% Computational Cost: add. (5458->418), mult. (14727->550), div. (0->0), fcn. (15388->8), ass. (0->204)
t195 = cos(pkin(8));
t196 = cos(pkin(7));
t285 = t195 * t196;
t156 = qJD(5) * t285 + qJD(1);
t193 = sin(pkin(8));
t198 = sin(qJ(1));
t284 = t195 * t198;
t194 = sin(pkin(7));
t200 = cos(qJ(1));
t286 = t194 * t200;
t145 = t193 * t286 + t284;
t197 = sin(qJ(5));
t199 = cos(qJ(5));
t280 = t196 * t200;
t112 = -t145 * t197 + t199 * t280;
t113 = t145 * t199 + t197 * t280;
t283 = t195 * t200;
t258 = t194 * t283;
t144 = t193 * t198 - t258;
t59 = Icges(6,5) * t113 + Icges(6,6) * t112 + Icges(6,3) * t144;
t294 = Icges(6,4) * t113;
t62 = Icges(6,2) * t112 + Icges(6,6) * t144 + t294;
t108 = Icges(6,4) * t112;
t65 = Icges(6,1) * t113 + Icges(6,5) * t144 + t108;
t14 = t112 * t62 + t113 * t65 + t144 * t59;
t146 = t193 * t200 + t194 * t284;
t287 = t194 * t198;
t147 = -t193 * t287 + t283;
t281 = t196 * t199;
t116 = t147 * t197 + t198 * t281;
t282 = t196 * t198;
t117 = -t147 * t199 + t197 * t282;
t61 = Icges(6,5) * t117 + Icges(6,6) * t116 - Icges(6,3) * t146;
t293 = Icges(6,4) * t117;
t63 = -Icges(6,2) * t116 + Icges(6,6) * t146 - t293;
t109 = Icges(6,4) * t116;
t66 = -Icges(6,1) * t117 + Icges(6,5) * t146 - t109;
t15 = -t112 * t63 - t113 * t66 + t144 * t61;
t226 = t14 * t144 - t146 * t15;
t142 = t193 * t196 * t197 + t194 * t199;
t143 = t193 * t281 - t194 * t197;
t92 = -Icges(6,5) * t143 + Icges(6,6) * t142 + Icges(6,3) * t285;
t292 = Icges(6,4) * t143;
t93 = Icges(6,2) * t142 + Icges(6,6) * t285 - t292;
t135 = Icges(6,4) * t142;
t94 = -Icges(6,1) * t143 + Icges(6,5) * t285 + t135;
t29 = t112 * t93 + t113 * t94 + t144 * t92;
t5 = qJD(5) * t226 + t29 * t156;
t268 = qJD(1) * t198;
t130 = -qJD(1) * t258 + t193 * t268;
t131 = t145 * qJD(1);
t265 = qJD(4) * t196;
t266 = qJD(3) * t194;
t317 = (t265 + t266) * t198;
t334 = -pkin(4) * t131 - pkin(6) * t130 - t317;
t332 = (-pkin(3) - qJ(2)) * t198;
t152 = pkin(3) * t200 - qJ(4) * t282;
t167 = t200 * t265;
t288 = qJ(3) * t194;
t227 = pkin(2) * t196 + t288;
t148 = t227 * t198;
t168 = t200 * t266;
t187 = t200 * qJ(2);
t157 = pkin(1) * t198 - t187;
t184 = qJD(2) * t198;
t269 = -qJD(1) * t157 + t184;
t237 = -qJD(1) * t148 + t168 + t269;
t267 = qJD(1) * t200;
t270 = qJ(2) * t267 + t184;
t253 = t168 + t270;
t273 = pkin(3) * t267 + t167;
t331 = -qJD(1) * t152 - t167 - t237 + t253 + t273;
t262 = qJD(5) * t146;
t70 = rSges(6,1) * t117 + rSges(6,2) * t116 - rSges(6,3) * t146;
t95 = -rSges(6,1) * t143 + rSges(6,2) * t142 + rSges(6,3) * t285;
t330 = -t156 * t70 - t95 * t262;
t16 = t116 * t62 + t117 * t65 - t146 * t59;
t17 = -t116 * t63 - t117 * t66 - t146 * t61;
t225 = t144 * t16 - t146 * t17;
t21 = -t142 * t63 + t143 * t66 + t285 * t61;
t151 = t198 * pkin(3) + qJ(4) * t280;
t159 = t200 * pkin(1) + t198 * qJ(2);
t236 = pkin(2) * t280 + qJ(3) * t286 + t159;
t313 = t236 + t151;
t321 = t145 * pkin(4) + pkin(6) * t144 + t313;
t320 = t145 * rSges(5,1) - t144 * rSges(5,2) + rSges(5,3) * t280 + t313;
t30 = t116 * t93 + t117 * t94 - t146 * t92;
t218 = rSges(3,1) * t280 - rSges(3,2) * t286 + t198 * rSges(3,3);
t316 = t159 + t218;
t217 = t198 * rSges(4,1) - rSges(4,2) * t280 + rSges(4,3) * t286;
t315 = t217 + t236;
t241 = -pkin(1) - t288;
t307 = -pkin(2) - qJ(4);
t314 = (-rSges(5,3) + t307) * t196 + t241;
t185 = qJD(2) * t200;
t207 = -t185 + t317;
t312 = t144 * (-Icges(6,1) * t112 + t294 + t62) - t146 * (-Icges(6,1) * t116 + t293 - t63);
t132 = t146 * qJD(1);
t250 = t196 * t267;
t76 = -qJD(5) * t117 - t131 * t197 + t199 * t250;
t77 = qJD(5) * t116 + t131 * t199 + t197 * t250;
t38 = rSges(6,1) * t77 + rSges(6,2) * t76 + rSges(6,3) * t130;
t133 = t147 * qJD(1);
t251 = t196 * t268;
t78 = -qJD(5) * t113 - t133 * t197 - t199 * t251;
t79 = qJD(5) * t112 + t133 * t199 - t197 * t251;
t39 = t79 * rSges(6,1) + t78 * rSges(6,2) + t132 * rSges(6,3);
t68 = t113 * rSges(6,1) + t112 * rSges(6,2) + t144 * rSges(6,3);
t203 = -t130 * t68 + t132 * t70 + t144 * t38 + t146 * t39;
t9 = t203 * qJD(5);
t311 = m(6) * t9;
t310 = t198 / 0.2e1;
t309 = -t200 / 0.2e1;
t308 = qJD(5) / 0.2e1;
t303 = rSges(4,1) * t200;
t302 = rSges(4,3) * t194;
t297 = -rSges(4,3) - qJ(3);
t296 = Icges(6,2) * t143 + t135 + t94;
t295 = Icges(6,1) * t142 + t292 - t93;
t141 = qJD(1) * t159 - t185;
t249 = t198 * t266;
t279 = -t227 * t267 - t141 - t249;
t278 = t133 * rSges(5,1) - t132 * rSges(5,2);
t261 = qJD(1) * qJD(2);
t277 = qJD(1) * (-pkin(1) * t268 + t270) + t198 * t261;
t276 = -t148 - t157;
t252 = t194 * t268;
t275 = rSges(3,2) * t252 + rSges(3,3) * t267;
t274 = rSges(4,1) * t267 + rSges(4,2) * t251;
t272 = t168 + t184;
t271 = rSges(3,2) * t287 + t200 * rSges(3,3);
t264 = qJD(4) * t198;
t263 = qJD(5) * t144;
t255 = t152 + t276;
t254 = t167 + t272;
t248 = -rSges(3,1) * t196 - pkin(1);
t247 = t130 * t308;
t246 = t132 * t308;
t244 = t263 / 0.2e1;
t243 = -t262 / 0.2e1;
t242 = t262 / 0.2e1;
t240 = t133 * pkin(4) + pkin(6) * t132;
t238 = t277 + (-t227 * t268 + 0.2e1 * t168) * qJD(1);
t234 = t152 + t187;
t233 = -pkin(1) + (rSges(4,2) - pkin(2)) * t196;
t231 = pkin(4) * t147 + pkin(6) * t146;
t229 = -rSges(5,1) * t131 + t130 * rSges(5,2);
t228 = rSges(5,1) * t147 - rSges(5,2) * t146;
t224 = t144 * t70 + t146 * t68;
t223 = t144 * (Icges(6,5) * t112 - Icges(6,6) * t113) - t146 * (Icges(6,5) * t116 - Icges(6,6) * t117);
t128 = t142 * qJD(5);
t129 = t143 * qJD(5);
t101 = rSges(6,1) * t128 + rSges(6,2) * t129;
t221 = -t101 * t146 + t130 * t95;
t220 = -t101 * t144 - t132 * t95;
t216 = -pkin(1) - t227;
t215 = -t249 + t279;
t214 = t238 + (-qJ(4) * t251 + t167 + t273) * qJD(1);
t209 = (-Icges(6,2) * t113 + t108 + t65) * t144 - (-Icges(6,2) * t117 + t109 - t66) * t146;
t208 = -rSges(5,3) * t282 + t228;
t32 = Icges(6,5) * t77 + Icges(6,6) * t76 + Icges(6,3) * t130;
t33 = Icges(6,5) * t79 + Icges(6,6) * t78 + Icges(6,3) * t132;
t34 = Icges(6,4) * t77 + Icges(6,2) * t76 + Icges(6,6) * t130;
t35 = Icges(6,4) * t79 + Icges(6,2) * t78 + Icges(6,6) * t132;
t36 = Icges(6,1) * t77 + Icges(6,4) * t76 + Icges(6,5) * t130;
t37 = Icges(6,1) * t79 + Icges(6,4) * t78 + Icges(6,5) * t132;
t206 = (t116 * t35 + t117 * t37 + t130 * t59 - t146 * t33 + t62 * t76 + t65 * t77) * t144 + t130 * t17 + t132 * t16 - t146 * (t116 * t34 + t117 * t36 + t130 * t61 - t146 * t32 - t63 * t76 - t66 * t77);
t205 = t130 * t15 + t132 * t14 + t144 * (t112 * t35 + t113 * t37 + t132 * t59 + t144 * t33 + t62 * t78 + t65 * t79) - t146 * (t112 * t34 + t113 * t36 + t132 * t61 + t144 * t32 - t63 * t78 - t66 * t79);
t20 = t142 * t62 - t143 * t65 + t285 * t59;
t7 = t128 * t65 + t129 * t62 + t142 * t35 - t143 * t37 + t285 * t33;
t8 = -t128 * t66 - t129 * t63 + t142 * t34 - t143 * t36 + t285 * t32;
t204 = t130 * t21 + t132 * t20 + t144 * t7 - t146 * t8;
t121 = qJD(1) * t151 + t196 * t264;
t179 = t200 * t261;
t12 = -t156 * t38 + t179 + t221 * qJD(5) + (-t121 + t279 + t334) * qJD(1);
t13 = qJD(1) * t240 + qJD(5) * t220 + t156 * t39 + t214;
t40 = qJD(1) * (-rSges(5,3) * t251 + t278) + t214;
t41 = t179 + (-t121 + (-rSges(5,3) * t267 - t264) * t196 + t215 + t229) * qJD(1);
t202 = m(5) * (t198 * t40 + t200 * t41) / 0.2e1 + m(6) * (t12 * t200 + t13 * t198) / 0.2e1;
t127 = t303 + (rSges(4,2) * t196 - t302) * t198;
t126 = rSges(3,1) * t282 - t271;
t104 = Icges(6,5) * t142 + Icges(6,6) * t143;
t103 = qJD(1) * t316 - t185;
t102 = t184 + (-t126 - t157) * qJD(1);
t100 = Icges(6,1) * t128 + Icges(6,4) * t129;
t99 = Icges(6,4) * t128 + Icges(6,2) * t129;
t98 = Icges(6,5) * t128 + Icges(6,6) * t129;
t91 = t179 + (-qJD(1) * t218 - t141) * qJD(1);
t90 = qJD(1) * (-rSges(3,1) * t251 + t275) + t277;
t89 = rSges(6,1) * t116 - rSges(6,2) * t117;
t88 = rSges(6,1) * t112 - rSges(6,2) * t113;
t81 = qJD(1) * t315 - t185 + t249;
t80 = (t127 + t276) * qJD(1) + t272;
t58 = t179 + (-qJD(1) * t217 + t215) * qJD(1);
t57 = qJD(1) * (-rSges(4,3) * t252 + t274) + t238;
t52 = (t208 + t255) * qJD(1) + t254;
t31 = -qJD(3) * t196 + qJD(4) * t194 + qJD(5) * t224;
t27 = t321 * qJD(1) + t156 * t68 - t95 * t263 + t207;
t26 = (t231 + t255) * qJD(1) + t254 + t330;
t19 = -t100 * t143 + t128 * t94 + t129 * t93 + t142 * t99 + t285 * t98;
t18 = t19 * t156;
t11 = t100 * t113 + t112 * t99 + t132 * t92 + t144 * t98 + t78 * t93 + t79 * t94;
t10 = t100 * t117 + t116 * t99 + t130 * t92 - t146 * t98 + t76 * t93 + t77 * t94;
t1 = [t18 + t5 * t242 + (t21 + t30) * t247 + (t20 + t29) * t246 + (t7 + t11) * t244 + (t13 * (t68 + t321) + (t216 * t198 + t231 + t234 - t70) * t12 + (t185 - t38 + ((t196 * t307 + t241) * t200 + t332) * qJD(1) + t334) * t26 + (t240 + t39 + t26 + ((-qJ(4) * t196 + t216) * t198 - t231) * qJD(1) - t330 + t331) * t27) * m(6) + (t40 * t320 + (t228 + t234 + (-rSges(5,3) * t196 + t216) * t198) * t41 + (t229 + (t200 * t314 + t332) * qJD(1) - t207) * t52 + (t278 + t52 + (t198 * t314 - t208) * qJD(1) + t331) * (t320 * qJD(1) + t207)) * m(5) + (t58 * (t187 + t303) + t80 * t185 + t57 * t315 + t81 * (t253 + t274) + (t58 * t233 + (-t80 * qJD(3) + t297 * t58) * t194) * t198 + (t80 * (t194 * t297 + t233) * t200 + (t80 * (-rSges(4,1) - qJ(2)) + t81 * (t216 - t302)) * t198) * qJD(1) - (qJD(1) * t127 + t237 - t80) * t81) * m(4) + (t91 * (t198 * t248 + t187 + t271) + t102 * t185 + t90 * t316 + t103 * (t270 + t275) + (t102 * (rSges(3,2) * t194 + t248) * t200 + (t102 * (-rSges(3,3) - qJ(2)) + t103 * t248) * t198) * qJD(1) - (-qJD(1) * t126 - t102 + t269) * t103) * m(3) + (t8 + t10 + t5) * t243; 0.2e1 * (t12 * t310 + t13 * t309) * m(6) + 0.2e1 * (t309 * t40 + t310 * t41) * m(5) + 0.2e1 * (t309 * t57 + t310 * t58) * m(4) + 0.2e1 * (t309 * t90 + t310 * t91) * m(3); -t196 * t311 + 0.2e1 * (m(4) * (t198 * t57 + t200 * t58) / 0.2e1 + t202) * t194; t194 * t311 + 0.2e1 * t196 * t202; t132 * t5 / 0.2e1 + t144 * (t205 * qJD(5) + t11 * t156) / 0.2e1 + (t285 * t29 + t226) * t246 + (t11 * t285 + t205) * t244 + t130 * (qJD(5) * t225 + t156 * t30) / 0.2e1 - t146 * (t206 * qJD(5) + t10 * t156) / 0.2e1 + (t285 * t30 + t225) * t247 + (t10 * t285 + t206) * t243 + (t204 * qJD(5) + t18) * t285 / 0.2e1 + t156 * (t19 * t285 + t204) / 0.2e1 - ((t144 * t104 + t112 * t296 + t113 * t295) * t156 + (t112 * t209 - t113 * t312 + t144 * t223) * qJD(5)) * t263 / 0.2e1 + ((-t146 * t104 + t116 * t296 + t117 * t295) * t156 + (t116 * t209 - t117 * t312 - t146 * t223) * qJD(5)) * t242 - t156 * ((t104 * t285 + t296 * t142 - t295 * t143) * t156 + (t209 * t142 + t143 * t312 + t223 * t285) * qJD(5)) / 0.2e1 + (t12 * (-t146 * t95 - t285 * t70) + t26 * (-t38 * t285 + t221) + t13 * (-t144 * t95 + t285 * t68) + t27 * (t39 * t285 + t220) + t9 * t224 + t31 * t203 - (-t26 * t89 + t27 * t88) * t156 - (t31 * (t144 * t89 + t146 * t88) + (-t144 * t27 - t146 * t26) * (rSges(6,1) * t142 + rSges(6,2) * t143)) * qJD(5)) * m(6);];
tauc = t1(:);
