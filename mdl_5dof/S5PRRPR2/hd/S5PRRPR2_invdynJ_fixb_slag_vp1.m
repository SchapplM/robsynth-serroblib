% Calculate vector of inverse dynamics joint torques for
% S5PRRPR2
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:18
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRRPR2_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR2_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR2_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR2_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR2_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR2_invdynJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR2_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPR2_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRPR2_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:17:19
% EndTime: 2019-12-05 16:17:28
% DurationCPUTime: 5.47s
% Computational Cost: add. (11786->410), mult. (9684->558), div. (0->0), fcn. (9248->8), ass. (0->207)
t190 = pkin(8) + qJ(2);
t188 = qJ(3) + t190;
t183 = sin(t188);
t184 = cos(t188);
t148 = t184 * pkin(3) + t183 * qJ(4);
t174 = qJD(4) * t184;
t191 = qJD(2) + qJD(3);
t107 = t148 * t191 - t174;
t192 = sin(pkin(9));
t247 = qJD(5) * t191;
t114 = (qJDD(5) * t183 + t184 * t247) * t192;
t193 = cos(pkin(9));
t194 = sin(qJ(5));
t195 = cos(qJ(5));
t139 = -rSges(6,3) * t193 + (rSges(6,1) * t195 - rSges(6,2) * t194) * t192;
t154 = (-rSges(6,1) * t194 - rSges(6,2) * t195) * t192;
t143 = qJD(5) * t154;
t189 = qJDD(2) + qJDD(3);
t170 = -qJDD(5) * t193 + t189;
t171 = -qJD(5) * t193 + t191;
t186 = sin(t190);
t187 = cos(t190);
t196 = qJD(2) ^ 2;
t207 = (-qJDD(2) * t186 - t187 * t196) * pkin(2);
t248 = qJD(4) * t191;
t200 = qJDD(4) * t183 + t184 * t248 + t207;
t281 = pkin(4) * t193;
t222 = pkin(7) * t192 + t281;
t246 = qJD(5) * t192;
t236 = t183 * t246;
t132 = t222 * t183;
t177 = t184 * qJ(4);
t146 = pkin(3) * t183 - t177;
t253 = -t146 - t132;
t264 = t184 * t191;
t259 = t193 * t195;
t129 = t183 * t259 - t184 * t194;
t260 = t193 * t194;
t130 = t183 * t195 - t184 * t260;
t90 = -qJD(5) * t129 + t130 * t191;
t128 = t183 * t260 + t184 * t195;
t131 = t183 * t194 + t184 * t259;
t91 = -qJD(5) * t128 + t131 * t191;
t223 = rSges(6,1) * t91 + rSges(6,2) * t90;
t261 = t191 * t192;
t238 = t184 * t261;
t51 = rSges(6,3) * t238 + t223;
t258 = -t129 * rSges(6,1) + t128 * rSges(6,2);
t265 = t183 * t192;
t83 = rSges(6,3) * t265 - t258;
t16 = t143 * t236 + t114 * t139 - t170 * t83 - t171 * t51 + t253 * t189 + (-t222 * t264 - t107) * t191 + t200;
t310 = t16 - g(1);
t309 = -t139 * t236 + t171 * t83;
t74 = Icges(6,5) * t129 - Icges(6,6) * t128 + Icges(6,3) * t265;
t118 = Icges(6,4) * t129;
t77 = -Icges(6,2) * t128 + Icges(6,6) * t265 + t118;
t117 = Icges(6,4) * t128;
t81 = -Icges(6,1) * t129 - Icges(6,5) * t265 + t117;
t30 = -(t194 * t77 + t195 * t81) * t192 - t193 * t74;
t307 = t177 + t258;
t306 = -t128 * t77 - t129 * t81;
t305 = t130 * t77 - t131 * t81;
t263 = t184 * t192;
t302 = t74 * t263;
t275 = pkin(2) * qJD(2);
t243 = t186 * t275;
t147 = rSges(4,1) * t183 + rSges(4,2) * t184;
t267 = t147 * t191;
t112 = -t243 - t267;
t115 = (qJDD(5) * t184 - t183 * t247) * t192;
t262 = t184 * t193;
t133 = pkin(4) * t262 + pkin(7) * t263;
t182 = pkin(2) * t187;
t282 = pkin(2) * t186;
t227 = qJDD(2) * t182 - t196 * t282;
t161 = qJ(4) * t264;
t173 = qJD(4) * t183;
t251 = t161 + t173;
t266 = t183 * t191;
t205 = t191 * (-pkin(3) * t266 + t251) + t189 * t148 + t183 * t248 + t227;
t239 = t183 * t261;
t88 = -qJD(5) * t131 + t128 * t191;
t89 = qJD(5) * t130 - t129 * t191;
t277 = t89 * rSges(6,1) + t88 * rSges(6,2);
t50 = -rSges(6,3) * t239 + t277;
t85 = t131 * rSges(6,1) + t130 * rSges(6,2) + rSges(6,3) * t263;
t17 = -t115 * t139 + t189 * t133 + t170 * t85 + t171 * t50 + (-t143 * t246 - qJDD(4)) * t184 - t191 ^ 2 * t132 + t205;
t301 = t17 - g(2);
t156 = rSges(5,2) * t238;
t276 = rSges(5,1) * t193;
t245 = t183 * t276;
t250 = rSges(5,2) * t265 + t184 * rSges(5,3);
t110 = t245 - t250;
t254 = -t146 - t110;
t293 = -rSges(5,1) * t262 - t183 * rSges(5,3);
t34 = t254 * t189 + (t191 * t293 - t107 + t156) * t191 + t200;
t300 = t34 - g(1);
t111 = -rSges(5,2) * t263 - t293;
t252 = rSges(5,2) * t239 + rSges(5,3) * t264;
t35 = -qJDD(4) * t184 + t189 * t111 + t191 * (-t191 * t245 + t252) + t205;
t299 = t35 - g(2);
t125 = rSges(4,1) * t264 - rSges(4,2) * t266;
t298 = -t125 * t191 - t147 * t189 - g(1) + t207;
t149 = t184 * rSges(4,1) - rSges(4,2) * t183;
t297 = t149 * t189 - t191 * t267 - g(2) + t227;
t26 = t302 + t305;
t296 = t26 - t302;
t103 = t111 + t148;
t295 = t103 * t191;
t294 = t191 * t110 + t252;
t138 = t191 * t146;
t292 = t161 + t138;
t291 = t251 - t173 + t138;
t76 = Icges(6,5) * t131 + Icges(6,6) * t130 + Icges(6,3) * t263;
t270 = Icges(6,4) * t131;
t79 = Icges(6,2) * t130 + Icges(6,6) * t263 + t270;
t119 = Icges(6,4) * t130;
t82 = Icges(6,1) * t131 + Icges(6,5) * t263 + t119;
t25 = -t128 * t79 + t129 * t82 + t76 * t265;
t206 = -t281 - pkin(3) + (-rSges(6,3) - pkin(7)) * t192;
t226 = t173 - t243;
t36 = t191 * t253 + t226 - t309;
t242 = t187 * t275;
t257 = t133 + t148;
t289 = -t184 * t139 * t246 + t171 * t85 + t191 * t257 - t174;
t37 = t242 + t289;
t290 = (t36 * t206 * t184 + (-t36 * qJ(4) + t37 * (-rSges(6,3) * t192 - pkin(3) - t222)) * t183) * t191 + t310 * t206 * t183;
t288 = t191 * t132 + t277 + t309;
t287 = t183 * (-Icges(6,2) * t129 - t117 - t81) + t184 * (-Icges(6,2) * t131 + t119 + t82);
t286 = t114 / 0.2e1;
t285 = t115 / 0.2e1;
t284 = t183 / 0.2e1;
t283 = -t184 / 0.2e1;
t268 = Icges(6,4) * t195;
t136 = -Icges(6,6) * t193 + (-Icges(6,2) * t194 + t268) * t192;
t269 = Icges(6,4) * t194;
t137 = -Icges(6,5) * t193 + (Icges(6,1) * t195 - t269) * t192;
t151 = (-Icges(6,5) * t194 - Icges(6,6) * t195) * t192;
t140 = qJD(5) * t151;
t152 = (-Icges(6,2) * t195 - t269) * t192;
t141 = qJD(5) * t152;
t153 = (-Icges(6,1) * t194 - t268) * t192;
t142 = qJD(5) * t153;
t42 = -t140 * t193 + (-t141 * t194 + t142 * t195 + (-t136 * t195 - t137 * t194) * qJD(5)) * t192;
t135 = -Icges(6,3) * t193 + (Icges(6,5) * t195 - Icges(6,6) * t194) * t192;
t65 = -t135 * t193 + (-t136 * t194 + t137 * t195) * t192;
t280 = t65 * t170 + t42 * t171;
t52 = -t128 * t136 + t129 * t137 + t135 * t265;
t274 = t171 * t52;
t24 = t265 * t74 + t306;
t273 = t183 * t24;
t272 = t30 * t114;
t31 = -t193 * t76 + (-t194 * t79 + t195 * t82) * t192;
t271 = t31 * t115;
t256 = -t136 + t153;
t255 = t137 + t152;
t27 = t130 * t79 + t131 * t82 + t76 * t263;
t240 = m(2) + m(3) + m(4) + m(5);
t234 = -pkin(3) - t276;
t233 = -t246 / 0.2e1;
t232 = t246 / 0.2e1;
t231 = t183 * t233;
t230 = t183 * t232;
t229 = t184 * t233;
t228 = t184 * t232;
t225 = -t174 + t242;
t158 = rSges(3,1) * t187 - rSges(3,2) * t186;
t157 = rSges(3,1) * t186 + rSges(3,2) * t187;
t220 = t36 * t183 - t184 * t37;
t219 = -t183 * t50 + t184 * t51;
t218 = -t183 * t85 + t184 * t83;
t217 = t183 * (-Icges(6,5) * t128 - Icges(6,6) * t129) + t184 * (Icges(6,5) * t130 - Icges(6,6) * t131);
t214 = t174 - t223;
t210 = t192 * t217;
t209 = (t184 * t25 + t273) * t192;
t208 = (t183 * t26 + t184 * t27) * t192;
t55 = t85 + t257;
t204 = (Icges(6,1) * t130 - t270 - t79) * t184 + (-Icges(6,1) * t128 - t118 - t77) * t183;
t102 = t183 * t234 + t177 + t250;
t53 = t130 * t136 + t131 * t137 + t135 * t263;
t43 = t53 * t171;
t10 = qJD(5) * t208 + t43;
t45 = Icges(6,5) * t91 + Icges(6,6) * t90 + Icges(6,3) * t238;
t47 = Icges(6,4) * t91 + Icges(6,2) * t90 + Icges(6,6) * t238;
t49 = Icges(6,1) * t91 + Icges(6,4) * t90 + Icges(6,5) * t238;
t13 = -t193 * t45 + (-t194 * t47 + t195 * t49 + (t194 * t81 - t195 * t77) * qJD(5)) * t192;
t44 = Icges(6,5) * t89 + Icges(6,6) * t88 - Icges(6,3) * t239;
t46 = Icges(6,4) * t89 + Icges(6,2) * t88 - Icges(6,6) * t239;
t48 = Icges(6,1) * t89 + Icges(6,4) * t88 - Icges(6,5) * t239;
t14 = -t193 * t44 + (-t194 * t46 + t195 * t48 + (-t194 * t82 - t195 * t79) * qJD(5)) * t192;
t21 = t130 * t141 + t131 * t142 + t136 * t88 + t137 * t89 + (-t135 * t266 + t140 * t184) * t192;
t22 = -t128 * t141 + t129 * t142 + t136 * t90 + t137 * t91 + (t135 * t264 + t140 * t183) * t192;
t9 = qJD(5) * t209 + t274;
t199 = (t43 + ((t24 + t27 - t306) * t184 + t296 * t183) * t246) * t231 + t272 / 0.2e1 + t271 / 0.2e1 + t52 * t286 + t53 * t285 + t280 + (t14 + t21) * t228 + (Icges(5,2) * t193 ^ 2 + (Icges(5,1) * t192 + 0.2e1 * Icges(5,4) * t193) * t192 + Icges(4,3)) * t189 + (t13 + t22 + t10) * t230 + (-t274 + ((-t25 + t296 - t305) * t184 - t273) * t246 + t9) * t229;
t66 = t191 * t254 + t226;
t67 = t225 + t295;
t198 = (t66 * t234 * t184 + (t66 * (-rSges(5,3) - qJ(4)) + t67 * t234) * t183) * t191;
t113 = t149 * t191 + t242;
t101 = rSges(6,1) * t130 - rSges(6,2) * t131;
t100 = -rSges(6,1) * t128 - rSges(6,2) * t129;
t40 = t218 * t246 + qJD(1);
t15 = -t114 * t85 + t115 * t83 + t219 * t246 + qJDD(1);
t6 = -t128 * t46 + t129 * t48 + t79 * t90 + t82 * t91 + (t183 * t44 + t264 * t76) * t192;
t5 = -t128 * t47 + t129 * t49 + t77 * t90 - t81 * t91 + (t183 * t45 + t264 * t74) * t192;
t4 = t130 * t46 + t131 * t48 + t79 * t88 + t82 * t89 + (t184 * t44 - t266 * t76) * t192;
t3 = t130 * t47 + t131 * t49 + t77 * t88 - t81 * t89 + (t184 * t45 - t266 * t74) * t192;
t1 = [m(6) * t15 + t240 * qJDD(1) + (-m(6) - t240) * g(3); Icges(3,3) * qJDD(2) + t199 + (t297 * (t149 + t182) + t298 * (-t147 - t282) + (-t125 - t242 + t113) * t112) * m(4) + ((t157 ^ 2 + t158 ^ 2) * qJDD(2) + g(1) * t157 - g(2) * t158) * m(3) + (t36 * (t214 - t242) + t301 * (t182 + t55) + (t36 + t288 + t292) * t37 + t290 + t310 * (t307 - t282)) * m(6) + (t66 * (t156 - t225) + t198 + t299 * (t103 + t182) + t300 * (t102 - t282) + (t66 + t292 + t294) * t67) * m(5); t199 + (t301 * t55 + (t288 + t291) * t37 + (t214 + t289) * t36 + t290 + t310 * t307) * m(6) + (t198 + (t291 + t294) * t67 + t299 * t103 + t300 * t102 + (t156 + t295) * t66) * m(5) + (-t112 * t125 - t113 * t267 + (t112 * t191 + t297) * t149 + (t113 * t191 - t298) * t147) * m(4); (-m(5) - m(6)) * (g(1) * t183 - g(2) * t184) + 0.2e1 * (t16 * t284 + t17 * t283) * m(6) + 0.2e1 * (t283 * t35 + t284 * t34) * m(5); -t10 * t239 / 0.2e1 + (-t193 * t53 + t208) * t285 + (-t193 * t21 + ((t191 * t26 + t4) * t184 + (-t191 * t27 + t3) * t183) * t192) * t228 + (t114 * t24 + t115 * t25 + t170 * t52 + t171 * t22 + (t183 * t5 + t184 * t6) * t246) * t265 / 0.2e1 + (-t193 * t52 + t209) * t286 + (-t193 * t22 + ((t191 * t24 + t6) * t184 + (-t191 * t25 + t5) * t183) * t192) * t230 - t193 * (t272 + t271 + (t13 * t183 + t14 * t184) * t246 + t280) / 0.2e1 + t170 * (-t193 * t65 + (t183 * t30 + t184 * t31) * t192) / 0.2e1 + t171 * (-t42 * t193 + ((t191 * t30 + t14) * t184 + (-t191 * t31 + t13) * t183) * t192) / 0.2e1 + ((t130 * t255 + t131 * t256 + t151 * t263) * t171 + (t287 * t130 + t204 * t131 + t184 * t210) * t246) * t229 + ((-t128 * t255 + t129 * t256 + t151 * t265) * t171 + (-t128 * t287 + t129 * t204 + t183 * t210) * t246) * t231 - t171 * (-t193 * t151 * t171 + ((-t194 * t255 + t195 * t256) * t171 + ((-t194 * t287 + t195 * t204) * t192 - t217 * t193) * qJD(5)) * t192) / 0.2e1 + (t114 * t26 + t115 * t27 + t170 * t53 + t171 * t21 + (t183 * t3 + t184 * t4) * t246 + t191 * t9) * t263 / 0.2e1 + ((t16 * t83 - t17 * t85 + t36 * t51 - t37 * t50) * t193 + (t15 * t218 + t40 * (-t264 * t85 - t266 * t83 + t219) + t220 * t143 + ((t191 * t36 - t17) * t184 + (t191 * t37 + t16) * t183) * t139) * t192 - (-t100 * t36 + t101 * t37) * t171 - (t40 * (t100 * t184 - t101 * t183) + t220 * t154) * t246 - g(1) * t101 - g(2) * t100 - g(3) * t154) * m(6);];
tau = t1;
