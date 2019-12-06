% Calculate vector of inverse dynamics joint torques for
% S5RPRPR3
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-12-05 17:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRPR3_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR3_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR3_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR3_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_invdynJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR3_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR3_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR3_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:51:14
% EndTime: 2019-12-05 17:51:27
% DurationCPUTime: 7.96s
% Computational Cost: add. (11923->442), mult. (9980->588), div. (0->0), fcn. (9392->10), ass. (0->219)
t193 = sin(pkin(9));
t194 = cos(pkin(9));
t195 = sin(qJ(5));
t197 = cos(qJ(5));
t140 = -rSges(6,3) * t194 + (rSges(6,1) * t197 - rSges(6,2) * t195) * t193;
t259 = qJD(5) * t193;
t323 = t140 * t259;
t199 = qJD(1) ^ 2;
t192 = qJ(1) + pkin(8);
t188 = qJ(3) + t192;
t182 = sin(t188);
t183 = cos(t188);
t191 = qJD(1) + qJD(3);
t260 = qJD(5) * t191;
t114 = (qJDD(5) * t183 - t182 * t260) * t193;
t190 = qJDD(1) + qJDD(3);
t172 = -qJDD(5) * t194 + t190;
t173 = -qJD(5) * t194 + t191;
t196 = sin(qJ(1));
t306 = pkin(1) * t196;
t185 = t199 * t306;
t187 = cos(t192);
t198 = cos(qJ(1));
t305 = pkin(1) * t198;
t240 = pkin(2) * t187 + t305;
t186 = sin(t192);
t304 = pkin(2) * t186;
t207 = -qJDD(1) * t240 + t199 * t304 + t185;
t205 = qJDD(4) * t183 + t207;
t176 = qJD(4) * t182;
t279 = t182 * t191;
t262 = -pkin(3) * t279 + t176;
t277 = t183 * t191;
t246 = -qJ(4) * t277 - t176 - t262;
t155 = (-rSges(6,1) * t195 - rSges(6,2) * t197) * t193;
t146 = qJD(5) * t155;
t252 = t146 * t259;
t255 = t194 * t279;
t274 = t191 * t193;
t256 = t182 * t274;
t263 = pkin(4) * t255 + pkin(7) * t256;
t303 = pkin(4) * t194;
t239 = pkin(7) * t193 + t303;
t132 = t239 * t183;
t281 = qJ(4) * t182;
t148 = pkin(3) * t183 + t281;
t264 = -t148 - t132;
t273 = t194 * t195;
t126 = t182 * t273 + t183 * t197;
t272 = t194 * t197;
t129 = t182 * t195 + t183 * t272;
t87 = -qJD(5) * t129 + t126 * t191;
t127 = t182 * t272 - t183 * t195;
t128 = -t182 * t197 + t183 * t273;
t88 = -qJD(5) * t128 - t127 * t191;
t49 = rSges(6,1) * t88 + rSges(6,2) * t87 - rSges(6,3) * t256;
t232 = t129 * rSges(6,1) - t128 * rSges(6,2);
t276 = t183 * t193;
t82 = -rSges(6,3) * t276 - t232;
t17 = t183 * t252 + t114 * t140 + t172 * t82 - t173 * t49 + t264 * t190 + (t246 + t263) * t191 + t205;
t322 = t17 - g(2);
t321 = t173 * t82 + t183 * t323;
t234 = rSges(4,1) * t182 + rSges(4,2) * t183;
t124 = t234 * t191;
t241 = t304 + t306;
t221 = t241 * qJD(1);
t102 = t221 + t124;
t319 = -t232 - t281;
t119 = Icges(6,4) * t129;
t75 = -Icges(6,2) * t128 + Icges(6,6) * t276 + t119;
t118 = Icges(6,4) * t128;
t79 = -Icges(6,1) * t129 - Icges(6,5) * t276 + t118;
t300 = t126 * t75 + t127 * t79;
t230 = -t128 * t75 - t129 * t79;
t72 = Icges(6,5) * t129 - Icges(6,6) * t128 + Icges(6,3) * t276;
t31 = -t194 * t72 + (-t195 * t75 - t197 * t79) * t193;
t113 = (-qJDD(5) * t182 - t183 * t260) * t193;
t177 = qJD(4) * t183;
t271 = t198 * t199;
t203 = (-qJDD(1) * t186 - t187 * t199) * pkin(2) + (-qJDD(1) * t196 - t271) * pkin(1);
t280 = qJ(4) * t183;
t231 = -pkin(3) * t182 + t280;
t261 = -t191 * t148 + t177;
t202 = qJDD(4) * t182 + t190 * t231 + t203 + (t177 + t261) * t191;
t254 = t183 * t274;
t89 = qJD(5) * t127 + t128 * t191;
t90 = qJD(5) * t126 - t129 * t191;
t294 = t90 * rSges(6,1) + t89 * rSges(6,2);
t50 = -rSges(6,3) * t254 + t294;
t270 = -t127 * rSges(6,1) + t126 * rSges(6,2);
t278 = t182 * t193;
t80 = -rSges(6,3) * t278 + t270;
t16 = -t113 * t140 + t172 * t80 + t173 * t50 - t191 ^ 2 * t132 + (-t190 * t239 + t252) * t182 + t202;
t318 = t16 - g(3);
t158 = rSges(5,2) * t254;
t290 = rSges(5,2) * t193;
t291 = rSges(5,1) * t194;
t233 = t290 - t291;
t289 = rSges(5,3) * t183;
t206 = t182 * t233 + t289;
t224 = -rSges(5,3) * t182 - t183 * t291;
t32 = t190 * t206 + (t191 * t224 + t158) * t191 + t202;
t317 = t32 - g(3);
t157 = rSges(5,1) * t255;
t167 = rSges(5,2) * t276;
t110 = -t167 - t224;
t265 = -t148 - t110;
t33 = t265 * t190 + (t157 + (-rSges(5,2) * t278 - t289) * t191 + t246) * t191 + t205;
t316 = t33 - g(2);
t165 = rSges(4,2) * t279;
t125 = -rSges(4,1) * t277 + t165;
t315 = t125 * t191 - t190 * t234 - g(3) + t203;
t292 = rSges(4,1) * t183;
t149 = -t182 * rSges(4,2) + t292;
t314 = t124 * t191 - t149 * t190 - g(2) + t207;
t311 = t240 * qJD(1);
t216 = -t303 - pkin(3) + (-rSges(6,3) - pkin(7)) * t193;
t302 = g(2) * t183;
t266 = -t191 * t231 - t176;
t209 = t221 + t266;
t310 = -t173 * t80 - t182 * t323 + t239 * t279;
t35 = t209 + t310;
t211 = t177 - t311;
t36 = t191 * t264 + t211 + t321;
t312 = t36 * (-t49 - t262 + t263) + (t17 * t216 + (-t36 * qJ(4) - t35 * (-rSges(6,3) * t193 - pkin(3) - t239)) * t191) * t183 - t35 * (-qJ(4) * t279 + t177 + t294) - t216 * t302;
t181 = t186 * rSges(3,2);
t293 = rSges(3,1) * t187;
t237 = -t293 - t305;
t145 = t181 + t237;
t235 = rSges(3,1) * t186 + rSges(3,2) * t187;
t215 = t235 + t306;
t71 = -Icges(6,5) * t127 + Icges(6,6) * t126 - Icges(6,3) * t278;
t284 = Icges(6,4) * t127;
t74 = Icges(6,2) * t126 - Icges(6,6) * t278 - t284;
t117 = Icges(6,4) * t126;
t77 = -Icges(6,1) * t127 - Icges(6,5) * t278 + t117;
t26 = -t128 * t74 + t129 * t77 + t71 * t276;
t309 = -t191 * t132 + t321;
t212 = t182 * (Icges(6,2) * t127 + t117 + t77) - t183 * (-Icges(6,2) * t129 - t118 - t79);
t213 = t182 * (-Icges(6,1) * t126 - t284 + t74) - t183 * (Icges(6,1) * t128 + t119 + t75);
t308 = t113 / 0.2e1;
t307 = t114 / 0.2e1;
t282 = Icges(6,4) * t197;
t135 = -Icges(6,6) * t194 + (-Icges(6,2) * t195 + t282) * t193;
t283 = Icges(6,4) * t195;
t136 = -Icges(6,5) * t194 + (Icges(6,1) * t197 - t283) * t193;
t152 = (-Icges(6,5) * t195 - Icges(6,6) * t197) * t193;
t141 = qJD(5) * t152;
t153 = (-Icges(6,2) * t197 - t283) * t193;
t142 = qJD(5) * t153;
t154 = (-Icges(6,1) * t195 - t282) * t193;
t143 = qJD(5) * t154;
t41 = -t141 * t194 + (-t142 * t195 + t143 * t197 + (-t135 * t197 - t136 * t195) * qJD(5)) * t193;
t134 = -Icges(6,3) * t194 + (Icges(6,5) * t197 - Icges(6,6) * t195) * t193;
t63 = -t134 * t194 + (-t135 * t195 + t136 * t197) * t193;
t301 = t63 * t172 + t41 * t173;
t299 = -t126 * t74 + t127 * t77;
t52 = -t128 * t135 + t129 * t136 + t134 * t276;
t288 = t173 * t52;
t30 = -t194 * t71 + (-t195 * t74 + t197 * t77) * t193;
t287 = t30 * t113;
t286 = t31 * t114;
t285 = rSges(5,3) + qJ(4);
t275 = t191 * t149;
t268 = t135 - t154;
t267 = t136 + t153;
t258 = m(3) + m(4) + m(5);
t250 = -pkin(3) - t291;
t249 = -t259 / 0.2e1;
t248 = t259 / 0.2e1;
t245 = t182 * t249;
t244 = t182 * t248;
t243 = t183 * t249;
t242 = t183 * t248;
t169 = rSges(2,1) * t198 - t196 * rSges(2,2);
t236 = rSges(2,1) * t196 + rSges(2,2) * t198;
t229 = -t182 * t35 + t183 * t36;
t228 = -t182 * t49 - t183 * t50;
t227 = t182 * t82 - t183 * t80;
t226 = t182 * (Icges(6,5) * t126 + Icges(6,6) * t127) - t183 * (-Icges(6,5) * t128 - Icges(6,6) * t129);
t225 = t278 * t71 + t299;
t25 = -t278 * t72 + t300;
t218 = (t182 * t225 + t183 * t25) * t193;
t27 = t276 * t72 + t230;
t217 = (-t182 * t26 + t183 * t27) * t193;
t210 = t311 - t261;
t100 = -t182 * t285 + t183 * t250 + t167;
t99 = t285 * t183 + (-pkin(3) + t233) * t182;
t54 = t182 * t216 + t270 + t280;
t10 = qJD(5) * t217 + t288;
t44 = Icges(6,5) * t90 + Icges(6,6) * t89 - Icges(6,3) * t254;
t46 = Icges(6,4) * t90 + Icges(6,2) * t89 - Icges(6,6) * t254;
t48 = Icges(6,1) * t90 + Icges(6,4) * t89 - Icges(6,5) * t254;
t13 = -t194 * t44 + (-t195 * t46 + t197 * t48 + (-t195 * t77 - t197 * t74) * qJD(5)) * t193;
t43 = Icges(6,5) * t88 + Icges(6,6) * t87 - Icges(6,3) * t256;
t45 = Icges(6,4) * t88 + Icges(6,2) * t87 - Icges(6,6) * t256;
t47 = Icges(6,1) * t88 + Icges(6,4) * t87 - Icges(6,5) * t256;
t14 = -t194 * t43 + (-t195 * t45 + t197 * t47 + (t195 * t79 - t197 * t75) * qJD(5)) * t193;
t21 = -t128 * t142 + t129 * t143 + t135 * t87 + t136 * t88 + (-t134 * t279 + t141 * t183) * t193;
t22 = t126 * t142 - t127 * t143 + t135 * t89 + t136 * t90 + (-t134 * t277 - t141 * t182) * t193;
t51 = t126 * t135 - t127 * t136 - t134 * t278;
t42 = t51 * t173;
t9 = qJD(5) * t218 + t42;
t204 = (t42 + (t300 * t183 + (t225 + t230 - t27) * t182) * t259) * t243 + t287 / 0.2e1 + t286 / 0.2e1 + t51 * t308 + t52 * t307 + t301 + (-t288 + (-(-t26 - t300) * t182 + t225 * t183 + (-t230 - t299) * t183 - t25 * t182 + (-t72 * t182 ^ 2 + (-t182 * t71 - t183 * t72) * t183) * t193) * t259 + t10) * t244 + (Icges(5,2) * t194 ^ 2 + (Icges(5,1) * t193 + 0.2e1 * Icges(5,4) * t194) * t193 + Icges(4,3)) * t190 + (t13 + t22) * t245 + (t14 + t21 + t9) * t242;
t104 = t191 * t206;
t64 = -t104 + t209;
t65 = t191 * t265 + t211;
t201 = t65 * (t157 - t262) + ((t285 * t64 - t290 * t65) * t182 + (-t250 * t64 - t285 * t65) * t183) * t191 - t64 * (t158 + t177);
t105 = t191 * t110;
t103 = -t311 - t275;
t98 = -rSges(6,1) * t128 - rSges(6,2) * t129;
t97 = rSges(6,1) * t126 + rSges(6,2) * t127;
t39 = t227 * t259 + qJD(2);
t15 = -t113 * t82 - t114 * t80 + t228 * t259 + qJDD(2);
t6 = t126 * t45 - t127 * t47 + t75 * t89 - t79 * t90 + (-t182 * t43 - t277 * t72) * t193;
t5 = t126 * t46 - t127 * t48 + t74 * t89 + t77 * t90 + (-t182 * t44 - t277 * t71) * t193;
t4 = -t128 * t45 + t129 * t47 + t75 * t87 - t79 * t88 + (t183 * t43 - t279 * t72) * t193;
t3 = -t128 * t46 + t129 * t48 + t74 * t87 + t77 * t88 + (t183 * t44 - t279 * t71) * t193;
t1 = [t204 + (Icges(2,3) + Icges(3,3)) * qJDD(1) + ((qJDD(1) * t145 + t199 * t235 - g(2) + t185) * t145 + (t271 * pkin(1) + t215 * qJDD(1) + g(3) + (-0.2e1 * t181 + t293 - t237 + t145) * t199) * t215) * m(3) + ((qJDD(1) * t236 + g(3)) * t236 + (qJDD(1) * t169 + g(2)) * t169) * m(2) + (-(t36 + t210 - t309) * t35 + (t240 * t35 + t241 * t36) * qJD(1) + t318 * (-t241 + t54) + t312 + t322 * (t319 - t240)) * m(6) + (-(t105 + t65 + t210) * t64 + (t240 * t64 + t241 * t65) * qJD(1) + t201 + t316 * (t100 - t240) + t317 * (-t241 + t99)) * m(5) + (t314 * (-t149 - t240) + t315 * (-t234 - t241) + (t191 * t292 - t165 - t275) * t102) * m(4); m(6) * t15 + t258 * qJDD(2) + (-m(6) - t258) * g(1); t204 + (-t36 * (t266 + t310) + t35 * (t261 + t309) + t318 * t54 + t312 + t322 * t319) * m(6) + (-t65 * (-t104 + t266) + t64 * (-t105 + t261) + t201 + t317 * t99 + t316 * t100) * m(5) + (-t102 * t125 + t103 * t124 + (-t102 * t191 - t314) * t149 - (t103 * t191 + t315) * t234) * m(4); (-m(5) - m(6)) * (g(3) * t182 + t302) + m(5) * (t182 * t32 + t183 * t33) + m(6) * (t16 * t182 + t17 * t183); -t194 * (t287 + t286 + (-t13 * t182 + t14 * t183) * t259 + t301) / 0.2e1 + t172 * (-t194 * t63 + (-t182 * t30 + t183 * t31) * t193) / 0.2e1 + t173 * (-t194 * t41 + ((-t191 * t30 + t14) * t183 + (-t191 * t31 - t13) * t182) * t193) / 0.2e1 - (-t113 * t225 + t114 * t25 + t172 * t51 + t173 * t22 + (-t182 * t5 + t183 * t6) * t259) * t278 / 0.2e1 + (-t194 * t51 + t218) * t308 + (-t194 * t22 + ((t191 * t225 + t6) * t183 + (-t191 * t25 - t5) * t182) * t193) * t245 + (t113 * t26 + t114 * t27 + t172 * t52 + t173 * t21 + (-t182 * t3 + t183 * t4) * t259) * t276 / 0.2e1 + (-t194 * t52 + t217) * t307 + (-t194 * t21 + ((-t191 * t26 + t4) * t183 + (-t191 * t27 - t3) * t182) * t193) * t242 - t173 * (-t194 * t152 * t173 + ((-t195 * t267 - t197 * t268) * t173 + ((t195 * t212 + t197 * t213) * t193 + t226 * t194) * qJD(5)) * t193) / 0.2e1 + ((t126 * t267 + t127 * t268 - t152 * t278) * t173 + (-t126 * t212 - t127 * t213 + t226 * t278) * t259) * t244 + ((-t128 * t267 - t129 * t268 + t152 * t276) * t173 + (t128 * t212 + t129 * t213 - t226 * t276) * t259) * t243 - (t182 * t10 + t183 * t9) * t274 / 0.2e1 + ((-t16 * t80 - t17 * t82 + t35 * t50 + t36 * t49) * t194 + (t15 * t227 + t39 * (t277 * t82 + t279 * t80 + t228) + t229 * t146 + ((-t35 * t191 + t17) * t183 + (-t36 * t191 + t16) * t182) * t140) * t193 - (-t35 * t97 - t36 * t98) * t173 - (t39 * (-t182 * t98 - t183 * t97) + t229 * t155) * t259 - g(1) * t155 - g(2) * t97 - g(3) * t98) * m(6);];
tau = t1;
