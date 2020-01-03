% Calculate vector of inverse dynamics joint torques for
% S5RPRPR2
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
% Datum: 2020-01-03 11:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRPR2_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR2_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR2_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR2_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR2_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR2_invdynJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR2_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR2_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR2_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:33:38
% EndTime: 2020-01-03 11:33:51
% DurationCPUTime: 4.95s
% Computational Cost: add. (9387->400), mult. (5723->492), div. (0->0), fcn. (4240->10), ass. (0->230)
t192 = qJ(1) + pkin(8);
t182 = sin(t192);
t184 = cos(t192);
t133 = rSges(3,1) * t184 - t182 * rSges(3,2);
t197 = cos(qJ(1));
t188 = t197 * pkin(1);
t339 = t133 + t188;
t191 = qJD(1) + qJD(3);
t185 = qJ(3) + t192;
t173 = sin(t185);
t174 = cos(t185);
t330 = t173 * rSges(4,1) + t174 * rSges(4,2);
t100 = t330 * t191;
t196 = sin(qJ(1));
t187 = t196 * pkin(1);
t328 = pkin(2) * t182 + t187;
t324 = t328 * qJD(1);
t74 = t324 + t100;
t190 = pkin(9) + qJ(5);
t181 = sin(t190);
t183 = cos(t190);
t167 = Icges(6,4) * t183;
t225 = -Icges(6,2) * t181 + t167;
t329 = Icges(6,1) * t181 + t167;
t275 = t329 + t225;
t294 = Icges(6,4) * t181;
t124 = Icges(6,2) * t183 + t294;
t127 = Icges(6,1) * t183 - t294;
t276 = t124 - t127;
t338 = (t181 * t275 + t183 * t276) * t191;
t259 = qJD(5) * t191;
t105 = -qJDD(5) * t173 - t174 * t259;
t306 = rSges(6,1) * t183;
t132 = -rSges(6,2) * t181 + t306;
t113 = t132 * qJD(5);
t130 = rSges(6,1) * t181 + rSges(6,2) * t183;
t281 = t174 * t191;
t145 = qJ(4) * t281;
t166 = t173 * pkin(3);
t189 = qJDD(1) + qJDD(3);
t157 = qJD(4) * t173;
t289 = pkin(1) * qJDD(1);
t177 = t197 * t289;
t198 = qJD(1) ^ 2;
t288 = pkin(2) * qJDD(1);
t206 = t184 * t288 - t198 * t328 + t177;
t201 = -qJDD(4) * t174 + t191 * t157 + t206;
t313 = pkin(3) * t174;
t116 = qJ(4) * t173 + t313;
t194 = cos(pkin(9));
t175 = pkin(4) * t194 + pkin(3);
t138 = t174 * t175;
t195 = -pkin(7) - qJ(4);
t71 = t313 - t138 + (qJ(4) + t195) * t173;
t282 = t174 * t183;
t254 = rSges(6,1) * t282;
t283 = t174 * t181;
t240 = -rSges(6,2) * t283 + t254;
t84 = rSges(6,3) * t173 + t240;
t256 = t116 - t71 + t84;
t263 = qJD(5) * t173;
t273 = t173 * t175 + t174 * t195;
t260 = qJD(5) * t183;
t250 = rSges(6,2) * t260;
t261 = qJD(5) * t181;
t280 = t181 * t191;
t252 = rSges(6,2) * t280;
t278 = rSges(6,3) * t281 + t173 * t252;
t284 = t173 * t191;
t55 = t174 * t250 + (t174 * t261 + t183 * t284) * rSges(6,1) - t278;
t271 = t145 + t157;
t86 = pkin(3) * t284 - t271;
t8 = -t113 * t263 + t105 * t130 + t256 * t189 + (-t86 - t145 - t55 + (t166 - t273) * t191) * t191 + t201;
t337 = -g(2) + t8;
t106 = -qJDD(5) * t174 + t173 * t259;
t118 = t175 * t281;
t114 = -qJ(4) * t174 + t166;
t172 = pkin(2) * t184;
t267 = t198 * t188 + t196 * t289;
t241 = t198 * t172 + t182 * t288 + t267;
t264 = qJD(4) * t174;
t272 = pkin(3) * t281 + qJ(4) * t284;
t205 = -qJDD(4) * t173 + t189 * t114 + t241 + t191 * (-t264 + t272);
t262 = qJD(5) * t174;
t279 = t191 * t195;
t285 = t173 * t183;
t286 = t173 * t181;
t83 = rSges(6,1) * t285 - rSges(6,2) * t286 - rSges(6,3) * t174;
t61 = t83 + t273;
t312 = -t114 + t61;
t251 = rSges(6,1) * t261;
t277 = rSges(6,3) * t284 + t191 * t254;
t56 = -t173 * t251 + (-t173 * t260 - t174 * t280) * rSges(6,2) + t277;
t7 = t113 * t262 - t106 * t130 + t312 * t189 + (-t173 * t279 + t118 - t264 - t272 + t56) * t191 + t205;
t336 = -g(3) + t7;
t307 = rSges(5,1) * t194;
t152 = t173 * t307;
t193 = sin(pkin(9));
t305 = rSges(5,2) * t193;
t253 = t191 * t305;
t274 = rSges(5,3) * t281 + t173 * t253;
t153 = t174 * t305;
t255 = t174 * t307;
t88 = rSges(5,3) * t173 - t153 + t255;
t296 = t116 + t88;
t24 = (-t152 * t191 + t274 - t86) * t191 + t296 * t189 + t201;
t335 = -g(2) + t24;
t117 = rSges(4,1) * t174 - t173 * rSges(4,2);
t334 = -t100 * t191 + t117 * t189 - g(2) + t206;
t207 = t191 * t255 + rSges(5,3) * t284 + (-qJD(4) - t253) * t174;
t239 = -t173 * t305 + t152;
t87 = -rSges(5,3) * t174 + t239;
t23 = t189 * t87 + t191 * t207 + t205;
t333 = -g(3) + t23;
t101 = rSges(4,1) * t281 - rSges(4,2) * t284;
t332 = t101 * t191 + t189 * t330 - g(3) + t241;
t331 = t191 * (t114 + t87);
t235 = -t182 * rSges(3,1) - t184 * rSges(3,2);
t107 = t191 * t116;
t220 = -t130 * t263 - t264;
t72 = t191 * t84;
t326 = t191 * t71 - t107 + t118 - t220 + t277 - t72;
t323 = -t191 * t88 - t107 + t207 + t272;
t78 = -Icges(6,6) * t174 + t173 * t225;
t80 = -Icges(6,5) * t174 + t127 * t173;
t40 = t181 * t80 + t183 * t78;
t218 = t225 * t191;
t320 = -Icges(6,6) * t191 + qJD(5) * t124;
t52 = -t320 * t173 + t174 * t218;
t219 = t127 * t191;
t317 = -Icges(6,5) * t191 + qJD(5) * t329;
t54 = -t317 * t173 + t174 * t219;
t123 = Icges(6,5) * t183 - Icges(6,6) * t181;
t76 = -Icges(6,3) * t174 + t123 * t173;
t322 = qJD(5) * t40 + t181 * t52 - t183 * t54 - t191 * t76;
t122 = Icges(6,5) * t181 + Icges(6,6) * t183;
t321 = -Icges(6,3) * t191 + qJD(5) * t122;
t111 = t225 * qJD(5);
t112 = t127 * qJD(5);
t224 = t124 * t183 + t181 * t329;
t319 = qJD(5) * t224 + t111 * t181 - t112 * t183 - t122 * t191;
t79 = Icges(6,4) * t282 - Icges(6,2) * t283 + Icges(6,6) * t173;
t136 = Icges(6,4) * t283;
t81 = Icges(6,1) * t282 + Icges(6,5) * t173 - t136;
t228 = t181 * t81 + t183 * t79;
t51 = t173 * t218 + t174 * t320;
t53 = t173 * t219 + t174 * t317;
t77 = Icges(6,5) * t282 - Icges(6,6) * t283 + Icges(6,3) * t173;
t318 = qJD(5) * t228 - t181 * t51 + t183 * t53 - t191 * t77;
t315 = t105 / 0.2e1;
t314 = t106 / 0.2e1;
t311 = t173 * t329 + t78;
t310 = t174 * t329 + t79;
t309 = -t124 * t173 + t80;
t308 = -Icges(6,2) * t282 - t136 + t81;
t304 = pkin(1) * qJD(1);
t303 = t181 * t78;
t302 = t181 * t79;
t301 = t183 * t80;
t265 = qJD(1) * t184;
t270 = pkin(2) * t265 + t197 * t304;
t28 = t191 * t256 + t220 + t270;
t300 = t191 * t28;
t287 = t122 * t173;
t48 = -t124 * t283 + t282 * t329 + t287;
t299 = t48 * t191;
t298 = rSges(5,3) + qJ(4);
t92 = t122 * t174;
t217 = t123 * t191;
t268 = t172 + t188;
t266 = qJD(1) * t182;
t258 = m(3) + m(4) + m(5);
t257 = t114 + t312;
t248 = t157 + t278;
t247 = t130 * t262;
t246 = pkin(3) + t307;
t245 = -t263 / 0.2e1;
t244 = t263 / 0.2e1;
t243 = -t262 / 0.2e1;
t242 = t262 / 0.2e1;
t75 = t117 * t191 + t270;
t236 = -t264 + t270;
t156 = rSges(2,1) * t197 - t196 * rSges(2,2);
t155 = rSges(2,1) * t196 + rSges(2,2) * t197;
t209 = t324 - t157;
t27 = t191 * t257 + t209 + t247;
t233 = -t173 * t28 + t174 * t27;
t63 = t80 * t285;
t29 = -t174 * t76 - t286 * t78 + t63;
t64 = t81 * t285;
t30 = t174 * t77 + t286 * t79 - t64;
t232 = -t30 * t173 - t29 * t174;
t65 = t78 * t283;
t31 = -t173 * t76 - t282 * t80 + t65;
t227 = -t183 * t81 + t302;
t32 = t173 * t77 - t227 * t174;
t231 = -t32 * t173 - t31 * t174;
t230 = t173 * t83 + t174 * t84;
t229 = t301 - t303;
t223 = -t124 * t181 + t183 * t329;
t214 = t181 * t309 + t183 * t311;
t213 = t181 * t308 + t183 * t310;
t211 = -t173 * t217 - t174 * t321 + t191 * t227;
t210 = t321 * t173 - t174 * t217 + t191 * t229;
t208 = -t123 * qJD(5) + t191 * t223;
t62 = t138 + (rSges(6,3) - t195) * t173 + t240;
t66 = -t174 * t298 + t166 + t239;
t67 = t173 * t298 + t174 * t246 - t153;
t203 = -t246 * t284 + t271 + t274;
t47 = t173 * t223 - t92;
t42 = t47 * t191;
t11 = qJD(5) * t232 + t42;
t12 = qJD(5) * t231 - t299;
t16 = qJD(5) * t229 + t181 * t54 + t183 * t52;
t17 = qJD(5) * t227 + t181 * t53 + t183 * t51;
t20 = t208 * t173 + t174 * t319;
t21 = -t319 * t173 + t208 * t174;
t200 = (t42 + ((t31 + t64 - t65 + (t76 - t302) * t173) * t173 + (-t63 - t32 + (-t227 + t76) * t174 + (t301 + t303) * t173) * t174) * qJD(5)) * t244 - t228 * t315 - t105 * t48 / 0.2e1 + (qJD(5) * t223 + t111 * t183 + t112 * t181) * t191 + (t40 + t47) * t314 + (t12 + t299 + ((-t30 + t65 + (t77 - t301) * t174) * t174 + (t29 - t63 + (t77 + t303) * t173) * t173) * qJD(5)) * t242 + (t17 + t20 + t11) * t245 + (t16 + t21) * t243 + (Icges(4,3) + t224 + Icges(5,2) * t194 ^ 2 + (Icges(5,1) * t193 + 0.2e1 * Icges(5,4) * t194) * t193) * t189;
t199 = (t28 * (-t250 - t251 - t279) + t27 * (-qJD(4) - t252)) * t174 + (-t27 * t130 * qJD(5) + (t28 * (-t175 - t306) - t27 * t195) * t191) * t173;
t98 = t130 * t174;
t97 = t130 * t173;
t44 = t191 * t296 + t236;
t43 = t209 + t331;
t39 = qJD(5) * t230 + qJD(2);
t13 = -t105 * t83 - t106 * t84 + qJDD(2) + (t173 * t56 - t174 * t55) * qJD(5);
t6 = t318 * t173 + t211 * t174;
t5 = -t322 * t173 + t210 * t174;
t4 = t211 * t173 - t174 * t318;
t3 = t210 * t173 + t174 * t322;
t1 = [t200 + (Icges(2,3) + Icges(3,3)) * qJDD(1) + (t334 * (t117 + t268) + t332 * (t328 + t330) + (-t75 + t101 + t270) * t74) * m(4) + ((qJDD(1) * t133 - g(2) + t177) * t339 + (t198 * t339 + qJDD(1) * t235 + g(3) - t267 + (-0.2e1 * rSges(3,1) * t265 + 0.2e1 * rSges(3,2) * t266 + qJD(1) * t133) * qJD(1)) * (t235 - t187)) * m(3) + ((t155 ^ 2 + t156 ^ 2) * qJDD(1) - g(2) * t156 - g(3) * t155) * m(2) + (t28 * (-pkin(2) * t266 - t196 * t304 + t248) + t199 + t337 * (t62 + t268) + t336 * (t61 + t328) + (t28 + t326) * t27) * m(6) + (t44 * (-t324 + t203) + t335 * (t67 + t268) + t333 * (t66 + t328) + (t270 - t236 + t44 + t323) * t43) * m(5); m(6) * t13 + t258 * qJDD(2) + (-m(6) - t258) * g(1); t200 + (t257 * t300 + t199 + t337 * t62 + t336 * t61 + (-t157 + t247 + t248) * t28 + t326 * t27) * m(6) + (t335 * t67 + t333 * t66 + (-t157 + t203 + t331) * t44 + (t264 + t323) * t43) * m(5) + (-t100 * t75 + t101 * t74 + (-t191 * t74 + t334) * t117 + (t191 * t75 + t332) * t330) * m(4); (-m(5) - m(6)) * (-g(2) * t174 - g(3) * t173) + m(5) * (-t173 * t23 - t174 * t24) + m(6) * (-t173 * t7 - t174 * t8); t189 * (t173 * t228 - t174 * t40) / 0.2e1 + t191 * ((t191 * t228 - t16) * t174 + (t191 * t40 - t17) * t173) / 0.2e1 + t11 * t284 / 0.2e1 - t174 * (t30 * t105 + t29 * t106 + t47 * t189 + t21 * t191 + (-t173 * t6 - t174 * t5) * qJD(5)) / 0.2e1 + t232 * t314 + ((-t191 * t30 - t5) * t174 + (t191 * t29 - t6) * t173) * t243 - t12 * t281 / 0.2e1 - t173 * (t105 * t32 + t106 * t31 - t189 * t48 + t191 * t20 + (-t173 * t4 - t3 * t174) * qJD(5)) / 0.2e1 + t231 * t315 + ((-t191 * t32 - t3) * t174 + (t191 * t31 - t4) * t173) * t245 - t191 * ((-t276 * t181 + t275 * t183) * t191 + ((t173 * t308 - t174 * t309) * t183 + (-t173 * t310 + t174 * t311) * t181) * qJD(5)) / 0.2e1 + ((-t262 * t287 - t217) * t174 + (-t338 + (-t213 * t173 + (t92 + t214) * t174) * qJD(5)) * t173) * t242 + ((t92 * t263 - t217) * t173 + (t338 + (-t214 * t174 + (-t287 + t213) * t173) * qJD(5)) * t174) * t244 + (t13 * t230 + t39 * ((t191 * t83 - t55) * t174 + (t56 - t72) * t173) + t233 * t113 + ((t7 - t300) * t174 + (-t191 * t27 - t8) * t173) * t130 - (-t27 * t97 - t28 * t98) * t191 - (t39 * (-t173 * t97 - t174 * t98) + t233 * t132) * qJD(5) - g(1) * t132 + g(2) * t97 - g(3) * t98) * m(6);];
tau = t1;
