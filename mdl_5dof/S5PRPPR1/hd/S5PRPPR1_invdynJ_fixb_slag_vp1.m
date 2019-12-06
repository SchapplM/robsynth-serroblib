% Calculate vector of inverse dynamics joint torques for
% S5PRPPR1
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
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3,theta4]';
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
% Datum: 2019-12-05 15:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRPPR1_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR1_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR1_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPPR1_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR1_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPPR1_invdynJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR1_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPPR1_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPPR1_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:21:47
% EndTime: 2019-12-05 15:22:00
% DurationCPUTime: 6.99s
% Computational Cost: add. (9215->458), mult. (9002->620), div. (0->0), fcn. (8494->8), ass. (0->217)
t190 = pkin(9) + qJ(5);
t186 = sin(t190);
t188 = cos(t190);
t193 = sin(pkin(8));
t195 = cos(pkin(8));
t115 = -rSges(6,3) * t195 + (rSges(6,1) * t188 - rSges(6,2) * t186) * t193;
t176 = -qJD(5) * t195 + qJD(2);
t191 = pkin(7) + qJ(2);
t187 = sin(t191);
t257 = qJD(5) * t193;
t236 = t187 * t257;
t189 = cos(t191);
t283 = t187 * t195;
t116 = t186 * t283 + t188 * t189;
t117 = -t189 * t186 + t188 * t283;
t277 = -t117 * rSges(6,1) + t116 * rSges(6,2);
t284 = t187 * t193;
t70 = rSges(6,3) * t284 - t277;
t330 = t115 * t236 - t176 * t70;
t61 = Icges(6,5) * t117 - Icges(6,6) * t116 + Icges(6,3) * t284;
t103 = Icges(6,4) * t117;
t64 = -Icges(6,2) * t116 + Icges(6,6) * t284 + t103;
t102 = Icges(6,4) * t116;
t68 = -Icges(6,1) * t117 - Icges(6,5) * t284 + t102;
t25 = -(t186 * t64 + t188 * t68) * t193 - t195 * t61;
t328 = -t116 * t64 - t117 * t68;
t280 = t189 * t195;
t118 = -t186 * t280 + t187 * t188;
t119 = t186 * t187 + t188 * t280;
t327 = t118 * t64 - t119 * t68;
t281 = t189 * t193;
t324 = t61 * t281;
t194 = cos(pkin(9));
t278 = t194 * t195;
t192 = sin(pkin(9));
t279 = t192 * t195;
t282 = t189 * t192;
t315 = (t187 * t279 + t189 * t194) * rSges(5,2) - (t187 * t278 - t282) * rSges(5,1);
t322 = -rSges(5,3) * t284 + t315;
t285 = t187 * t192;
t321 = -(t189 * t278 + t285) * rSges(5,1) - (t187 * t194 - t189 * t279) * rSges(5,2);
t320 = -m(2) - m(3);
t18 = t324 + t327;
t317 = t18 - t324;
t143 = pkin(3) * t280 + qJ(4) * t281;
t152 = t189 * pkin(2) + t187 * qJ(3);
t271 = t143 + t152;
t80 = rSges(5,3) * t281 - t321;
t316 = t80 + t271;
t63 = Icges(6,5) * t119 + Icges(6,6) * t118 + Icges(6,3) * t281;
t290 = Icges(6,4) * t119;
t66 = Icges(6,2) * t118 + Icges(6,6) * t281 + t290;
t104 = Icges(6,4) * t118;
t69 = Icges(6,1) * t119 + Icges(6,5) * t281 + t104;
t17 = -t116 * t66 + t117 * t69 + t63 * t284;
t101 = rSges(4,1) * t280 - rSges(4,2) * t281 + t187 * rSges(4,3);
t314 = t187 * (-Icges(6,2) * t117 - t102 - t68) + t189 * (-Icges(6,2) * t119 + t104 + t69);
t312 = -m(5) - m(6);
t254 = qJD(2) * qJD(5);
t132 = (qJDD(5) * t187 + t189 * t254) * t193;
t311 = t132 / 0.2e1;
t133 = (qJDD(5) * t189 - t187 * t254) * t193;
t310 = t133 / 0.2e1;
t309 = t187 / 0.2e1;
t308 = -t189 / 0.2e1;
t307 = pkin(3) * t195;
t306 = pkin(4) * t192;
t305 = g(1) * t187;
t142 = (-rSges(6,1) * t186 - rSges(6,2) * t188) * t193;
t127 = qJD(5) * t142;
t174 = -qJDD(5) * t195 + qJDD(2);
t196 = -pkin(6) - qJ(4);
t230 = (qJ(4) + t196) * t193;
t185 = pkin(4) * t194 + pkin(3);
t233 = (pkin(3) - t185) * t195;
t252 = qJDD(4) * t193;
t255 = qJD(2) * qJD(3);
t265 = qJDD(3) * t187 + t189 * t255;
t242 = t189 * t252 + t265;
t287 = qJ(4) * t193;
t215 = t287 + t307;
t141 = t215 * t187;
t181 = t189 * qJ(3);
t150 = pkin(2) * t187 - t181;
t272 = -t141 - t150;
t224 = pkin(4) * t282 - t185 * t283;
t77 = (t230 + t307) * t187 + t224;
t245 = t77 + t272;
t247 = qJD(2) * t306;
t258 = qJD(4) * t193;
t261 = qJD(2) * t189;
t179 = qJD(3) * t189;
t121 = qJD(2) * t152 - t179;
t237 = t187 * t258;
t291 = -t215 * t261 - t121 - t237;
t83 = qJD(2) * t118 - qJD(5) * t117;
t84 = qJD(2) * t119 - qJD(5) * t116;
t223 = -rSges(6,1) * t84 - rSges(6,2) * t83;
t260 = qJD(2) * t193;
t238 = t189 * t260;
t40 = rSges(6,3) * t238 - t223;
t7 = t127 * t236 + t132 * t115 - t174 * t70 - t176 * t40 + t245 * qJDD(2) + ((-t247 - t258) * t187 + (t233 + t230) * t261 + t291) * qJD(2) + t242;
t304 = t7 * t187;
t162 = t189 * t258;
t262 = qJD(2) * t187;
t178 = qJD(3) * t187;
t264 = qJ(3) * t261 + t178;
t243 = qJD(2) * (-pkin(2) * t262 + t264) + qJDD(2) * t152 + t187 * t255;
t201 = qJDD(2) * t143 + t187 * t252 + t243 + (-t215 * t262 + 0.2e1 * t162) * qJD(2);
t239 = t187 * t260;
t269 = t189 * t247 + t196 * t239;
t81 = qJD(2) * t116 - qJD(5) * t119;
t82 = -qJD(2) * t117 + qJD(5) * t118;
t299 = t82 * rSges(6,1) + t81 * rSges(6,2);
t39 = -rSges(6,3) * t239 + t299;
t72 = t119 * rSges(6,1) + t118 * rSges(6,2) + rSges(6,3) * t281;
t209 = pkin(4) * t285 + t185 * t280 - t196 * t281;
t78 = t209 - t143;
t8 = qJDD(2) * t78 - t133 * t115 + t174 * t72 + t176 * t39 + (-t127 * t257 - qJDD(3)) * t189 + ((t233 + t287) * t262 + t269) * qJD(2) + t201;
t303 = t8 * t189;
t288 = Icges(6,4) * t188;
t113 = -Icges(6,6) * t195 + (-Icges(6,2) * t186 + t288) * t193;
t289 = Icges(6,4) * t186;
t114 = -Icges(6,5) * t195 + (Icges(6,1) * t188 - t289) * t193;
t138 = (-Icges(6,5) * t186 - Icges(6,6) * t188) * t193;
t124 = qJD(5) * t138;
t139 = (-Icges(6,2) * t188 - t289) * t193;
t125 = qJD(5) * t139;
t140 = (-Icges(6,1) * t186 - t288) * t193;
t126 = qJD(5) * t140;
t28 = -t124 * t195 + (-t125 * t186 + t126 * t188 + (-t113 * t188 - t114 * t186) * qJD(5)) * t193;
t112 = -Icges(6,3) * t195 + (Icges(6,5) * t188 - Icges(6,6) * t186) * t193;
t50 = -t112 * t195 + (-t113 * t186 + t114 * t188) * t193;
t302 = t50 * t174 + t28 * t176;
t16 = t284 * t61 + t328;
t298 = t16 * t187;
t31 = t112 * t284 - t113 * t116 + t114 * t117;
t297 = t176 * t31;
t267 = t162 + t178;
t23 = qJD(2) * t245 + t267 + t330;
t296 = t189 * t23;
t295 = t25 * t132;
t26 = -t195 * t63 + (-t186 * t66 + t188 * t69) * t193;
t294 = t26 * t133;
t293 = -rSges(5,3) - qJ(4);
t292 = -rSges(6,3) + t196;
t286 = t185 * t195;
t96 = t101 + t152;
t276 = t315 * qJD(2);
t275 = -t113 + t140;
t274 = t114 + t139;
t251 = rSges(4,1) * t283;
t266 = rSges(4,2) * t284 + t189 * rSges(4,3);
t100 = t251 - t266;
t270 = -t150 - t100;
t268 = rSges(4,2) * t239 + rSges(4,3) * t261;
t263 = -qJD(2) * t150 + t178;
t259 = qJD(4) * t187;
t256 = -m(4) + t312;
t253 = qJDD(3) * t189;
t19 = t118 * t66 + t119 * t69 + t63 * t281;
t244 = t272 + t322;
t241 = t162 + t264;
t240 = -pkin(2) - t307;
t234 = -rSges(4,1) * t195 - pkin(2);
t232 = -t257 / 0.2e1;
t231 = t257 / 0.2e1;
t229 = -qJD(2) * t141 + t162 + t263;
t175 = -qJDD(4) * t195 + qJDD(1);
t228 = t187 * t232;
t227 = t187 * t231;
t226 = t189 * t232;
t225 = t189 * t231;
t33 = Icges(6,5) * t82 + Icges(6,6) * t81 - Icges(6,3) * t239;
t34 = Icges(6,5) * t84 + Icges(6,6) * t83 + Icges(6,3) * t238;
t35 = Icges(6,4) * t82 + Icges(6,2) * t81 - Icges(6,6) * t239;
t36 = Icges(6,4) * t84 + Icges(6,2) * t83 + Icges(6,6) * t238;
t37 = Icges(6,1) * t82 + Icges(6,4) * t81 - Icges(6,5) * t239;
t38 = Icges(6,1) * t84 + Icges(6,4) * t83 + Icges(6,5) * t238;
t222 = (t118 * t36 + t119 * t38 + t64 * t81 - t68 * t82 + (t189 * t34 - t262 * t61) * t193) * t187 + t189 * (t118 * t35 + t119 * t37 + t66 * t81 + t69 * t82 + (t189 * t33 - t262 * t63) * t193);
t221 = t187 * (-t116 * t36 + t117 * t38 + t64 * t83 - t68 * t84 + (t187 * t34 + t261 * t61) * t193) + t189 * (-t116 * t35 + t117 * t37 + t66 * t83 + t69 * t84 + (t187 * t33 + t261 * t63) * t193);
t220 = t292 * t193 - pkin(2);
t10 = -t195 * t33 + (-t186 * t35 + t188 * t37 + (-t186 * t69 - t188 * t66) * qJD(5)) * t193;
t9 = -t195 * t34 + (-t186 * t36 + t188 * t38 + (t186 * t68 - t188 * t64) * qJD(5)) * t193;
t219 = t10 * t189 + t187 * t9;
t153 = rSges(3,1) * t189 - rSges(3,2) * t187;
t151 = rSges(3,1) * t187 + rSges(3,2) * t189;
t218 = t321 * qJD(2);
t24 = t176 * t72 - t179 + (-qJD(5) * t115 * t189 + t259) * t193 + (t78 + t271) * qJD(2);
t214 = t187 * t23 - t189 * t24;
t213 = -t187 * t39 + t189 * t40;
t212 = -t187 * t72 + t189 * t70;
t211 = t187 * (-Icges(6,5) * t116 - Icges(6,6) * t117) + t189 * (Icges(6,5) * t118 - Icges(6,6) * t119);
t204 = t193 * t211;
t203 = (t17 * t189 + t298) * t193;
t202 = (t18 * t187 + t189 * t19) * t193;
t200 = t293 * t193 + t240;
t199 = (Icges(6,1) * t118 - t290 - t66) * t189 + (-Icges(6,1) * t116 - t103 - t64) * t187;
t94 = qJD(2) * t96 - t179;
t93 = qJD(2) * t270 + t178;
t92 = rSges(6,1) * t118 - rSges(6,2) * t119;
t91 = -rSges(6,1) * t116 - rSges(6,2) * t117;
t47 = qJD(2) * t316 - t179 + t237;
t46 = qJD(2) * t244 + t267;
t45 = -t253 + qJDD(2) * t101 + qJD(2) * (-qJD(2) * t251 + t268) + t243;
t44 = t270 * qJDD(2) + (-qJD(2) * t101 - t121) * qJD(2) + t265;
t32 = t112 * t281 + t113 * t118 + t114 * t119;
t30 = t32 * t176;
t29 = -qJD(4) * t195 + t212 * t257 + qJD(1);
t15 = -t253 + qJDD(2) * t80 + qJD(2) * (-rSges(5,3) * t239 + t276) + t201;
t14 = t244 * qJDD(2) + ((-rSges(5,3) * t261 - t259) * t193 + t218 + t291) * qJD(2) + t242;
t13 = t113 * t83 + t114 * t84 - t116 * t125 + t117 * t126 + (t112 * t261 + t124 * t187) * t193;
t12 = t113 * t81 + t114 * t82 + t118 * t125 + t119 * t126 + (-t112 * t262 + t124 * t189) * t193;
t11 = -t132 * t72 + t133 * t70 + t213 * t257 + t175;
t6 = qJD(5) * t202 + t30;
t5 = qJD(5) * t203 + t297;
t1 = [m(5) * t175 + m(6) * t11 + (m(4) - t320) * qJDD(1) + (t256 + t320) * g(3); (t30 + ((t16 + t19 - t328) * t189 + t317 * t187) * t257) * t228 + t294 / 0.2e1 + t295 / 0.2e1 + t32 * t310 + t31 * t311 + t302 - m(3) * (-g(1) * t151 + g(2) * t153) + (t10 + t12) * t225 + (t13 + t9 + t6) * t227 + (-t297 + ((-t17 + t317 - t327) * t189 - t298) * t257 + t5) * t226 + (-t220 * t305 + t23 * (t179 + t223) + t24 * (t241 + t269 + t299) + (-t7 * pkin(2) + (-t23 * qJD(4) + t292 * t7) * t193) * t187 + ((t220 - t286) * t296 + (t23 * (-qJ(3) - t306) + t24 * (-rSges(6,3) * t193 - pkin(2) - t286)) * t187) * qJD(2) - (qJD(2) * t77 + t229 - t23 + t330) * t24 + (-g(2) + t8) * (t209 + t72 + t152) + (-g(1) + t7) * (t181 + t224 + t277)) * m(6) + (-t200 * t305 + t46 * (t179 + t218) + t47 * (t241 + t276) + (t14 * t240 + (-t46 * qJD(4) + t14 * t293) * t193) * t187 + (t46 * t200 * t189 + (-t46 * qJ(3) + t47 * (-rSges(5,3) * t193 - pkin(2) - t215)) * t187) * qJD(2) - (qJD(2) * t322 + t229 - t46) * t47 + (-g(2) + t15) * t316 + (-g(1) + t14) * (t181 + t315)) * m(5) + (t93 * t179 + t94 * (t264 + t268) + (t93 * (rSges(4,2) * t193 + t234) * t189 + (t93 * (-rSges(4,3) - qJ(3)) + t94 * t234) * t187) * qJD(2) - (-qJD(2) * t100 + t263 - t93) * t94 + (t45 - g(2)) * t96 + (t44 - g(1)) * (t234 * t187 + t181 + t266)) * m(4) + (m(3) * (t151 ^ 2 + t153 ^ 2) + Icges(3,3) + (Icges(4,2) + Icges(5,3)) * t195 ^ 2 + ((Icges(4,1) + Icges(5,1) * t194 ^ 2 + (-0.2e1 * Icges(5,4) * t194 + Icges(5,2) * t192) * t192) * t193 + 0.2e1 * (-Icges(5,5) * t194 + Icges(5,6) * t192 + Icges(4,4)) * t195) * t193) * qJDD(2); t256 * (-g(2) * t189 + t305) + 0.2e1 * (t304 / 0.2e1 - t303 / 0.2e1) * m(6) + 0.2e1 * (t14 * t309 + t15 * t308) * m(5) + 0.2e1 * (t308 * t45 + t309 * t44) * m(4); t312 * (-g(3) * t195 + (g(1) * t189 + g(2) * t187) * t193) + m(5) * (t14 * t281 + t15 * t284 - t175 * t195) + m(6) * (-t11 * t195 + t281 * t7 + t284 * t8); -t6 * t239 / 0.2e1 + (-t195 * t32 + t202) * t310 + (-t12 * t195 + ((t18 * t189 - t19 * t187) * qJD(2) + t222) * t193) * t225 + (t13 * t176 + t132 * t16 + t133 * t17 + t174 * t31 + t221 * t257) * t284 / 0.2e1 + (-t195 * t31 + t203) * t311 + (-t13 * t195 + ((t16 * t189 - t17 * t187) * qJD(2) + t221) * t193) * t227 - t195 * (t219 * t257 + t294 + t295 + t302) / 0.2e1 + t174 * (-t195 * t50 + (t187 * t25 + t189 * t26) * t193) / 0.2e1 + t176 * (-t195 * t28 + ((-t187 * t26 + t25 * t189) * qJD(2) + t219) * t193) / 0.2e1 + ((t118 * t274 + t119 * t275 + t138 * t281) * t176 + (t118 * t314 + t199 * t119 + t189 * t204) * t257) * t226 + ((-t116 * t274 + t117 * t275 + t138 * t284) * t176 + (-t116 * t314 + t117 * t199 + t187 * t204) * t257) * t228 - t176 * (-t195 * t138 * t176 + ((-t186 * t274 + t188 * t275) * t176 + ((-t186 * t314 + t188 * t199) * t193 - t211 * t195) * qJD(5)) * t193) / 0.2e1 + (qJD(2) * t5 + t12 * t176 + t132 * t18 + t133 * t19 + t174 * t32 + t222 * t257) * t281 / 0.2e1 + ((t23 * t40 - t24 * t39 + t7 * t70 - t72 * t8) * t195 + (t11 * t212 + t29 * (-t261 * t72 - t262 * t70 + t213) + t214 * t127 + (t304 - t303 + (t187 * t24 + t296) * qJD(2)) * t115) * t193 - (-t23 * t91 + t24 * t92) * t176 - (t29 * (-t187 * t92 + t189 * t91) + t214 * t142) * t257 - g(1) * t92 - g(2) * t91 - g(3) * t142) * m(6);];
tau = t1;
