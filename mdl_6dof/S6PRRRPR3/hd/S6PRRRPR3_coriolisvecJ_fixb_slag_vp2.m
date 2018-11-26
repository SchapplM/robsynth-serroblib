% Calculate vector of centrifugal and coriolis load on the joints for
% S6PRRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 15:23
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6PRRRPR3_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR3_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR3_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR3_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR3_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR3_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR3_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:23:29
% EndTime: 2018-11-23 15:23:39
% DurationCPUTime: 9.79s
% Computational Cost: add. (6262->523), mult. (15695->699), div. (0->0), fcn. (11183->10), ass. (0->260)
t212 = sin(qJ(2));
t207 = sin(pkin(6));
t287 = qJD(1) * t207;
t269 = t212 * t287;
t211 = sin(qJ(3));
t283 = qJD(3) * t211;
t353 = pkin(3) * t283 - t269;
t367 = -pkin(9) - pkin(8);
t186 = t367 * t211;
t328 = cos(qJ(3));
t187 = t367 * t328;
t210 = sin(qJ(4));
t327 = cos(qJ(4));
t149 = t210 * t186 - t187 * t327;
t177 = t210 * t328 + t211 * t327;
t180 = qJD(3) * t186;
t225 = qJD(3) * t187;
t214 = cos(qJ(2));
t266 = t214 * t287;
t373 = -qJD(4) * t149 + t177 * t266 - t210 * t180 + t327 * t225;
t206 = qJD(3) + qJD(4);
t290 = t210 * t211;
t248 = t206 * t290;
t252 = t327 * t328;
t139 = -t206 * t252 + t248;
t381 = qJ(5) * t139 - qJD(5) * t177 + t353;
t380 = mrSges(6,2) - mrSges(5,1);
t379 = -t139 * pkin(5) - t373;
t262 = qJD(4) * t327;
t208 = cos(pkin(6));
t271 = t208 * t328;
t193 = qJD(1) * t271;
t182 = qJD(2) * pkin(8) + t269;
t241 = (-pkin(9) * qJD(2) - t182) * t211;
t141 = t193 + t241;
t286 = qJD(1) * t211;
t151 = t182 * t328 + t208 * t286;
t264 = qJD(2) * t328;
t223 = -pkin(9) * t264 - t151;
t291 = t210 * t223;
t73 = t141 * t327 + t291;
t378 = pkin(3) * t262 - t73;
t140 = t206 * t177;
t342 = pkin(4) + pkin(10);
t377 = -t140 * t342 - t381;
t376 = -Ifges(5,5) + Ifges(6,4);
t375 = -Ifges(5,6) + Ifges(6,5);
t238 = qJD(2) * t252;
t285 = qJD(2) * t211;
t168 = t210 * t285 - t238;
t209 = sin(qJ(6));
t213 = cos(qJ(6));
t147 = t168 * t209 + t206 * t213;
t159 = Ifges(5,4) * t168;
t169 = t177 * qJD(2);
t162 = qJD(6) + t169;
t362 = t162 * Ifges(7,3);
t146 = t168 * t213 - t206 * t209;
t363 = t146 * Ifges(7,6);
t374 = t169 * Ifges(5,1) + t206 * Ifges(5,5) + t147 * Ifges(7,5) - t159 + t362 + t363;
t359 = -qJD(5) - t378;
t125 = t140 * qJD(2);
t124 = qJD(2) * t248 - t206 * t238;
t284 = qJD(2) * t212;
t261 = qJD(1) * t284;
t188 = t207 * t261;
t279 = qJD(2) * qJD(3);
t260 = t211 * t279;
t163 = pkin(3) * t260 + t188;
t229 = qJ(5) * t124 - qJD(5) * t169 + t163;
t24 = t125 * t342 + t229;
t323 = t169 * pkin(5);
t137 = qJD(3) * pkin(3) + t141;
t70 = -t327 * t137 - t291;
t237 = t70 + t323;
t370 = qJD(5) + t237;
t35 = -t206 * t342 + t370;
t202 = -pkin(3) * t328 - pkin(2);
t160 = qJD(2) * t202 - t266;
t222 = -t169 * qJ(5) + t160;
t67 = t168 * t342 + t222;
t17 = t209 * t35 + t213 * t67;
t303 = qJD(6) * t17;
t293 = t207 * t214;
t267 = qJD(2) * t293;
t253 = t211 * t267;
t218 = -qJD(1) * t253 + qJD(3) * t223;
t136 = t327 * t223;
t71 = t210 * t137 - t136;
t272 = t207 * t328;
t256 = t214 * t272;
t239 = qJD(2) * t256;
t288 = qJD(1) * t239 + qJD(3) * t193;
t97 = qJD(3) * t241 + t288;
t21 = qJD(4) * t71 + t210 * t97 - t327 * t218;
t7 = -t124 * pkin(5) + t21;
t2 = -t209 * t24 + t213 * t7 - t303;
t16 = -t209 * t67 + t213 * t35;
t243 = t16 * t209 - t17 * t213;
t1 = qJD(6) * t16 + t209 * t7 + t213 * t24;
t325 = t1 * t209;
t221 = m(7) * (-qJD(6) * t243 + t2 * t213 + t325);
t68 = qJD(6) * t146 + t125 * t209;
t30 = -mrSges(7,1) * t124 - mrSges(7,3) * t68;
t69 = -qJD(6) * t147 + t125 * t213;
t31 = mrSges(7,2) * t124 + mrSges(7,3) * t69;
t372 = t209 * t31 + t213 * t30 + t221;
t247 = mrSges(7,1) * t213 - mrSges(7,2) * t209;
t324 = t168 * pkin(5);
t64 = -t206 * qJ(5) - t71;
t44 = -t64 - t324;
t312 = t147 * Ifges(7,4);
t61 = t146 * Ifges(7,2) + t162 * Ifges(7,6) + t312;
t371 = -t213 * t61 / 0.2e1 + t247 * t44;
t89 = t168 * pkin(4) + t222;
t369 = t160 * mrSges(5,1) - t89 * mrSges(6,2);
t368 = t160 * mrSges(5,2) - t17 * mrSges(7,2) - t89 * mrSges(6,3);
t344 = t68 / 0.2e1;
t343 = t69 / 0.2e1;
t341 = -t124 / 0.2e1;
t148 = -t327 * t186 - t187 * t210;
t110 = pkin(5) * t177 + t148;
t176 = -t252 + t290;
t230 = -t177 * qJ(5) + t202;
t90 = t176 * t342 + t230;
t36 = t110 * t213 - t209 * t90;
t366 = qJD(6) * t36 + t209 * t379 - t213 * t377;
t37 = t110 * t209 + t213 * t90;
t365 = -qJD(6) * t37 + t209 * t377 + t213 * t379;
t358 = t323 - t359;
t357 = pkin(4) * t140 + t381;
t356 = -qJD(5) - t70;
t320 = mrSges(5,3) * t168;
t152 = -mrSges(5,2) * t206 - t320;
t321 = mrSges(6,1) * t168;
t154 = -mrSges(6,3) * t206 + t321;
t355 = t152 - t154;
t319 = mrSges(5,3) * t169;
t354 = -mrSges(6,1) * t169 - t206 * t380 - t319;
t114 = -t182 * t283 + t288;
t263 = qJD(3) * t328;
t115 = -t182 * t263 + (-qJD(3) * t208 - t267) * t286;
t352 = t328 * t114 - t115 * t211;
t349 = -t124 * t148 - t125 * t149 + t177 * t21;
t63 = -pkin(4) * t206 - t356;
t348 = m(6) * t63 - t354;
t347 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t81 = -mrSges(7,1) * t146 + mrSges(7,2) * t147;
t346 = -m(5) * t71 + m(6) * t64 - m(7) * t44 - t355 - t81;
t216 = qJD(2) ^ 2;
t345 = Ifges(7,1) * t344 + Ifges(7,4) * t343 + Ifges(7,5) * t341;
t339 = -t146 / 0.2e1;
t338 = -t147 / 0.2e1;
t337 = t147 / 0.2e1;
t336 = -t162 / 0.2e1;
t335 = -t168 / 0.2e1;
t334 = t168 / 0.2e1;
t333 = -t169 / 0.2e1;
t332 = t169 / 0.2e1;
t330 = -t206 / 0.2e1;
t329 = t206 / 0.2e1;
t326 = pkin(3) * t210;
t318 = mrSges(7,3) * t213;
t317 = Ifges(4,4) * t211;
t316 = Ifges(7,4) * t209;
t315 = Ifges(7,4) * t213;
t164 = -t211 * t207 * t212 + t271;
t165 = t208 * t211 + t212 * t272;
t112 = -t164 * t327 + t165 * t210;
t314 = t112 * t21;
t313 = t124 * mrSges(6,1);
t311 = t148 * t21;
t310 = t169 * Ifges(5,4);
t309 = t169 * Ifges(6,6);
t308 = t169 * t64;
t304 = t154 - t81;
t297 = t140 * t209;
t296 = t169 * t209;
t295 = t176 * t209;
t294 = t176 * t213;
t289 = t213 * t140;
t282 = qJD(4) * t210;
t281 = qJD(6) * t209;
t280 = qJD(6) * t213;
t278 = Ifges(7,5) * t68 + Ifges(7,6) * t69 - Ifges(7,3) * t124;
t277 = t327 * pkin(3);
t204 = pkin(3) * t285;
t275 = Ifges(4,4) * t328;
t274 = mrSges(4,3) * t285;
t268 = t207 * t284;
t258 = -t281 / 0.2e1;
t72 = t141 * t210 - t136;
t255 = mrSges(4,3) * t264;
t254 = t211 * t266;
t201 = -t277 - pkin(4);
t251 = qJD(2) * t263;
t246 = Ifges(7,1) * t209 + t315;
t245 = Ifges(7,2) * t213 + t316;
t244 = Ifges(7,5) * t209 + Ifges(7,6) * t213;
t127 = pkin(4) * t169 + qJ(5) * t168;
t94 = -mrSges(7,2) * t162 + mrSges(7,3) * t146;
t95 = mrSges(7,1) * t162 - mrSges(7,3) * t147;
t242 = -t209 * t95 + t213 * t94;
t240 = qJD(1) * t256;
t100 = t127 + t204;
t236 = mrSges(4,1) * t211 + mrSges(4,2) * t328;
t235 = Ifges(4,5) * t328 - Ifges(4,6) * t211;
t234 = -t112 * t209 + t213 * t293;
t86 = t112 * t213 + t209 * t293;
t233 = t176 * t280 + t297;
t232 = t176 * t281 - t289;
t20 = t137 * t262 + t210 * t218 + t223 * t282 + t327 * t97;
t74 = -t327 * t180 - t186 * t262 - t187 * t282 - t210 * t225;
t228 = t211 * (Ifges(4,1) * t328 - t317);
t227 = qJD(3) * t236;
t226 = (Ifges(4,2) * t328 + t317) * qJD(2);
t15 = -qJD(5) * t206 - t20;
t220 = qJD(3) * t165 + t253;
t219 = qJD(3) * t164 + t239;
t10 = t68 * Ifges(7,4) + t69 * Ifges(7,2) - t124 * Ifges(7,6);
t105 = t206 * Ifges(6,5) + t168 * Ifges(6,3) - t309;
t158 = Ifges(6,6) * t168;
t106 = t206 * Ifges(6,4) - t169 * Ifges(6,2) + t158;
t107 = -t168 * Ifges(5,2) + t206 * Ifges(5,6) + t310;
t6 = -pkin(5) * t125 - t15;
t143 = Ifges(7,4) * t146;
t62 = t147 * Ifges(7,1) + t162 * Ifges(7,5) + t143;
t217 = t6 * (mrSges(7,1) * t209 + mrSges(7,2) * t213) - t209 * t10 / 0.2e1 - t20 * mrSges(5,2) - t15 * mrSges(6,3) + t70 * t320 + t63 * t321 + (Ifges(7,5) * t213 - Ifges(7,6) * t209) * t341 + (-Ifges(7,2) * t209 + t315) * t343 + (Ifges(7,1) * t213 - t316) * t344 + t213 * t345 + (-t296 / 0.2e1 + t258) * t62 + (t158 + t106) * t335 + (-t310 + t105) * t333 + (t309 + t107) * t332 + t380 * t21 + (-Ifges(5,1) * t333 - Ifges(7,5) * t338 + Ifges(6,2) * t332 - Ifges(7,6) * t339 - Ifges(7,3) * t336 + t330 * t376 + t368) * t168 + (Ifges(6,3) * t335 - t17 * t318 + t244 * t336 + t245 * t339 + t246 * t338 + t330 * t375 - t369 + t371) * t169 + (mrSges(7,1) * t168 + (t281 + t296) * mrSges(7,3)) * t16 + t375 * t125 + t376 * t124 + t371 * qJD(6) + (-Ifges(5,2) * t169 - t159 + t374) * t334 - (t146 * t245 + t147 * t246 + t162 * t244) * qJD(6) / 0.2e1;
t203 = Ifges(4,4) * t264;
t199 = qJ(5) + t326;
t185 = -qJD(3) * mrSges(4,2) + t255;
t184 = qJD(3) * mrSges(4,1) - t274;
t183 = -qJD(2) * pkin(2) - t266;
t173 = qJD(2) * t227;
t167 = Ifges(4,1) * t285 + Ifges(4,5) * qJD(3) + t203;
t166 = Ifges(4,6) * qJD(3) + t226;
t161 = t169 * pkin(10);
t150 = -t182 * t211 + t193;
t132 = t176 * pkin(4) + t230;
t130 = -mrSges(6,2) * t168 - mrSges(6,3) * t169;
t129 = mrSges(5,1) * t168 + mrSges(5,2) * t169;
t113 = t210 * t164 + t165 * t327;
t111 = -t176 * pkin(5) + t149;
t88 = t127 + t161;
t80 = t100 + t161;
t59 = mrSges(5,1) * t125 - mrSges(5,2) * t124;
t58 = -mrSges(6,2) * t125 + mrSges(6,3) * t124;
t51 = t72 - t324;
t50 = t71 - t324;
t42 = -pkin(5) * t140 - t74;
t33 = t164 * t282 + t165 * t262 + t210 * t219 + t220 * t327;
t29 = pkin(4) * t125 + t229;
t27 = -mrSges(7,1) * t69 + mrSges(7,2) * t68;
t26 = t209 * t50 + t213 * t88;
t25 = -t209 * t88 + t213 * t50;
t23 = t209 * t51 + t213 * t80;
t22 = -t209 * t80 + t213 * t51;
t14 = qJD(6) * t234 - t209 * t268 + t213 * t33;
t13 = qJD(6) * t86 + t209 * t33 + t213 * t268;
t3 = [m(4) * (-t207 ^ 2 * t214 * t261 + t114 * t165 + t115 * t164 - t150 * t220 + t151 * t219) + m(5) * (t314 + t113 * t20 + (t160 * t284 - t163 * t214) * t207) + m(6) * (t314 - t113 * t15 + (-t214 * t29 + t284 * t89) * t207) + t219 * t185 - t220 * t184 + t113 * t27 + t13 * t94 + t14 * t95 + t86 * t30 - t234 * t31 + m(7) * (-t1 * t234 + t113 * t6 + t13 * t17 + t14 * t16 + t2 * t86) + (m(5) * t70 + t348) * t33 + (-mrSges(3,1) * t212 - mrSges(3,2) * t214) * t207 * t216 + (-t164 * t251 - t165 * t260) * mrSges(4,3) + t346 * (-t164 * t262 + t165 * t282 + t210 * t220 - t219 * t327) + (-t173 - t59 - t58) * t293 + (m(4) * t183 + t129 + t130 + (-mrSges(4,1) * t328 + t211 * mrSges(4,2)) * qJD(2)) * t268 + (mrSges(5,3) + mrSges(6,1)) * (-t112 * t124 - t113 * t125); t146 * (Ifges(7,4) * t233 - Ifges(7,2) * t232) / 0.2e1 + (t124 * t176 + t140 * t333) * Ifges(6,6) + (-Ifges(4,2) * t211 + t275) * t251 + (Ifges(7,1) * t233 - Ifges(7,4) * t232) * t337 + (-t140 * t71 + t349) * mrSges(5,3) + (t140 * t64 + t349) * mrSges(6,1) + t62 * t297 / 0.2e1 + t61 * t289 / 0.2e1 + (t1 * t294 - t16 * t233 - t17 * t232 - t2 * t295) * mrSges(7,3) + t365 * t95 + (t1 * t37 + t111 * t6 + t16 * t365 + t17 * t366 + t2 * t36 + t42 * t44) * m(7) + t366 * t94 + qJD(3) ^ 2 * t235 / 0.2e1 + t357 * t130 + t44 * (mrSges(7,1) * t232 + mrSges(7,2) * t233) + t162 * (Ifges(7,5) * t233 - Ifges(7,6) * t232) / 0.2e1 + (t132 * t29 - t149 * t15 + t357 * t89 - t373 * t63 + t64 * t74 + t311) * m(6) + (t149 * t20 + t160 * t353 + t163 * t202 - t373 * t70 - t71 * t74 + t311) * m(5) + (t163 * mrSges(5,1) + t15 * mrSges(6,1) - t29 * mrSges(6,2) - t20 * mrSges(5,3) + t244 * t341 + t245 * t343 + t246 * t344 - t6 * t247 + t61 * t258 + (Ifges(6,3) + Ifges(5,2)) * t125 + (t124 / 0.2e1 - t341) * Ifges(5,4)) * t176 + (t334 * Ifges(6,6) - t70 * mrSges(5,3) - t63 * mrSges(6,1) + Ifges(6,2) * t333 - t362 / 0.2e1 - t363 / 0.2e1 - t16 * mrSges(7,1) + t106 / 0.2e1 - Ifges(5,1) * t332 - Ifges(5,4) * t335 - Ifges(7,5) * t337 + t376 * t329 - t368 - t374 / 0.2e1) * t139 - t355 * t74 + (-(t183 * t212 + (-t150 * t211 + t151 * t328) * t214) * t287 - pkin(2) * t188 + ((-t150 * t328 - t151 * t211) * qJD(3) + t352) * pkin(8)) * m(4) + (-t150 * t263 - t151 * t283 + t352) * mrSges(4,3) + t353 * t129 + (-t184 * t263 - t185 * t283) * pkin(8) + (qJD(2) * (Ifges(4,1) * t211 + t275) + t167) * t263 / 0.2e1 - (t226 + t166) * t283 / 0.2e1 + (qJD(6) * t62 + t10) * t294 / 0.2e1 + t346 * (-t210 * t254 + t240 * t327) + t183 * t227 + t202 * t59 - pkin(2) * t173 + (-t29 * mrSges(6,3) + t163 * mrSges(5,2) + Ifges(7,6) * t343 + Ifges(7,5) * t344 + (Ifges(7,3) + Ifges(5,1)) * t341 + t347 + t278 / 0.2e1 + (-Ifges(6,2) - Ifges(5,1) / 0.2e1) * t124 + (-Ifges(6,6) - Ifges(5,4)) * t125) * t177 + t132 * t58 + t111 * t27 + t42 * t81 + (Ifges(6,3) * t334 + t105 / 0.2e1 - t107 / 0.2e1 - Ifges(5,4) * t332 - Ifges(5,2) * t335 + t375 * t329 + t369) * t140 + t36 * t30 + t37 * t31 + t373 * t354 + t184 * t254 + t228 * t279 + t295 * t345 - t185 * t240; t217 - t129 * t204 + Ifges(4,5) * t251 + t166 * t285 / 0.2e1 - t235 * t279 / 0.2e1 + t358 * t81 + (-t16 * t22 - t17 * t23 + t199 * t6 + t358 * t44) * m(7) - t216 * t228 / 0.2e1 + (m(7) * (t16 * t213 + t17 * t209) + t213 * t95 + t209 * t94 + t348) * pkin(3) * t282 + t354 * t72 - t201 * t313 + t359 * t154 + (-t100 * t89 - t15 * t199 + t201 * t21 + t359 * t64 - t63 * t72) * m(6) - (-Ifges(4,2) * t285 + t167 + t203) * t264 / 0.2e1 + (t274 + t184) * t151 + (t255 - t185) * t150 + (-t17 * t280 - t325) * mrSges(7,3) - t183 * t236 * qJD(2) + (t124 * t277 - t125 * t326) * mrSges(5,3) + (-t125 * t199 - t308) * mrSges(6,1) + ((-t327 * t21 + t20 * t210 + (t210 * t70 + t327 * t71) * qJD(4)) * pkin(3) - t160 * t204 - t70 * t72 - t71 * t73) * m(5) + t199 * t27 - t100 * t130 - t114 * mrSges(4,2) + t115 * mrSges(4,1) - t23 * t94 - t22 * t95 + (t280 * t94 - t281 * t95 + t372) * (-pkin(10) + t201) + t378 * t152 + t71 * t319 - t2 * t318 - Ifges(4,6) * t260; t217 - t342 * t221 + (t319 + t354) * t71 + t355 * t70 + (pkin(4) * t124 - qJ(5) * t125 - t308) * mrSges(6,1) - t304 * qJD(5) + (-(qJD(6) * t94 + t30) * t342 + (-t2 - t303) * mrSges(7,3)) * t213 + (-t1 * mrSges(7,3) - (-qJD(6) * t95 + t31) * t342) * t209 - t127 * t130 - t26 * t94 - t25 * t95 + t237 * t81 + qJ(5) * t27 + (qJ(5) * t6 - t16 * t25 - t17 * t26 + t370 * t44) * m(7) + (-pkin(4) * t21 - qJ(5) * t15 - t127 * t89 + t356 * t64 - t63 * t71) * m(6); -t313 + t304 * t206 + t242 * qJD(6) + (t130 + t242) * t169 - m(7) * (t169 * t243 + t44 * t206) + (t169 * t89 + t206 * t64 + t21) * m(6) + t372; -t44 * (mrSges(7,1) * t147 + mrSges(7,2) * t146) + (Ifges(7,1) * t146 - t312) * t338 + t61 * t337 + (Ifges(7,5) * t146 - Ifges(7,6) * t147) * t336 - t16 * t94 + t17 * t95 + (t146 * t16 + t147 * t17) * mrSges(7,3) + t278 + (-Ifges(7,2) * t147 + t143 + t62) * t339 + t347;];
tauc  = t3(:);
