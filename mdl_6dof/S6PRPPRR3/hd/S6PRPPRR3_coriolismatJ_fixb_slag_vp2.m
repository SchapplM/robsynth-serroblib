% Calculate matrix of centrifugal and coriolis load on the joints for
% S6PRPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta4]';
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
% Cq [6x6]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6PRPPRR3_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR3_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR3_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR3_coriolismatJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPPRR3_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPPRR3_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPPRR3_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:21:34
% EndTime: 2019-03-08 19:21:40
% DurationCPUTime: 3.38s
% Computational Cost: add. (5709->361), mult. (13494->541), div. (0->0), fcn. (13545->10), ass. (0->226)
t227 = sin(qJ(6));
t360 = t227 / 0.2e1;
t230 = cos(qJ(6));
t357 = -t230 / 0.2e1;
t338 = t230 * mrSges(7,2);
t340 = t227 * mrSges(7,1);
t190 = t338 + t340;
t231 = cos(qJ(5));
t170 = t231 * t190;
t228 = sin(qJ(5));
t310 = t227 * t228;
t295 = mrSges(7,3) * t310;
t183 = -mrSges(7,2) * t231 + t295;
t307 = t228 * t230;
t294 = mrSges(7,3) * t307;
t185 = mrSges(7,1) * t231 + t294;
t383 = t183 * t357 + t185 * t360;
t243 = t170 / 0.2e1 - t383;
t189 = -mrSges(7,1) * t230 + mrSges(7,2) * t227;
t374 = m(7) * pkin(5);
t385 = -mrSges(6,1) + t189 - t374;
t218 = Ifges(7,4) * t230;
t192 = Ifges(7,1) * t227 + t218;
t219 = t227 ^ 2;
t221 = t230 ^ 2;
t300 = t219 + t221;
t337 = t231 * mrSges(6,2);
t339 = t228 * mrSges(6,1);
t382 = -t337 / 0.2e1 - t339 / 0.2e1;
t223 = sin(pkin(11));
t225 = cos(pkin(11));
t232 = -pkin(2) - pkin(3);
t301 = t225 * qJ(3) + t223 * t232;
t182 = -pkin(8) + t301;
t354 = pkin(5) * t228;
t197 = pkin(9) * t231 - t354;
t103 = t182 * t310 + t197 * t230;
t104 = -t182 * t307 + t197 * t227;
t381 = -t103 * t227 + t104 * t230;
t274 = Ifges(7,2) * t227 - t218;
t304 = t230 * t231;
t167 = t223 * t227 + t225 * t304;
t322 = t167 * t230;
t309 = t227 * t231;
t165 = t223 * t230 - t225 * t309;
t324 = t165 * t227;
t245 = (-t324 / 0.2e1 + t322 / 0.2e1) * mrSges(7,3);
t279 = -t223 * qJ(3) + t225 * t232;
t181 = pkin(4) - t279;
t352 = t231 * pkin(5);
t353 = t228 * pkin(9);
t134 = t181 + t352 + t353;
t91 = t134 * t230 - t182 * t309;
t92 = t134 * t227 + t182 * t304;
t272 = t227 * t91 - t230 * t92;
t318 = t182 * t231;
t248 = t272 + 0.2e1 * t318;
t266 = t322 - t324;
t253 = t266 * pkin(9);
t164 = -t223 * t309 - t225 * t230;
t166 = t223 * t304 - t227 * t225;
t186 = -mrSges(7,1) * t228 + mrSges(7,3) * t304;
t363 = -t186 / 0.2e1;
t184 = mrSges(7,2) * t228 + mrSges(7,3) * t309;
t365 = -t184 / 0.2e1;
t256 = t164 * t363 + t166 * t365;
t269 = t103 * t164 + t104 * t166;
t169 = t228 * t190;
t368 = t169 / 0.2e1;
t375 = m(7) / 0.2e1;
t376 = -m(7) / 0.2e1;
t12 = (-t225 * mrSges(6,2) + t223 * t368) * t231 + t245 + t269 * t376 + t253 * t375 + ((t189 / 0.2e1 - mrSges(6,1) - t374 / 0.2e1) * t225 + (t248 * t376 + t243) * t223) * t228 + t256;
t224 = sin(pkin(6));
t355 = cos(qJ(2));
t288 = t224 * t355;
t229 = sin(qJ(2));
t314 = t224 * t229;
t148 = -t223 * t288 + t225 * t314;
t332 = cos(pkin(6));
t112 = t148 * t228 + t231 * t332;
t113 = t148 * t231 - t228 * t332;
t147 = (t223 * t229 + t225 * t355) * t224;
t70 = t113 * t230 + t147 * t227;
t335 = t70 * t230;
t69 = -t113 * t227 + t147 * t230;
t336 = t69 * t227;
t263 = t113 - t335 + t336;
t323 = t166 * t230;
t325 = t164 * t227;
t267 = -t323 + t325;
t20 = (t267 * t112 + (t112 * t231 + t228 * t263) * t223) * t375;
t220 = t228 ^ 2;
t222 = t231 ^ 2;
t298 = -0.1e1 + t300;
t54 = -t267 * t231 + (-t220 * t298 - t222) * t223;
t371 = t54 / 0.2e1;
t315 = t223 * t231;
t316 = t223 * t228;
t68 = (t267 + t315) * t316;
t380 = (t68 * qJD(3) + qJD(4) * t371) * m(7) + t20 * qJD(1) - t12 * qJD(2);
t379 = -m(5) * t301 + t228 * t169 - mrSges(5,2);
t277 = t231 * mrSges(6,1) - t228 * mrSges(6,2);
t378 = m(5) * t279 - m(6) * t181 - mrSges(5,1) - t277;
t377 = m(6) / 0.2e1;
t373 = -mrSges(7,1) / 0.2e1;
t372 = mrSges(7,2) / 0.2e1;
t370 = t70 / 0.2e1;
t168 = t228 * t189;
t369 = -t168 / 0.2e1;
t285 = t168 / 0.2e1;
t367 = -t169 / 0.2e1;
t366 = -t183 / 0.2e1;
t364 = -t185 / 0.2e1;
t362 = t190 / 0.2e1;
t361 = -t227 / 0.2e1;
t358 = -t228 / 0.2e1;
t356 = -t231 / 0.2e1;
t351 = t20 * qJD(5);
t348 = m(7) * qJD(5);
t350 = t348 * t371;
t349 = m(7) * qJD(2);
t347 = mrSges(7,3) * t228;
t345 = Ifges(7,4) * t227;
t217 = Ifges(7,5) * t230;
t344 = Ifges(7,5) * t231;
t342 = Ifges(7,6) * t227;
t341 = Ifges(7,6) * t231;
t82 = -t147 * t309 - t148 * t230;
t334 = t82 * t227;
t83 = t147 * t304 - t148 * t227;
t333 = t83 * t230;
t329 = t112 * t228;
t328 = t147 * t222;
t327 = t147 * t223;
t326 = t148 * t225;
t236 = (t223 * t285 + (t323 / 0.2e1 - t325 / 0.2e1) * mrSges(7,3)) * t228 + t164 * t183 / 0.2e1 + t166 * t364;
t260 = t165 * t373 + t167 * t372;
t17 = t236 + t260;
t321 = t17 * qJD(2);
t320 = t182 * t220;
t319 = t182 * t228;
t317 = t223 * t225;
t313 = t225 * t228;
t311 = t227 * t186;
t305 = t230 * t184;
t278 = mrSges(7,3) * (t221 / 0.2e1 + t219 / 0.2e1);
t38 = t168 * t356 + (t183 * t361 + t185 * t357) * t228 + t220 * t278;
t303 = t38 * qJD(2);
t302 = Ifges(7,5) * t310 + Ifges(7,6) * t307;
t122 = t298 * t231 * t228;
t299 = t122 * qJD(4);
t254 = m(7) * t381;
t27 = (t220 - t222) * t182 * t375 + (t272 * t376 + t243) * t231 + (t254 / 0.2e1 + t305 / 0.2e1 - t311 / 0.2e1 + t367) * t228;
t292 = qJD(3) * t375;
t73 = (-t225 * t231 + t266) * t228;
t297 = t27 * qJD(5) + t38 * qJD(6) + t73 * t292;
t273 = t333 - t334;
t41 = (-t147 * t231 + t273) * t228;
t296 = m(7) * t41 * qJD(1);
t293 = t349 / 0.2e1;
t291 = qJD(4) * t376;
t287 = t316 / 0.2e1;
t280 = t300 * t231;
t276 = -t337 - t339;
t275 = Ifges(7,1) * t230 - t345;
t25 = (-t231 * t263 - t298 * t329) * t375;
t270 = -t25 * qJD(1) - t27 * qJD(2);
t268 = t113 * t231 + t329;
t264 = t82 * mrSges(7,1) / 0.2e1 - t83 * mrSges(7,2) / 0.2e1;
t262 = pkin(5) * t369 + t231 * t217 / 0.4e1;
t261 = t103 * t373 + t104 * t372;
t259 = t338 / 0.2e1 + t340 / 0.2e1;
t149 = t228 * t274 + t341;
t172 = t228 * t192;
t258 = -t149 / 0.4e1 + t172 / 0.4e1 + pkin(9) * t366;
t208 = Ifges(7,4) * t310;
t151 = -Ifges(7,1) * t307 + t208 + t344;
t171 = Ifges(7,2) * t307 + t208;
t257 = pkin(9) * t364 + t151 / 0.4e1 + t171 / 0.4e1;
t191 = Ifges(7,2) * t230 + t345;
t255 = t191 * t360 + t192 * t357;
t14 = m(6) * (-t148 + t268) * t147 + (t147 * t329 + t69 * t82 + t70 * t83) * m(7);
t252 = -t14 * qJD(1) + t291 * t41;
t251 = t217 / 0.2e1 - t342 / 0.2e1 - Ifges(6,4);
t234 = (t333 / 0.2e1 - t334 / 0.2e1) * mrSges(7,3) + (pkin(9) * t273 - t147 * t354) * t375 + t147 * t285;
t239 = t113 * t368 + t363 * t69 + t365 * t70;
t242 = t103 * t69 + t104 * t70 + t113 * t319;
t249 = t272 + t318;
t1 = t242 * t376 + (t249 * t376 + t243) * t112 + t234 + t239;
t150 = -Ifges(7,6) * t228 + t231 * t274;
t152 = -Ifges(7,5) * t228 - t231 * t275;
t5 = t103 * t185 + t91 * t186 + m(7) * (t103 * t91 + t104 * t92) + t104 * t183 + t92 * t184 + (-t181 * mrSges(6,2) + t149 * t360 + t151 * t357 - t182 * t169 - t251 * t231) * t231 + (-t181 * mrSges(6,1) + t152 * t357 + t150 * t360 - t182 * t170 + t251 * t228 + (m(7) * t182 ^ 2 + Ifges(6,1) - Ifges(6,2) - Ifges(7,3)) * t231) * t228;
t250 = -t1 * qJD(1) + t5 * qJD(2) + t27 * qJD(4);
t247 = -t190 / 0.2e1 + t259;
t11 = -t92 * t185 + t231 * t302 / 0.2e1 + t91 * t183 + (t182 * t168 + (t92 * mrSges(7,3) - t172 / 0.2e1 + t149 / 0.2e1) * t230 + (-t91 * mrSges(7,3) + t151 / 0.2e1 + t171 / 0.2e1) * t227) * t228;
t240 = t112 * t369 + t185 * t370 + t366 * t69;
t3 = (-t335 / 0.2e1 + t336 / 0.2e1) * t347 + t240 + t264;
t246 = -t3 * qJD(1) + t11 * qJD(2) + t38 * qJD(4);
t21 = m(7) * t263 * t112;
t244 = t21 * qJD(1) + t20 * qJD(3) + t25 * qJD(4);
t22 = -mrSges(4,3) - t165 * t185 - t167 * t183 - m(4) * qJ(3) - m(7) * (t165 * t91 + t167 * t92) + t378 * t223 + ((-m(6) * t182 + mrSges(6,3)) * (t220 + t222) - m(7) * t320 + t379) * t225;
t235 = -m(6) * (t225 * t268 + t327) / 0.2e1 + (t112 * t313 + t165 * t69 + t167 * t70) * t376;
t117 = t220 * t327;
t237 = (t222 * t327 + t117 + t326) * t377 + (t164 * t82 + t166 * t83 + t117) * t375;
t7 = t235 + t237;
t241 = t7 * qJD(1) + t22 * qJD(2) + t291 * t73;
t110 = t247 * t231;
t23 = t247 * t112;
t44 = t247 * t316;
t61 = pkin(5) * t190 - t274 * t357 + t275 * t361 + t255;
t233 = pkin(9) * t278 + t182 * t362 + (t192 - t274) * t227 / 0.4e1 + (-t275 / 0.4e1 + t191 / 0.4e1) * t230;
t9 = (t344 / 0.2e1 + t257) * t230 + (-0.3e1 / 0.4e1 * t341 + t258) * t227 + (Ifges(7,3) / 0.2e1 + t233) * t228 + t261 + t262;
t238 = t23 * qJD(1) - t9 * qJD(2) + t44 * qJD(3) - t110 * qJD(4) + t61 * qJD(5);
t195 = t220 * t317;
t111 = t190 * t356 - t231 * t259;
t99 = t147 * t320;
t45 = t190 * t287 + t259 * t316;
t24 = (t259 + t362) * t112;
t16 = t236 - t260;
t13 = -t170 * t287 + t315 * t367 + t245 - t256 + (-pkin(5) * t313 + t253 + t269) * t375 + (t248 * t375 + t383) * t316 + (-t276 / 0.2e1 + t285 + t382) * t225;
t10 = -Ifges(7,5) * t304 / 0.2e1 + Ifges(7,6) * t309 / 0.2e1 + Ifges(7,3) * t358 + (-t341 / 0.4e1 + t258) * t227 + t257 * t230 + t233 * t228 - t261 + t262;
t8 = t25 * qJD(5) + t293 * t41;
t6 = m(4) * t314 + m(5) * (t326 + t327) - t235 + t237;
t4 = t294 * t370 - t69 * t295 / 0.2e1 - t240 + t264;
t2 = t242 * t375 + t234 - t239 + (t276 / 0.2e1 + t382) * t147 + (t249 * t375 - t243) * t112;
t15 = [t14 * qJD(2) + t21 * qJD(5), t6 * qJD(3) + t2 * qJD(5) + t4 * qJD(6) - t252 + (t83 * t183 + t82 * t185 + m(7) * (t91 * t82 + t92 * t83 + t99) + m(6) * (t182 * t328 + t99) - mrSges(6,3) * t328 + m(4) * (-pkin(2) * t229 + qJ(3) * t355) * t224 + t378 * t148 + (-t220 * mrSges(6,3) - t379) * t147 + (-mrSges(3,1) - mrSges(4,1)) * t314 + (-mrSges(3,2) + mrSges(4,3)) * t288) * qJD(2), qJD(2) * t6 + t351, t8, t2 * qJD(2) + (t385 * t113 + (mrSges(6,2) - (m(7) * pkin(9) + mrSges(7,3)) * t300) * t112) * qJD(5) + t24 * qJD(6) + t244, t4 * qJD(2) + t24 * qJD(5) + (-mrSges(7,1) * t70 - mrSges(7,2) * t69) * qJD(6); -t7 * qJD(3) - t1 * qJD(5) - t3 * qJD(6) + t252, -qJD(3) * t22 + qJD(5) * t5 + qJD(6) * t11, 0.2e1 * ((t164 * t165 + t166 * t167 + t195) * t375 + (t195 + (t222 - 0.1e1) * t317) * t377) * qJD(3) + t13 * qJD(5) + t16 * qJD(6) - t241, -t296 / 0.2e1 + t297, t13 * qJD(3) + t10 * qJD(6) + t250 + (mrSges(6,2) * t319 + (Ifges(7,5) * t227 + Ifges(7,6) * t230) * t358 + t152 * t360 + t230 * t150 / 0.2e1 + pkin(5) * t170 + Ifges(6,6) * t228 + (t254 + t305 - t311) * pkin(9) + (t182 * t385 - Ifges(6,5) + t255) * t231 + t381 * mrSges(7,3)) * qJD(5), t16 * qJD(3) + t10 * qJD(5) + (-mrSges(7,1) * t92 - mrSges(7,2) * t91 + t302) * qJD(6) + t246; qJD(2) * t7 + t351, -t12 * qJD(5) + t17 * qJD(6) + t241, t68 * t348, -t73 * t349 / 0.2e1 + t350, t45 * qJD(6) + (t231 * t189 + m(7) * (-t300 * t353 - t352) - t300 * t347 - t277) * qJD(5) * t223 + t380, t321 + t45 * qJD(5) + (-mrSges(7,1) * t166 - mrSges(7,2) * t164) * qJD(6); t8, t296 / 0.2e1 + t297, t293 * t73 + t350, t122 * t348, t54 * t292 + m(7) * t299 + (t168 + m(7) * (pkin(9) * t280 - t354) + mrSges(7,3) * t280 + t276) * qJD(5) + t111 * qJD(6) - t270, t111 * qJD(5) + t168 * qJD(6) + t303; qJD(2) * t1 - qJD(6) * t23 - t244, qJD(3) * t12 + qJD(6) * t9 - t250, -t44 * qJD(6) - t380, t110 * qJD(6) + (-t54 * qJD(3) / 0.2e1 - t299) * m(7) + t270, -t61 * qJD(6) (pkin(9) * t189 + t217 - t342) * qJD(6) - t238; t3 * qJD(2) + t23 * qJD(5), -qJD(3) * t17 - qJD(5) * t9 - t246, t44 * qJD(5) - t321, -t110 * qJD(5) - t303, t238, 0;];
Cq  = t15;
