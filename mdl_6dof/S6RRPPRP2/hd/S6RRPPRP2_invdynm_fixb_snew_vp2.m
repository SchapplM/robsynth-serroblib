% Calculate vector of cutting torques with Newton-Euler for
% S6RRPPRP2
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
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
% m [3x7]
%   vector of cutting torques (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-06 09:13
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRPPRP2_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP2_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP2_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRP2_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP2_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP2_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP2_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRP2_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRP2_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 09:09:03
% EndTime: 2019-05-06 09:09:20
% DurationCPUTime: 7.28s
% Computational Cost: add. (103825->385), mult. (239402->450), div. (0->0), fcn. (157701->8), ass. (0->138)
t341 = sin(qJ(2));
t344 = cos(qJ(2));
t376 = qJD(1) * qJD(2);
t322 = t341 * qJDD(1) + t344 * t376;
t342 = sin(qJ(1));
t345 = cos(qJ(1));
t329 = -t345 * g(1) - t342 * g(2);
t347 = qJD(1) ^ 2;
t317 = -t347 * pkin(1) + qJDD(1) * pkin(7) + t329;
t387 = t341 * t317;
t391 = pkin(2) * t347;
t250 = qJDD(2) * pkin(2) - t322 * qJ(3) - t387 + (qJ(3) * t376 + t341 * t391 - g(3)) * t344;
t296 = -t341 * g(3) + t344 * t317;
t323 = t344 * qJDD(1) - t341 * t376;
t381 = qJD(1) * t341;
t325 = qJD(2) * pkin(2) - qJ(3) * t381;
t338 = t344 ^ 2;
t251 = t323 * qJ(3) - qJD(2) * t325 - t338 * t391 + t296;
t339 = sin(pkin(9));
t388 = cos(pkin(9));
t311 = (t339 * t344 + t388 * t341) * qJD(1);
t216 = -0.2e1 * qJD(3) * t311 + t388 * t250 - t339 * t251;
t380 = qJD(1) * t344;
t310 = t339 * t381 - t388 * t380;
t378 = qJD(3) * t310;
t305 = -0.2e1 * t378;
t385 = t339 * t250 + t388 * t251;
t217 = t305 + t385;
t266 = Ifges(4,4) * t311 - Ifges(4,2) * t310 + Ifges(4,6) * qJD(2);
t274 = -t310 * mrSges(5,2) - t311 * mrSges(5,3);
t289 = t339 * t322 - t388 * t323;
t290 = t388 * t322 + t339 * t323;
t300 = t310 * mrSges(5,1) - qJD(2) * mrSges(5,3);
t272 = t310 * pkin(3) - t311 * qJ(4);
t346 = qJD(2) ^ 2;
t397 = -2 * qJD(4);
t394 = t346 * pkin(3) - qJDD(2) * qJ(4) + qJD(2) * t397 + t310 * t272 - t385;
t210 = 0.2e1 * t378 + t394;
t301 = t311 * mrSges(5,1) + qJD(2) * mrSges(5,2);
t302 = t311 * pkin(4) - qJD(2) * pkin(8);
t309 = t310 ^ 2;
t207 = -t289 * pkin(4) - t309 * pkin(8) + qJD(2) * t302 + t305 - t394;
t340 = sin(qJ(5));
t343 = cos(qJ(5));
t293 = t343 * qJD(2) + t340 * t310;
t242 = -t293 * qJD(5) - t340 * qJDD(2) + t343 * t289;
t292 = -t340 * qJD(2) + t343 * t310;
t243 = t292 * qJD(5) + t343 * qJDD(2) + t340 * t289;
t308 = qJD(5) + t311;
t258 = -t308 * mrSges(7,2) + t292 * mrSges(7,3);
t259 = -t308 * mrSges(6,2) + t292 * mrSges(6,3);
t262 = t308 * mrSges(6,1) - t293 * mrSges(6,3);
t260 = t308 * pkin(5) - t293 * qJ(6);
t291 = t292 ^ 2;
t202 = -t242 * pkin(5) - t291 * qJ(6) + t293 * t260 + qJDD(6) + t207;
t261 = t308 * mrSges(7,1) - t293 * mrSges(7,3);
t372 = m(7) * t202 + t243 * mrSges(7,2) + t293 * t261;
t392 = -m(6) * t207 - t243 * mrSges(6,2) + (mrSges(6,1) + mrSges(7,1)) * t242 - t293 * t262 + (t258 + t259) * t292 - t372;
t353 = -m(5) * t210 + qJDD(2) * mrSges(5,3) + qJD(2) * t301 - t392;
t212 = -qJDD(2) * pkin(3) - t346 * qJ(4) + t311 * t272 + qJDD(4) - t216;
t379 = qJD(2) * t310;
t205 = (t310 * t311 - qJDD(2)) * pkin(8) + (t290 + t379) * pkin(4) + t212;
t328 = t342 * g(1) - t345 * g(2);
t365 = -qJDD(1) * pkin(1) - t328;
t257 = -t323 * pkin(2) + qJDD(3) + t325 * t381 + (-qJ(3) * t338 - pkin(7)) * t347 + t365;
t350 = (-t290 + t379) * qJ(4) + t257 + (pkin(3) * qJD(2) + t397) * t311;
t209 = -t309 * pkin(4) - t311 * t302 + (pkin(3) + pkin(8)) * t289 + t350;
t199 = t340 * t205 + t343 * t209;
t223 = Ifges(7,5) * t293 + Ifges(7,6) * t292 + Ifges(7,3) * t308;
t224 = Ifges(6,5) * t293 + Ifges(6,6) * t292 + Ifges(6,3) * t308;
t228 = Ifges(6,1) * t293 + Ifges(6,4) * t292 + Ifges(6,5) * t308;
t288 = qJDD(5) + t290;
t196 = -t291 * pkin(5) + t242 * qJ(6) + 0.2e1 * qJD(6) * t292 - t308 * t260 + t199;
t227 = Ifges(7,1) * t293 + Ifges(7,4) * t292 + Ifges(7,5) * t308;
t362 = -mrSges(7,1) * t202 + mrSges(7,3) * t196 + Ifges(7,4) * t243 + Ifges(7,2) * t242 + Ifges(7,6) * t288 + t308 * t227;
t252 = -t292 * mrSges(7,1) + t293 * mrSges(7,2);
t373 = m(7) * t196 + t242 * mrSges(7,3) + t292 * t252;
t167 = Ifges(6,4) * t243 + Ifges(6,2) * t242 + Ifges(6,6) * t288 + t308 * t228 - mrSges(6,1) * t207 + mrSges(6,3) * t199 - pkin(5) * (-t242 * mrSges(7,1) - t292 * t258 + t372) + qJ(6) * (-t288 * mrSges(7,2) - t308 * t261 + t373) + (-t224 - t223) * t293 + t362;
t198 = t343 * t205 - t340 * t209;
t193 = -0.2e1 * qJD(6) * t293 + (t292 * t308 - t243) * qJ(6) + (t292 * t293 + t288) * pkin(5) + t198;
t374 = m(7) * t193 + t288 * mrSges(7,1) + t308 * t258;
t190 = -t243 * mrSges(7,3) - t293 * t252 + t374;
t225 = Ifges(7,4) * t293 + Ifges(7,2) * t292 + Ifges(7,6) * t308;
t226 = Ifges(6,4) * t293 + Ifges(6,2) * t292 + Ifges(6,6) * t308;
t360 = mrSges(7,2) * t202 - mrSges(7,3) * t193 + Ifges(7,1) * t243 + Ifges(7,4) * t242 + Ifges(7,5) * t288 + t292 * t223;
t176 = mrSges(6,2) * t207 - mrSges(6,3) * t198 + Ifges(6,1) * t243 + Ifges(6,4) * t242 + Ifges(6,5) * t288 - qJ(6) * t190 + t292 * t224 + (-t225 - t226) * t308 + t360;
t253 = -t292 * mrSges(6,1) + t293 * mrSges(6,2);
t183 = m(6) * t198 + t288 * mrSges(6,1) + t308 * t259 + (-t252 - t253) * t293 + (-mrSges(6,3) - mrSges(7,3)) * t243 + t374;
t185 = m(6) * t199 + t242 * mrSges(6,3) + t292 * t253 + (-t261 - t262) * t308 + (-mrSges(6,2) - mrSges(7,2)) * t288 + t373;
t177 = t343 * t183 + t340 * t185;
t263 = Ifges(5,5) * qJD(2) - Ifges(5,6) * t311 + Ifges(5,3) * t310;
t357 = -mrSges(5,2) * t212 + mrSges(5,3) * t210 - Ifges(5,1) * qJDD(2) + Ifges(5,4) * t290 - Ifges(5,5) * t289 + pkin(8) * t177 + t340 * t167 - t343 * t176 + t311 * t263;
t358 = -m(5) * t212 - t290 * mrSges(5,1) - t311 * t274 - t177;
t265 = Ifges(5,4) * qJD(2) - Ifges(5,2) * t311 + Ifges(5,6) * t310;
t382 = Ifges(4,1) * t311 - Ifges(4,4) * t310 + Ifges(4,5) * qJD(2) - t265;
t399 = -mrSges(4,2) * t217 + pkin(3) * (-qJDD(2) * mrSges(5,2) - qJD(2) * t300 + t358) + qJ(4) * (-t289 * mrSges(5,1) - t310 * t274 + t353) + mrSges(4,1) * t216 + t311 * t266 - Ifges(4,6) * t289 + Ifges(4,5) * t290 + Ifges(4,3) * qJDD(2) - t357 + t382 * t310;
t361 = -mrSges(7,1) * t193 + mrSges(7,2) * t196 - Ifges(7,5) * t243 - Ifges(7,6) * t242 - Ifges(7,3) * t288 - t293 * t225;
t398 = mrSges(6,1) * t198 - mrSges(6,2) * t199 + Ifges(6,5) * t243 + Ifges(6,6) * t242 + Ifges(6,3) * t288 + pkin(5) * t190 + t293 * t226 - t361 - (t228 + t227) * t292;
t214 = t289 * pkin(3) + t350;
t396 = mrSges(5,1) * t212 - mrSges(5,3) * t214 + pkin(4) * t177 + t398;
t273 = t310 * mrSges(4,1) + t311 * mrSges(4,2);
t298 = -qJD(2) * mrSges(4,2) - t310 * mrSges(4,3);
t171 = m(4) * t216 - t290 * mrSges(4,3) - t311 * t273 + (mrSges(4,1) - mrSges(5,2)) * qJDD(2) + (t298 - t300) * qJD(2) + t358;
t299 = qJD(2) * mrSges(4,1) - t311 * mrSges(4,3);
t181 = t353 - qJDD(2) * mrSges(4,2) - qJD(2) * t299 + m(4) * t217 + (-t273 - t274) * t310 + (-mrSges(4,3) - mrSges(5,1)) * t289;
t166 = t388 * t171 + t339 * t181;
t295 = -t344 * g(3) - t387;
t313 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t341 + Ifges(3,2) * t344) * qJD(1);
t314 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t341 + Ifges(3,4) * t344) * qJD(1);
t395 = mrSges(3,1) * t295 - mrSges(3,2) * t296 + Ifges(3,5) * t322 + Ifges(3,6) * t323 + Ifges(3,3) * qJDD(2) + pkin(2) * t166 + (t341 * t313 - t344 * t314) * qJD(1) + t399;
t389 = -Ifges(5,6) - Ifges(4,4);
t178 = -t340 * t183 + t343 * t185;
t267 = Ifges(5,1) * qJD(2) - Ifges(5,4) * t311 + Ifges(5,5) * t310;
t383 = -Ifges(4,5) * t311 + Ifges(4,6) * t310 - Ifges(4,3) * qJD(2) - t267;
t321 = (-mrSges(3,1) * t344 + mrSges(3,2) * t341) * qJD(1);
t327 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t380;
t164 = m(3) * t295 + qJDD(2) * mrSges(3,1) - t322 * mrSges(3,3) + qJD(2) * t327 - t321 * t381 + t166;
t326 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t381;
t368 = -t339 * t171 + t388 * t181;
t165 = m(3) * t296 - qJDD(2) * mrSges(3,2) + t323 * mrSges(3,3) - qJD(2) * t326 + t321 * t380 + t368;
t369 = -t341 * t164 + t344 * t165;
t172 = m(5) * t214 - t289 * mrSges(5,2) - t290 * mrSges(5,3) - t310 * t300 - t311 * t301 + t178;
t356 = -mrSges(5,1) * t210 + mrSges(5,2) * t214 - pkin(4) * t392 - pkin(8) * t178 - t343 * t167 - t340 * t176;
t158 = -mrSges(4,1) * t257 + mrSges(4,3) * t217 - pkin(3) * t172 + t383 * t311 - t389 * t290 + (-Ifges(4,2) - Ifges(5,3)) * t289 + (Ifges(4,6) - Ifges(5,5)) * qJDD(2) + t382 * qJD(2) + t356;
t159 = t389 * t289 + (Ifges(4,5) - Ifges(5,4)) * qJDD(2) + (-t266 + t263) * qJD(2) + t383 * t310 + mrSges(4,2) * t257 - mrSges(4,3) * t216 - qJ(4) * t172 + (Ifges(4,1) + Ifges(5,2)) * t290 + t396;
t312 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t341 + Ifges(3,6) * t344) * qJD(1);
t316 = -t347 * pkin(7) + t365;
t354 = m(4) * t257 + t289 * mrSges(4,1) + t290 * mrSges(4,2) + t310 * t298 + t311 * t299 + t172;
t154 = -mrSges(3,1) * t316 + mrSges(3,3) * t296 + Ifges(3,4) * t322 + Ifges(3,2) * t323 + Ifges(3,6) * qJDD(2) - pkin(2) * t354 + qJ(3) * t368 + qJD(2) * t314 + t388 * t158 + t339 * t159 - t312 * t381;
t156 = mrSges(3,2) * t316 - mrSges(3,3) * t295 + Ifges(3,1) * t322 + Ifges(3,4) * t323 + Ifges(3,5) * qJDD(2) - qJ(3) * t166 - qJD(2) * t313 - t339 * t158 + t388 * t159 + t312 * t380;
t349 = -m(3) * t316 + t323 * mrSges(3,1) - t322 * mrSges(3,2) - t326 * t381 + t327 * t380 - t354;
t359 = mrSges(2,1) * t328 - mrSges(2,2) * t329 + Ifges(2,3) * qJDD(1) + pkin(1) * t349 + pkin(7) * t369 + t344 * t154 + t341 * t156;
t168 = m(2) * t328 + qJDD(1) * mrSges(2,1) - t347 * mrSges(2,2) + t349;
t162 = t344 * t164 + t341 * t165;
t160 = m(2) * t329 - t347 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t369;
t157 = mrSges(2,1) * g(3) + mrSges(2,3) * t329 + t347 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t162 - t395;
t152 = -mrSges(2,2) * g(3) - mrSges(2,3) * t328 + Ifges(2,5) * qJDD(1) - t347 * Ifges(2,6) - pkin(7) * t162 - t341 * t154 + t344 * t156;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t345 * t152 - t342 * t157 - pkin(6) * (t342 * t160 + t345 * t168), t152, t156, t159, -t310 * t265 - t357, t176, -t308 * t225 + t360; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t342 * t152 + t345 * t157 + pkin(6) * (t345 * t160 - t342 * t168), t157, t154, t158, Ifges(5,4) * qJDD(2) - Ifges(5,2) * t290 + Ifges(5,6) * t289 - qJD(2) * t263 + t310 * t267 - t396, t167, -t293 * t223 + t362; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t359, t359, t395, t399, Ifges(5,5) * qJDD(2) - Ifges(5,6) * t290 + Ifges(5,3) * t289 + qJD(2) * t265 + t311 * t267 - t356, t398, -t292 * t227 - t361;];
m_new  = t1;
