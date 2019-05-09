% Calculate vector of cutting torques with Newton-Euler for
% S6RPRPPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2]';
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
% Datum: 2019-05-05 16:43
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPRPPR3_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR3_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR3_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPPR3_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR3_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR3_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR3_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR3_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR3_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 16:41:31
% EndTime: 2019-05-05 16:41:39
% DurationCPUTime: 3.99s
% Computational Cost: add. (43145->350), mult. (85677->407), div. (0->0), fcn. (40783->8), ass. (0->135)
t337 = sin(qJ(1));
t340 = cos(qJ(1));
t308 = t337 * g(1) - g(2) * t340;
t288 = qJDD(1) * pkin(1) + t308;
t309 = -g(1) * t340 - g(2) * t337;
t342 = qJD(1) ^ 2;
t293 = -pkin(1) * t342 + t309;
t332 = sin(pkin(9));
t333 = cos(pkin(9));
t236 = t332 * t288 + t333 * t293;
t224 = -pkin(2) * t342 + qJDD(1) * pkin(7) + t236;
t330 = -g(3) + qJDD(2);
t336 = sin(qJ(3));
t339 = cos(qJ(3));
t219 = -t336 * t224 + t339 * t330;
t220 = t339 * t224 + t336 * t330;
t255 = Ifges(5,6) * qJD(3) + (Ifges(5,5) * t336 - Ifges(5,3) * t339) * qJD(1);
t259 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t336 + Ifges(4,2) * t339) * qJD(1);
t290 = (-mrSges(5,1) * t339 - mrSges(5,3) * t336) * qJD(1);
t370 = qJD(1) * qJD(3);
t366 = t339 * t370;
t295 = qJDD(1) * t336 + t366;
t365 = t336 * t370;
t296 = qJDD(1) * t339 - t365;
t371 = qJD(1) * t339;
t305 = -qJD(3) * mrSges(6,1) + mrSges(6,3) * t371;
t341 = qJD(3) ^ 2;
t289 = (-pkin(3) * t339 - qJ(4) * t336) * qJD(1);
t372 = qJD(1) * t336;
t361 = t289 * t372 + qJDD(4) - t219;
t217 = -qJDD(3) * pkin(3) - t341 * qJ(4) + t361;
t377 = t339 * t342;
t369 = qJD(1) * qJD(5);
t392 = -0.2e1 * t336 * t369 + (-t295 + t366) * qJ(5);
t211 = (-t336 * t377 - qJDD(3)) * pkin(4) + t217 + t392;
t292 = (mrSges(6,1) * t336 - mrSges(6,2) * t339) * qJD(1);
t301 = -qJD(3) * pkin(4) - qJ(5) * t372;
t235 = t333 * t288 - t332 * t293;
t223 = -qJDD(1) * pkin(2) - t342 * pkin(7) - t235;
t360 = -t296 * pkin(3) + t223 + (-t295 - t366) * qJ(4);
t378 = t339 ^ 2 * t342;
t388 = 2 * qJD(4);
t350 = -qJ(5) * t378 + qJDD(5) - t360 + (t301 + t388) * t372;
t386 = pkin(4) + pkin(8);
t387 = -pkin(3) - pkin(8);
t200 = t295 * pkin(5) + (pkin(5) * t339 + t336 * t387) * t370 + t386 * t296 + t350;
t294 = (pkin(5) * t336 + pkin(8) * t339) * qJD(1);
t203 = (-pkin(5) - qJ(4)) * t341 + (-pkin(4) * t377 - qJD(1) * t294) * t336 + (-pkin(3) - t386) * qJDD(3) + t361 + t392;
t335 = sin(qJ(6));
t338 = cos(qJ(6));
t196 = t200 * t338 - t203 * t335;
t286 = -qJD(3) * t338 + t335 * t371;
t233 = qJD(6) * t286 - qJDD(3) * t335 - t296 * t338;
t287 = -qJD(3) * t335 - t338 * t371;
t237 = -mrSges(7,1) * t286 + mrSges(7,2) * t287;
t312 = qJD(6) + t372;
t238 = -mrSges(7,2) * t312 + mrSges(7,3) * t286;
t284 = qJDD(6) + t295;
t193 = m(7) * t196 + mrSges(7,1) * t284 - mrSges(7,3) * t233 - t237 * t287 + t238 * t312;
t197 = t200 * t335 + t203 * t338;
t232 = -qJD(6) * t287 - qJDD(3) * t338 + t296 * t335;
t239 = mrSges(7,1) * t312 - mrSges(7,3) * t287;
t194 = m(7) * t197 - mrSges(7,2) * t284 + mrSges(7,3) * t232 + t237 * t286 - t239 * t312;
t376 = -t335 * t193 + t338 * t194;
t362 = -m(6) * t211 + t292 * t372 - t376;
t355 = -qJDD(3) * mrSges(6,2) - qJD(3) * t305 + t362;
t174 = -t295 * mrSges(6,3) - t355;
t317 = qJD(3) * t388;
t384 = t341 * pkin(3);
t394 = qJDD(3) * qJ(4) + t289 * t371 + t220;
t215 = t317 - t384 + t394;
t357 = pkin(4) * t378 + t296 * qJ(5) - t394;
t202 = qJDD(3) * pkin(5) + qJD(3) * t301 + t317 + t387 * t341 + (-0.2e1 * qJD(5) - t294) * t371 - t357;
t225 = Ifges(7,5) * t287 + Ifges(7,6) * t286 + Ifges(7,3) * t312;
t227 = Ifges(7,1) * t287 + Ifges(7,4) * t286 + Ifges(7,5) * t312;
t185 = -mrSges(7,1) * t202 + mrSges(7,3) * t197 + Ifges(7,4) * t233 + Ifges(7,2) * t232 + Ifges(7,6) * t284 - t225 * t287 + t227 * t312;
t226 = Ifges(7,4) * t287 + Ifges(7,2) * t286 + Ifges(7,6) * t312;
t186 = mrSges(7,2) * t202 - mrSges(7,3) * t196 + Ifges(7,1) * t233 + Ifges(7,4) * t232 + Ifges(7,5) * t284 + t225 * t286 - t226 * t312;
t389 = -2 * qJD(4);
t210 = 0.2e1 * t339 * t369 + t384 + (t389 - t301) * qJD(3) + t357;
t257 = -Ifges(6,6) * qJD(3) + (-Ifges(6,4) * t339 - Ifges(6,2) * t336) * qJD(1);
t260 = -Ifges(6,5) * qJD(3) + (-Ifges(6,1) * t339 - Ifges(6,4) * t336) * qJD(1);
t359 = -m(7) * t202 + t232 * mrSges(7,1) - t233 * mrSges(7,2) + t286 * t238 - t287 * t239;
t354 = mrSges(6,1) * t210 - mrSges(6,2) * t211 - Ifges(6,5) * t296 - Ifges(6,6) * t295 - Ifges(6,3) * qJDD(3) + pkin(5) * t359 + pkin(8) * t376 + t338 * t185 + t335 * t186 - t257 * t371 + t260 * t372;
t347 = -mrSges(5,1) * t217 + mrSges(5,3) * t215 + Ifges(5,4) * t295 + Ifges(5,2) * qJDD(3) - Ifges(5,6) * t296 - pkin(4) * t174 - t354;
t304 = -qJD(3) * mrSges(5,1) + mrSges(5,2) * t372;
t302 = qJD(3) * mrSges(6,2) - mrSges(6,3) * t372;
t352 = -m(6) * t210 + qJDD(3) * mrSges(6,1) - t296 * mrSges(6,3) + qJD(3) * t302 - t359;
t348 = m(5) * t215 + qJDD(3) * mrSges(5,3) + qJD(3) * t304 + t290 * t371 + t352;
t367 = t292 * t371;
t261 = Ifges(5,4) * qJD(3) + (Ifges(5,1) * t336 - Ifges(5,5) * t339) * qJD(1);
t373 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t336 + Ifges(4,4) * t339) * qJD(1) + t261;
t307 = mrSges(5,2) * t371 + qJD(3) * mrSges(5,3);
t393 = -m(5) * t217 + qJDD(3) * mrSges(5,1) + qJD(3) * t307;
t395 = -((t255 - t259) * t336 + t339 * t373) * qJD(1) + mrSges(4,1) * t219 - mrSges(4,2) * t220 + Ifges(4,5) * t295 + Ifges(4,6) * t296 + Ifges(4,3) * qJDD(3) + pkin(3) * (-t290 * t372 + (-mrSges(5,2) + mrSges(6,3)) * t295 + t355 + t393) + qJ(4) * (t296 * mrSges(5,2) + t348 - t367) + t347;
t254 = -Ifges(6,3) * qJD(3) + (-Ifges(6,5) * t339 - Ifges(6,6) * t336) * qJD(1);
t391 = Ifges(6,4) * t296 + Ifges(6,2) * t295 - t254 * t371;
t383 = mrSges(4,3) + mrSges(5,2);
t382 = Ifges(5,6) - Ifges(6,5);
t291 = (-mrSges(4,1) * t339 + mrSges(4,2) * t336) * qJD(1);
t306 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t371;
t170 = m(4) * t219 + (mrSges(4,1) - mrSges(6,2)) * qJDD(3) + (-t305 + t306) * qJD(3) + (-t290 - t291) * t372 + (mrSges(6,3) - t383) * t295 + t362 + t393;
t303 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t372;
t181 = -qJD(3) * t303 + m(4) * t220 - qJDD(3) * mrSges(4,2) + t383 * t296 + (t291 - t292) * t371 + t348;
t363 = -t170 * t336 + t339 * t181;
t163 = m(3) * t236 - mrSges(3,1) * t342 - qJDD(1) * mrSges(3,2) + t363;
t177 = t338 * t193 + t335 * t194;
t207 = -pkin(3) * t365 + t296 * pkin(4) + t350;
t173 = -m(6) * t207 - t295 * mrSges(6,1) + t296 * mrSges(6,2) - t302 * t372 + t305 * t371 - t177;
t212 = (pkin(3) * qJD(3) + t389) * t372 + t360;
t171 = m(5) * t212 - mrSges(5,1) * t296 - t295 * mrSges(5,3) - t304 * t372 - t307 * t371 + t173;
t345 = -m(4) * t223 + t296 * mrSges(4,1) - mrSges(4,2) * t295 - t303 * t372 + t306 * t371 - t171;
t167 = m(3) * t235 + qJDD(1) * mrSges(3,1) - mrSges(3,2) * t342 + t345;
t160 = t332 * t163 + t333 * t167;
t165 = t339 * t170 + t336 * t181;
t258 = Ifges(5,2) * qJD(3) + (Ifges(5,4) * t336 - Ifges(5,6) * t339) * qJD(1);
t375 = -t254 + t258;
t364 = t333 * t163 - t167 * t332;
t256 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t336 + Ifges(4,6) * t339) * qJD(1);
t353 = -mrSges(6,2) * t207 + mrSges(6,3) * t210 + Ifges(6,1) * t296 + Ifges(6,4) * t295 + pkin(8) * t177 - qJD(3) * t257 + t335 * t185 - t338 * t186;
t346 = mrSges(5,1) * t212 - mrSges(5,2) * t215 + pkin(4) * t173 + qJ(5) * (t352 - t367) - t353;
t154 = (Ifges(4,2) + Ifges(5,3)) * t296 + (Ifges(4,4) - Ifges(5,5)) * t295 + (Ifges(4,6) - t382) * qJDD(3) + t373 * qJD(3) - mrSges(4,1) * t223 + (-t256 - t375) * t372 + mrSges(4,3) * t220 - pkin(3) * t171 - t346;
t356 = mrSges(7,1) * t196 - mrSges(7,2) * t197 + Ifges(7,5) * t233 + Ifges(7,6) * t232 + Ifges(7,3) * t284 + t287 * t226 - t286 * t227;
t349 = mrSges(6,1) * t207 - mrSges(6,3) * t211 + Ifges(6,6) * qJDD(3) + pkin(5) * t177 + qJD(3) * t260 + t356;
t344 = mrSges(5,2) * t217 - mrSges(5,3) * t212 + Ifges(5,1) * t295 + Ifges(5,4) * qJDD(3) - Ifges(5,5) * t296 - qJ(5) * t174 + qJD(3) * t255 + t258 * t371 + t349;
t156 = -qJD(3) * t259 + mrSges(4,2) * t223 + (-t254 + t256) * t371 - mrSges(4,3) * t219 - qJ(4) * t171 + Ifges(4,5) * qJDD(3) + t344 + (Ifges(6,4) + Ifges(4,4)) * t296 + (Ifges(4,1) + Ifges(6,2)) * t295;
t358 = mrSges(3,1) * t235 - mrSges(3,2) * t236 + Ifges(3,3) * qJDD(1) + pkin(2) * t345 + pkin(7) * t363 + t339 * t154 + t336 * t156;
t351 = mrSges(2,1) * t308 - mrSges(2,2) * t309 + Ifges(2,3) * qJDD(1) + pkin(1) * t160 + t358;
t158 = m(2) * t309 - mrSges(2,1) * t342 - qJDD(1) * mrSges(2,2) + t364;
t157 = m(2) * t308 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t342 + t160;
t152 = -mrSges(3,1) * t330 + mrSges(3,3) * t236 + t342 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t165 - t395;
t151 = mrSges(3,2) * t330 - mrSges(3,3) * t235 + Ifges(3,5) * qJDD(1) - Ifges(3,6) * t342 - pkin(7) * t165 - t154 * t336 + t156 * t339;
t150 = -mrSges(2,2) * g(3) - mrSges(2,3) * t308 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t342 - qJ(2) * t160 + t151 * t333 - t152 * t332;
t149 = Ifges(2,6) * qJDD(1) + t342 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t309 + t332 * t151 + t333 * t152 - pkin(1) * (m(3) * t330 + t165) + qJ(2) * t364;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t340 * t150 - t337 * t149 - pkin(6) * (t157 * t340 + t158 * t337), t150, t151, t156, t344 + t391, -Ifges(6,5) * qJDD(3) - t254 * t372 - t353, t186; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t337 * t150 + t340 * t149 + pkin(6) * (-t157 * t337 + t158 * t340), t149, t152, t154, (-t336 * t255 - t339 * t261) * qJD(1) + t347, -t349 - t391, t185; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t351, t351, t358, t395, Ifges(5,5) * t295 - Ifges(5,3) * t296 - qJD(3) * t261 + qJDD(3) * t382 + t372 * t375 + t346, t354, t356;];
m_new  = t1;
