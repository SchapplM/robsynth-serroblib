% Calculate vector of cutting torques with Newton-Euler for
% S6PRRPRR4
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
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
% Datum: 2019-05-05 05:12
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6PRRPRR4_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR4_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR4_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRR4_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR4_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR4_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR4_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR4_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR4_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 05:06:34
% EndTime: 2019-05-05 05:06:51
% DurationCPUTime: 9.91s
% Computational Cost: add. (178255->341), mult. (353047->429), div. (0->0), fcn. (227069->12), ass. (0->141)
t326 = sin(pkin(11));
t328 = cos(pkin(11));
t301 = g(1) * t326 - g(2) * t328;
t302 = -g(1) * t328 - g(2) * t326;
t324 = -g(3) + qJDD(1);
t338 = cos(qJ(2));
t329 = cos(pkin(6));
t334 = sin(qJ(2));
t366 = t329 * t334;
t327 = sin(pkin(6));
t368 = t327 * t334;
t252 = t301 * t366 + t338 * t302 + t324 * t368;
t340 = qJD(2) ^ 2;
t247 = -pkin(2) * t340 + qJDD(2) * pkin(8) + t252;
t333 = sin(qJ(3));
t240 = t333 * t247;
t262 = -t301 * t327 + t324 * t329;
t337 = cos(qJ(3));
t364 = t337 * t262;
t230 = -t240 + t364;
t231 = t337 * t247 + t333 * t262;
t271 = Ifges(5,6) * qJD(3) + (Ifges(5,5) * t333 - Ifges(5,3) * t337) * qJD(2);
t274 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t333 + Ifges(4,2) * t337) * qJD(2);
t294 = (-mrSges(5,1) * t337 - mrSges(5,3) * t333) * qJD(2);
t359 = qJD(2) * qJD(3);
t358 = t337 * t359;
t296 = qJDD(2) * t333 + t358;
t357 = t333 * t359;
t297 = qJDD(2) * t337 - t357;
t293 = (-pkin(3) * t337 - qJ(4) * t333) * qJD(2);
t339 = qJD(3) ^ 2;
t360 = qJD(2) * t337;
t372 = 2 * qJD(4);
t219 = -pkin(3) * t339 + qJDD(3) * qJ(4) + qJD(3) * t372 + t293 * t360 + t231;
t361 = qJD(2) * t333;
t309 = -qJD(3) * pkin(4) - pkin(9) * t361;
t369 = t337 ^ 2 * t340;
t213 = -pkin(4) * t369 - pkin(9) * t297 + qJD(3) * t309 + t219;
t353 = -t339 * qJ(4) + t293 * t361 + qJDD(4) + t240;
t214 = -t296 * pkin(9) + (-pkin(3) - pkin(4)) * qJDD(3) + (-pkin(4) * t333 * t340 + pkin(9) * t359 - t262) * t337 + t353;
t332 = sin(qJ(5));
t336 = cos(qJ(5));
t209 = t336 * t213 + t332 * t214;
t278 = (-t332 * t337 + t333 * t336) * qJD(2);
t242 = -qJD(5) * t278 - t296 * t332 - t297 * t336;
t277 = (t332 * t333 + t336 * t337) * qJD(2);
t255 = mrSges(6,1) * t277 + mrSges(6,2) * t278;
t317 = -qJD(3) + qJD(5);
t261 = mrSges(6,1) * t317 - mrSges(6,3) * t278;
t316 = -qJDD(3) + qJDD(5);
t256 = pkin(5) * t277 - pkin(10) * t278;
t315 = t317 ^ 2;
t205 = -pkin(5) * t315 + pkin(10) * t316 - t256 * t277 + t209;
t365 = t329 * t338;
t367 = t327 * t338;
t251 = t301 * t365 - t334 * t302 + t324 * t367;
t246 = -qJDD(2) * pkin(2) - t340 * pkin(8) - t251;
t351 = -t297 * pkin(3) + t246 + (-t296 - t358) * qJ(4);
t216 = -pkin(3) * t357 + t297 * pkin(4) - pkin(9) * t369 - t351 + (t309 + t372) * t361;
t243 = -qJD(5) * t277 + t296 * t336 - t297 * t332;
t210 = (t277 * t317 - t243) * pkin(10) + (t278 * t317 - t242) * pkin(5) + t216;
t331 = sin(qJ(6));
t335 = cos(qJ(6));
t202 = -t205 * t331 + t210 * t335;
t257 = -t278 * t331 + t317 * t335;
t225 = qJD(6) * t257 + t243 * t335 + t316 * t331;
t258 = t278 * t335 + t317 * t331;
t232 = -mrSges(7,1) * t257 + mrSges(7,2) * t258;
t239 = qJDD(6) - t242;
t268 = qJD(6) + t277;
t244 = -mrSges(7,2) * t268 + mrSges(7,3) * t257;
t198 = m(7) * t202 + mrSges(7,1) * t239 - mrSges(7,3) * t225 - t232 * t258 + t244 * t268;
t203 = t205 * t335 + t210 * t331;
t224 = -qJD(6) * t258 - t243 * t331 + t316 * t335;
t245 = mrSges(7,1) * t268 - mrSges(7,3) * t258;
t199 = m(7) * t203 - mrSges(7,2) * t239 + mrSges(7,3) * t224 + t232 * t257 - t245 * t268;
t355 = -t198 * t331 + t335 * t199;
t185 = m(6) * t209 - mrSges(6,2) * t316 + mrSges(6,3) * t242 - t255 * t277 - t261 * t317 + t355;
t208 = -t213 * t332 + t214 * t336;
t260 = -mrSges(6,2) * t317 - mrSges(6,3) * t277;
t204 = -pkin(5) * t316 - pkin(10) * t315 + t256 * t278 - t208;
t349 = -m(7) * t204 + t224 * mrSges(7,1) - mrSges(7,2) * t225 + t257 * t244 - t245 * t258;
t194 = m(6) * t208 + mrSges(6,1) * t316 - mrSges(6,3) * t243 - t255 * t278 + t260 * t317 + t349;
t179 = t332 * t185 + t336 * t194;
t221 = -qJDD(3) * pkin(3) + t353 - t364;
t226 = Ifges(7,5) * t258 + Ifges(7,6) * t257 + Ifges(7,3) * t268;
t228 = Ifges(7,1) * t258 + Ifges(7,4) * t257 + Ifges(7,5) * t268;
t192 = -mrSges(7,1) * t204 + mrSges(7,3) * t203 + Ifges(7,4) * t225 + Ifges(7,2) * t224 + Ifges(7,6) * t239 - t226 * t258 + t228 * t268;
t227 = Ifges(7,4) * t258 + Ifges(7,2) * t257 + Ifges(7,6) * t268;
t193 = mrSges(7,2) * t204 - mrSges(7,3) * t202 + Ifges(7,1) * t225 + Ifges(7,4) * t224 + Ifges(7,5) * t239 + t226 * t257 - t227 * t268;
t249 = Ifges(6,4) * t278 - Ifges(6,2) * t277 + Ifges(6,6) * t317;
t250 = Ifges(6,1) * t278 - Ifges(6,4) * t277 + Ifges(6,5) * t317;
t347 = mrSges(6,1) * t208 - mrSges(6,2) * t209 + Ifges(6,5) * t243 + Ifges(6,6) * t242 + Ifges(6,3) * t316 + pkin(5) * t349 + pkin(10) * t355 + t335 * t192 + t331 * t193 + t278 * t249 + t277 * t250;
t343 = -mrSges(5,1) * t221 + mrSges(5,3) * t219 + Ifges(5,4) * t296 + Ifges(5,2) * qJDD(3) - Ifges(5,6) * t297 - pkin(4) * t179 - t347;
t306 = mrSges(5,2) * t360 + qJD(3) * mrSges(5,3);
t348 = -m(5) * t221 + qJDD(3) * mrSges(5,1) + qJD(3) * t306 - t179;
t180 = t336 * t185 - t332 * t194;
t304 = -qJD(3) * mrSges(5,1) + mrSges(5,2) * t361;
t352 = m(5) * t219 + qJDD(3) * mrSges(5,3) + qJD(3) * t304 + t294 * t360 + t180;
t275 = Ifges(5,4) * qJD(3) + (Ifges(5,1) * t333 - Ifges(5,5) * t337) * qJD(2);
t362 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t333 + Ifges(4,4) * t337) * qJD(2) + t275;
t374 = -((t271 - t274) * t333 + t337 * t362) * qJD(2) + mrSges(4,1) * t230 - mrSges(4,2) * t231 + Ifges(4,5) * t296 + Ifges(4,6) * t297 + Ifges(4,3) * qJDD(3) + pkin(3) * (-t296 * mrSges(5,2) - t294 * t361 + t348) + qJ(4) * (mrSges(5,2) * t297 + t352) + t343;
t295 = (-mrSges(4,1) * t337 + mrSges(4,2) * t333) * qJD(2);
t303 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t361;
t370 = mrSges(4,3) + mrSges(5,2);
t175 = m(4) * t231 - qJDD(3) * mrSges(4,2) - qJD(3) * t303 + t295 * t360 + t297 * t370 + t352;
t305 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t360;
t176 = m(4) * t230 + qJDD(3) * mrSges(4,1) + qJD(3) * t305 - t370 * t296 + (-t294 - t295) * t361 + t348;
t356 = t337 * t175 - t176 * t333;
t166 = m(3) * t252 - mrSges(3,1) * t340 - qJDD(2) * mrSges(3,2) + t356;
t188 = t335 * t198 + t331 * t199;
t186 = -m(6) * t216 + t242 * mrSges(6,1) - t243 * mrSges(6,2) - t277 * t260 - t278 * t261 - t188;
t222 = (pkin(3) * qJD(3) - (2 * qJD(4))) * t361 + t351;
t183 = m(5) * t222 - mrSges(5,1) * t297 - t296 * mrSges(5,3) - t304 * t361 - t306 * t360 + t186;
t342 = -m(4) * t246 + t297 * mrSges(4,1) - mrSges(4,2) * t296 - t303 * t361 + t305 * t360 - t183;
t182 = m(3) * t251 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t340 + t342;
t163 = t338 * t166 - t182 * t334;
t371 = pkin(7) * t163;
t168 = t333 * t175 + t337 * t176;
t167 = m(3) * t262 + t168;
t159 = t166 * t366 - t167 * t327 + t182 * t365;
t272 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t333 + Ifges(4,6) * t337) * qJD(2);
t273 = Ifges(5,2) * qJD(3) + (Ifges(5,4) * t333 - Ifges(5,6) * t337) * qJD(2);
t248 = Ifges(6,5) * t278 - Ifges(6,6) * t277 + Ifges(6,3) * t317;
t170 = mrSges(6,2) * t216 - mrSges(6,3) * t208 + Ifges(6,1) * t243 + Ifges(6,4) * t242 + Ifges(6,5) * t316 - pkin(10) * t188 - t192 * t331 + t193 * t335 - t248 * t277 - t249 * t317;
t344 = mrSges(7,1) * t202 - mrSges(7,2) * t203 + Ifges(7,5) * t225 + Ifges(7,6) * t224 + Ifges(7,3) * t239 + t227 * t258 - t228 * t257;
t171 = -mrSges(6,1) * t216 + mrSges(6,3) * t209 + Ifges(6,4) * t243 + Ifges(6,2) * t242 + Ifges(6,6) * t316 - pkin(5) * t188 - t248 * t278 + t250 * t317 - t344;
t345 = -mrSges(5,1) * t222 + mrSges(5,2) * t219 - pkin(4) * t186 - pkin(9) * t180 - t332 * t170 - t336 * t171;
t155 = -mrSges(4,1) * t246 + mrSges(4,3) * t231 - pkin(3) * t183 + (Ifges(4,2) + Ifges(5,3)) * t297 + (Ifges(4,4) - Ifges(5,5)) * t296 + (Ifges(4,6) - Ifges(5,6)) * qJDD(3) + t362 * qJD(3) + (-t272 - t273) * t361 + t345;
t346 = mrSges(5,2) * t221 - mrSges(5,3) * t222 + Ifges(5,1) * t296 + Ifges(5,4) * qJDD(3) - Ifges(5,5) * t297 - pkin(9) * t179 + qJD(3) * t271 + t336 * t170 - t171 * t332 + t273 * t360;
t160 = mrSges(4,2) * t246 - mrSges(4,3) * t230 + Ifges(4,1) * t296 + Ifges(4,4) * t297 + Ifges(4,5) * qJDD(3) - qJ(4) * t183 - qJD(3) * t274 + t272 * t360 + t346;
t150 = mrSges(3,1) * t251 - mrSges(3,2) * t252 + Ifges(3,3) * qJDD(2) + pkin(2) * t342 + pkin(8) * t356 + t337 * t155 + t333 * t160;
t152 = mrSges(3,2) * t262 - mrSges(3,3) * t251 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t340 - pkin(8) * t168 - t155 * t333 + t160 * t337;
t154 = -mrSges(3,1) * t262 + mrSges(3,3) * t252 + t340 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t168 - t374;
t350 = mrSges(2,1) * t301 - mrSges(2,2) * t302 + pkin(1) * t159 + t329 * t150 + t152 * t368 + t154 * t367 + t327 * t371;
t161 = m(2) * t302 + t163;
t158 = t329 * t167 + (t166 * t334 + t182 * t338) * t327;
t156 = m(2) * t301 + t159;
t148 = mrSges(2,2) * t324 - mrSges(2,3) * t301 + t338 * t152 - t334 * t154 + (-t158 * t327 - t159 * t329) * pkin(7);
t147 = -mrSges(2,1) * t324 + mrSges(2,3) * t302 - pkin(1) * t158 - t327 * t150 + (t152 * t334 + t154 * t338 + t371) * t329;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t328 * t148 - t326 * t147 - qJ(1) * (t156 * t328 + t161 * t326), t148, t152, t160, t346, t170, t193; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t326 * t148 + t328 * t147 + qJ(1) * (-t156 * t326 + t161 * t328), t147, t154, t155, t343 + (-t333 * t271 - t337 * t275) * qJD(2), t171, t192; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t350, t350, t150, t374, Ifges(5,5) * t296 + Ifges(5,6) * qJDD(3) - Ifges(5,3) * t297 - qJD(3) * t275 + t273 * t361 - t345, t347, t344;];
m_new  = t1;
