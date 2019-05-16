% Calculate vector of cutting torques with Newton-Euler for
% S6RPRPRP6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
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
% Datum: 2019-05-05 17:55
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPRPRP6_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP6_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP6_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP6_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP6_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP6_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP6_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP6_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP6_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 17:52:32
% EndTime: 2019-05-05 17:52:46
% DurationCPUTime: 6.94s
% Computational Cost: add. (95607->363), mult. (229379->421), div. (0->0), fcn. (160739->8), ass. (0->138)
t337 = sin(qJ(1));
t339 = cos(qJ(1));
t316 = -g(1) * t339 - g(2) * t337;
t341 = qJD(1) ^ 2;
t309 = -pkin(1) * t341 + qJDD(1) * qJ(2) + t316;
t333 = sin(pkin(9));
t334 = cos(pkin(9));
t378 = qJD(1) * qJD(2);
t370 = -t334 * g(3) - 0.2e1 * t333 * t378;
t390 = pkin(2) * t334;
t258 = (-pkin(7) * qJDD(1) + t341 * t390 - t309) * t333 + t370;
t290 = -g(3) * t333 + (t309 + 0.2e1 * t378) * t334;
t375 = qJDD(1) * t334;
t327 = t334 ^ 2;
t386 = t327 * t341;
t259 = -pkin(2) * t386 + pkin(7) * t375 + t290;
t336 = sin(qJ(3));
t391 = cos(qJ(3));
t217 = t258 * t391 - t259 * t336;
t218 = t258 * t336 + t259 * t391;
t371 = t334 * t391;
t381 = qJD(1) * t333;
t307 = -qJD(1) * t371 + t336 * t381;
t358 = t333 * t391 + t334 * t336;
t308 = t358 * qJD(1);
t263 = Ifges(4,4) * t308 - Ifges(4,2) * t307 + Ifges(4,6) * qJD(3);
t276 = -mrSges(5,2) * t307 - mrSges(5,3) * t308;
t376 = qJDD(1) * t333;
t380 = qJD(3) * t308;
t287 = -qJDD(1) * t371 + t336 * t376 + t380;
t379 = t307 * qJD(3);
t288 = qJDD(1) * t358 - t379;
t298 = mrSges(5,1) * t307 - qJD(3) * mrSges(5,3);
t274 = pkin(3) * t307 - qJ(4) * t308;
t340 = qJD(3) ^ 2;
t397 = -2 * qJD(4);
t210 = pkin(3) * t340 - qJDD(3) * qJ(4) + qJD(3) * t397 + t274 * t307 - t218;
t299 = mrSges(5,1) * t308 + qJD(3) * mrSges(5,2);
t300 = pkin(4) * t308 - qJD(3) * pkin(8);
t306 = t307 ^ 2;
t207 = -pkin(4) * t287 - pkin(8) * t306 + qJD(3) * t300 - t210;
t335 = sin(qJ(5));
t338 = cos(qJ(5));
t293 = qJD(3) * t338 + t307 * t335;
t240 = -qJD(5) * t293 - qJDD(3) * t335 + t287 * t338;
t292 = -qJD(3) * t335 + t307 * t338;
t241 = qJD(5) * t292 + qJDD(3) * t338 + t287 * t335;
t304 = qJD(5) + t308;
t251 = -mrSges(7,2) * t304 + mrSges(7,3) * t292;
t252 = -mrSges(6,2) * t304 + mrSges(6,3) * t292;
t255 = mrSges(6,1) * t304 - mrSges(6,3) * t293;
t253 = pkin(5) * t304 - qJ(6) * t293;
t291 = t292 ^ 2;
t200 = -pkin(5) * t240 - qJ(6) * t291 + t253 * t293 + qJDD(6) + t207;
t254 = mrSges(7,1) * t304 - mrSges(7,3) * t293;
t372 = m(7) * t200 + mrSges(7,2) * t241 + t254 * t293;
t392 = -m(6) * t207 - t241 * mrSges(6,2) + (mrSges(6,1) + mrSges(7,1)) * t240 - t293 * t255 + (t251 + t252) * t292 - t372;
t347 = -m(5) * t210 + qJDD(3) * mrSges(5,3) + qJD(3) * t299 - t392;
t326 = t333 ^ 2;
t315 = t337 * g(1) - g(2) * t339;
t365 = qJDD(2) - t315;
t286 = (-pkin(1) - t390) * qJDD(1) + (-qJ(2) + (-t326 - t327) * pkin(7)) * t341 + t365;
t343 = pkin(3) * t380 + t308 * t397 + (-t288 + t379) * qJ(4) + t286;
t202 = -t306 * pkin(4) - t308 * t300 + (pkin(3) + pkin(8)) * t287 + t343;
t212 = -qJDD(3) * pkin(3) - t340 * qJ(4) + t274 * t308 + qJDD(4) - t217;
t205 = (t307 * t308 - qJDD(3)) * pkin(8) + (t288 + t379) * pkin(4) + t212;
t197 = t202 * t338 + t205 * t335;
t221 = Ifges(7,5) * t293 + Ifges(7,6) * t292 + Ifges(7,3) * t304;
t222 = Ifges(6,5) * t293 + Ifges(6,6) * t292 + Ifges(6,3) * t304;
t226 = Ifges(6,1) * t293 + Ifges(6,4) * t292 + Ifges(6,5) * t304;
t285 = qJDD(5) + t288;
t194 = -pkin(5) * t291 + qJ(6) * t240 + 0.2e1 * qJD(6) * t292 - t253 * t304 + t197;
t225 = Ifges(7,1) * t293 + Ifges(7,4) * t292 + Ifges(7,5) * t304;
t356 = -mrSges(7,1) * t200 + mrSges(7,3) * t194 + Ifges(7,4) * t241 + Ifges(7,2) * t240 + Ifges(7,6) * t285 + t225 * t304;
t245 = -mrSges(7,1) * t292 + mrSges(7,2) * t293;
t373 = m(7) * t194 + mrSges(7,3) * t240 + t245 * t292;
t165 = Ifges(6,4) * t241 + Ifges(6,2) * t240 + Ifges(6,6) * t285 + t304 * t226 - mrSges(6,1) * t207 + mrSges(6,3) * t197 - pkin(5) * (-t240 * mrSges(7,1) - t292 * t251 + t372) + qJ(6) * (-t285 * mrSges(7,2) - t304 * t254 + t373) + (-t222 - t221) * t293 + t356;
t196 = -t335 * t202 + t205 * t338;
t191 = -0.2e1 * qJD(6) * t293 + (t292 * t304 - t241) * qJ(6) + (t292 * t293 + t285) * pkin(5) + t196;
t374 = m(7) * t191 + mrSges(7,1) * t285 + t251 * t304;
t188 = -t241 * mrSges(7,3) - t293 * t245 + t374;
t223 = Ifges(7,4) * t293 + Ifges(7,2) * t292 + Ifges(7,6) * t304;
t224 = Ifges(6,4) * t293 + Ifges(6,2) * t292 + Ifges(6,6) * t304;
t354 = mrSges(7,2) * t200 - mrSges(7,3) * t191 + Ifges(7,1) * t241 + Ifges(7,4) * t240 + Ifges(7,5) * t285 + t221 * t292;
t174 = mrSges(6,2) * t207 - mrSges(6,3) * t196 + Ifges(6,1) * t241 + Ifges(6,4) * t240 + Ifges(6,5) * t285 - qJ(6) * t188 + t292 * t222 + (-t223 - t224) * t304 + t354;
t246 = -mrSges(6,1) * t292 + mrSges(6,2) * t293;
t178 = m(6) * t196 + t285 * mrSges(6,1) + t304 * t252 + (-t245 - t246) * t293 + (-mrSges(6,3) - mrSges(7,3)) * t241 + t374;
t180 = m(6) * t197 + t240 * mrSges(6,3) + t292 * t246 + (-t254 - t255) * t304 + (-mrSges(6,2) - mrSges(7,2)) * t285 + t373;
t175 = t338 * t178 + t335 * t180;
t260 = Ifges(5,5) * qJD(3) - Ifges(5,6) * t308 + Ifges(5,3) * t307;
t351 = -mrSges(5,2) * t212 + mrSges(5,3) * t210 - Ifges(5,1) * qJDD(3) + Ifges(5,4) * t288 - Ifges(5,5) * t287 + pkin(8) * t175 + t335 * t165 - t174 * t338 + t260 * t308;
t352 = -m(5) * t212 - mrSges(5,1) * t288 - t276 * t308 - t175;
t262 = Ifges(5,4) * qJD(3) - Ifges(5,2) * t308 + Ifges(5,6) * t307;
t382 = Ifges(4,1) * t308 - Ifges(4,4) * t307 + Ifges(4,5) * qJD(3) - t262;
t399 = t382 * t307 - mrSges(4,2) * t218 + pkin(3) * (-qJDD(3) * mrSges(5,2) - qJD(3) * t298 + t352) + qJ(4) * (-t287 * mrSges(5,1) - t307 * t276 + t347) + mrSges(4,1) * t217 + t308 * t263 - Ifges(4,6) * t287 + Ifges(4,5) * t288 + Ifges(4,3) * qJDD(3) - t351;
t355 = -mrSges(7,1) * t191 + mrSges(7,2) * t194 - Ifges(7,5) * t241 - Ifges(7,6) * t240 - Ifges(7,3) * t285 - t223 * t293;
t398 = mrSges(6,1) * t196 - mrSges(6,2) * t197 + Ifges(6,5) * t241 + Ifges(6,6) * t240 + Ifges(6,3) * t285 + pkin(5) * t188 + t293 * t224 - t355 - (t226 + t225) * t292;
t209 = t287 * pkin(3) + t343;
t396 = mrSges(5,1) * t212 - mrSges(5,3) * t209 + pkin(4) * t175 + t398;
t275 = mrSges(4,1) * t307 + mrSges(4,2) * t308;
t296 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t307;
t169 = m(4) * t217 - t288 * mrSges(4,3) - t308 * t275 + (mrSges(4,1) - mrSges(5,2)) * qJDD(3) + (t296 - t298) * qJD(3) + t352;
t297 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t308;
t183 = -qJDD(3) * mrSges(4,2) + (-t275 - t276) * t307 + t347 + (-mrSges(4,3) - mrSges(5,1)) * t287 - qJD(3) * t297 + m(4) * t218;
t164 = t169 * t391 + t183 * t336;
t289 = -t333 * t309 + t370;
t363 = Ifges(3,4) * t333 + Ifges(3,2) * t334;
t364 = Ifges(3,1) * t333 + Ifges(3,4) * t334;
t394 = qJD(1) * t334;
t395 = -mrSges(3,1) * t289 + mrSges(3,2) * t290 - pkin(2) * t164 - (t363 * t381 - t364 * t394) * qJD(1) - t399;
t388 = -Ifges(5,6) - Ifges(4,4);
t387 = mrSges(3,2) * t333;
t176 = -t178 * t335 + t180 * t338;
t264 = Ifges(5,1) * qJD(3) - Ifges(5,4) * t308 + Ifges(5,5) * t307;
t383 = -Ifges(4,5) * t308 + Ifges(4,6) * t307 - Ifges(4,3) * qJD(3) - t264;
t360 = mrSges(3,3) * qJDD(1) + t341 * (-mrSges(3,1) * t334 + t387);
t162 = m(3) * t289 - t333 * t360 + t164;
t366 = -t336 * t169 + t183 * t391;
t163 = m(3) * t290 + t334 * t360 + t366;
t367 = -t162 * t333 + t163 * t334;
t362 = Ifges(3,5) * t333 + Ifges(3,6) * t334;
t170 = m(5) * t209 - mrSges(5,2) * t287 - mrSges(5,3) * t288 - t298 * t307 - t299 * t308 + t176;
t350 = -mrSges(5,1) * t210 + mrSges(5,2) * t209 - pkin(4) * t392 - pkin(8) * t176 - t338 * t165 - t335 * t174;
t156 = -mrSges(4,1) * t286 + mrSges(4,3) * t218 - pkin(3) * t170 + t383 * t308 - t388 * t288 + (-Ifges(4,2) - Ifges(5,3)) * t287 + (Ifges(4,6) - Ifges(5,5)) * qJDD(3) + t382 * qJD(3) + t350;
t157 = t383 * t307 + (Ifges(4,1) + Ifges(5,2)) * t288 + t388 * t287 + (Ifges(4,5) - Ifges(5,4)) * qJDD(3) + (t260 - t263) * qJD(3) + mrSges(4,2) * t286 - mrSges(4,3) * t217 - qJ(4) * t170 + t396;
t305 = -qJDD(1) * pkin(1) - t341 * qJ(2) + t365;
t311 = t362 * qJD(1);
t348 = m(4) * t286 + mrSges(4,1) * t287 + t288 * mrSges(4,2) + t296 * t307 + t308 * t297 + t170;
t152 = -mrSges(3,1) * t305 + mrSges(3,3) * t290 - pkin(2) * t348 + pkin(7) * t366 + qJDD(1) * t363 + t156 * t391 + t336 * t157 - t311 * t381;
t154 = mrSges(3,2) * t305 - mrSges(3,3) * t289 - pkin(7) * t164 + qJDD(1) * t364 - t336 * t156 + t157 * t391 + t311 * t394;
t346 = -m(3) * t305 + mrSges(3,1) * t375 - t348 + (t326 * t341 + t386) * mrSges(3,3);
t353 = -mrSges(2,2) * t316 + qJ(2) * t367 + t334 * t152 + t333 * t154 + pkin(1) * (-mrSges(3,2) * t376 + t346) + mrSges(2,1) * t315 + Ifges(2,3) * qJDD(1);
t166 = (mrSges(2,1) - t387) * qJDD(1) + t346 - t341 * mrSges(2,2) + m(2) * t315;
t160 = t162 * t334 + t163 * t333;
t158 = m(2) * t316 - mrSges(2,1) * t341 - qJDD(1) * mrSges(2,2) + t367;
t155 = (Ifges(2,6) - t362) * qJDD(1) + mrSges(2,1) * g(3) + t341 * Ifges(2,5) + mrSges(2,3) * t316 - pkin(1) * t160 + t395;
t150 = -mrSges(2,2) * g(3) - mrSges(2,3) * t315 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t341 - qJ(2) * t160 - t152 * t333 + t154 * t334;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t339 * t150 - t337 * t155 - pkin(6) * (t158 * t337 + t166 * t339), t150, t154, t157, -t307 * t262 - t351, t174, -t223 * t304 + t354; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t337 * t150 + t339 * t155 + pkin(6) * (t158 * t339 - t166 * t337), t155, t152, t156, Ifges(5,4) * qJDD(3) - Ifges(5,2) * t288 + Ifges(5,6) * t287 - qJD(3) * t260 + t307 * t264 - t396, t165, -t293 * t221 + t356; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t353, t353, qJDD(1) * t362 - t395, t399, Ifges(5,5) * qJDD(3) - Ifges(5,6) * t288 + Ifges(5,3) * t287 + qJD(3) * t262 + t308 * t264 - t350, t398, -t292 * t225 - t355;];
m_new  = t1;
