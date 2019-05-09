% Calculate vector of cutting torques with Newton-Euler for
% S6RPRPPR5
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta5]';
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
% Datum: 2019-05-05 17:01
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPRPPR5_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR5_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR5_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPPR5_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR5_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR5_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR5_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR5_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR5_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 16:57:48
% EndTime: 2019-05-05 16:58:10
% DurationCPUTime: 13.65s
% Computational Cost: add. (205172->366), mult. (509795->441), div. (0->0), fcn. (368539->10), ass. (0->147)
t338 = sin(qJ(1));
t340 = cos(qJ(1));
t315 = -t340 * g(1) - t338 * g(2);
t342 = qJD(1) ^ 2;
t308 = -t342 * pkin(1) + qJDD(1) * qJ(2) + t315;
t333 = sin(pkin(9));
t335 = cos(pkin(9));
t374 = qJD(1) * qJD(2);
t369 = -t335 * g(3) - 0.2e1 * t333 * t374;
t383 = pkin(2) * t335;
t250 = (-pkin(7) * qJDD(1) + t342 * t383 - t308) * t333 + t369;
t287 = -t333 * g(3) + (t308 + 0.2e1 * t374) * t335;
t371 = qJDD(1) * t335;
t326 = t335 ^ 2;
t380 = t326 * t342;
t259 = -pkin(2) * t380 + pkin(7) * t371 + t287;
t337 = sin(qJ(3));
t384 = cos(qJ(3));
t231 = t384 * t250 - t337 * t259;
t232 = t337 * t250 + t384 * t259;
t370 = t335 * t384;
t377 = qJD(1) * t333;
t306 = -qJD(1) * t370 + t337 * t377;
t358 = t384 * t333 + t335 * t337;
t307 = t358 * qJD(1);
t263 = Ifges(4,4) * t307 - Ifges(4,2) * t306 + Ifges(4,6) * qJD(3);
t272 = -t306 * mrSges(5,2) - t307 * mrSges(5,3);
t372 = qJDD(1) * t333;
t375 = t307 * qJD(3);
t284 = -qJDD(1) * t370 + t337 * t372 + t375;
t376 = t306 * qJD(3);
t285 = qJDD(1) * t358 - t376;
t298 = t306 * mrSges(5,1) - qJD(3) * mrSges(5,3);
t270 = t306 * pkin(3) - t307 * qJ(4);
t341 = qJD(3) ^ 2;
t388 = -2 * qJD(4);
t216 = t341 * pkin(3) - qJDD(3) * qJ(4) + qJD(3) * t388 + t306 * t270 - t232;
t297 = t307 * pkin(4) - qJD(3) * qJ(5);
t305 = t306 ^ 2;
t212 = -t284 * pkin(4) - t305 * qJ(5) + qJD(3) * t297 + qJDD(5) - t216;
t332 = sin(pkin(10));
t334 = cos(pkin(10));
t291 = -t332 * qJD(3) + t334 * t306;
t254 = -t307 * mrSges(6,2) + t291 * mrSges(6,3);
t292 = t334 * qJD(3) + t332 * t306;
t255 = t307 * mrSges(6,1) - t292 * mrSges(6,3);
t257 = -t332 * qJDD(3) + t334 * t284;
t258 = t334 * qJDD(3) + t332 * t284;
t256 = t307 * pkin(5) - t292 * pkin(8);
t290 = t291 ^ 2;
t205 = -t257 * pkin(5) - t290 * pkin(8) + t292 * t256 + t212;
t336 = sin(qJ(6));
t339 = cos(qJ(6));
t240 = t336 * t291 + t339 * t292;
t222 = -t240 * qJD(6) + t339 * t257 - t336 * t258;
t239 = t339 * t291 - t336 * t292;
t223 = t239 * qJD(6) + t336 * t257 + t339 * t258;
t303 = qJD(6) + t307;
t233 = -t303 * mrSges(7,2) + t239 * mrSges(7,3);
t234 = t303 * mrSges(7,1) - t240 * mrSges(7,3);
t356 = m(7) * t205 - t222 * mrSges(7,1) + t223 * mrSges(7,2) - t239 * t233 + t240 * t234;
t195 = -m(6) * t212 + t257 * mrSges(6,1) - t258 * mrSges(6,2) + t291 * t254 - t292 * t255 - t356;
t299 = t307 * mrSges(5,1) + qJD(3) * mrSges(5,2);
t349 = -m(5) * t216 + qJDD(3) * mrSges(5,3) + qJD(3) * t299 - t195;
t325 = t333 ^ 2;
t314 = t338 * g(1) - t340 * g(2);
t364 = qJDD(2) - t314;
t283 = (-pkin(1) - t383) * qJDD(1) + (-qJ(2) + (-t325 - t326) * pkin(7)) * t342 + t364;
t345 = pkin(3) * t375 + t307 * t388 + (-t285 + t376) * qJ(4) + t283;
t207 = -t305 * pkin(4) - t307 * t297 + (pkin(3) + qJ(5)) * t284 + t345;
t221 = -qJDD(3) * pkin(3) - t341 * qJ(4) + t307 * t270 + qJDD(4) - t231;
t210 = (t306 * t307 - qJDD(3)) * qJ(5) + (t285 + t376) * pkin(4) + t221;
t202 = -0.2e1 * qJD(5) * t292 - t332 * t207 + t334 * t210;
t199 = (t291 * t307 - t258) * pkin(8) + (t291 * t292 + t285) * pkin(5) + t202;
t203 = 0.2e1 * qJD(5) * t291 + t334 * t207 + t332 * t210;
t200 = -t290 * pkin(5) + t257 * pkin(8) - t307 * t256 + t203;
t198 = t336 * t199 + t339 * t200;
t224 = Ifges(7,5) * t240 + Ifges(7,6) * t239 + Ifges(7,3) * t303;
t226 = Ifges(7,1) * t240 + Ifges(7,4) * t239 + Ifges(7,5) * t303;
t282 = qJDD(6) + t285;
t183 = -mrSges(7,1) * t205 + mrSges(7,3) * t198 + Ifges(7,4) * t223 + Ifges(7,2) * t222 + Ifges(7,6) * t282 - t240 * t224 + t303 * t226;
t197 = t339 * t199 - t336 * t200;
t225 = Ifges(7,4) * t240 + Ifges(7,2) * t239 + Ifges(7,6) * t303;
t184 = mrSges(7,2) * t205 - mrSges(7,3) * t197 + Ifges(7,1) * t223 + Ifges(7,4) * t222 + Ifges(7,5) * t282 + t239 * t224 - t303 * t225;
t235 = Ifges(6,5) * t292 + Ifges(6,6) * t291 + Ifges(6,3) * t307;
t237 = Ifges(6,1) * t292 + Ifges(6,4) * t291 + Ifges(6,5) * t307;
t228 = -t239 * mrSges(7,1) + t240 * mrSges(7,2);
t191 = m(7) * t197 + t282 * mrSges(7,1) - t223 * mrSges(7,3) - t240 * t228 + t303 * t233;
t192 = m(7) * t198 - t282 * mrSges(7,2) + t222 * mrSges(7,3) + t239 * t228 - t303 * t234;
t365 = -t336 * t191 + t339 * t192;
t165 = -mrSges(6,1) * t212 + mrSges(6,3) * t203 + Ifges(6,4) * t258 + Ifges(6,2) * t257 + Ifges(6,6) * t285 - pkin(5) * t356 + pkin(8) * t365 + t339 * t183 + t336 * t184 - t292 * t235 + t307 * t237;
t182 = t339 * t191 + t336 * t192;
t236 = Ifges(6,4) * t292 + Ifges(6,2) * t291 + Ifges(6,6) * t307;
t167 = mrSges(6,2) * t212 - mrSges(6,3) * t202 + Ifges(6,1) * t258 + Ifges(6,4) * t257 + Ifges(6,5) * t285 - pkin(8) * t182 - t336 * t183 + t339 * t184 + t291 * t235 - t307 * t236;
t243 = -t291 * mrSges(6,1) + t292 * mrSges(6,2);
t179 = m(6) * t202 + t285 * mrSges(6,1) - t258 * mrSges(6,3) - t292 * t243 + t307 * t254 + t182;
t180 = m(6) * t203 - t285 * mrSges(6,2) + t257 * mrSges(6,3) + t291 * t243 - t307 * t255 + t365;
t175 = t334 * t179 + t332 * t180;
t260 = Ifges(5,5) * qJD(3) - Ifges(5,6) * t307 + Ifges(5,3) * t306;
t352 = -mrSges(5,2) * t221 + mrSges(5,3) * t216 - Ifges(5,1) * qJDD(3) + Ifges(5,4) * t285 - Ifges(5,5) * t284 + qJ(5) * t175 + t332 * t165 - t334 * t167 + t307 * t260;
t354 = -m(5) * t221 - t285 * mrSges(5,1) - t307 * t272 - t175;
t262 = Ifges(5,4) * qJD(3) - Ifges(5,2) * t307 + Ifges(5,6) * t306;
t378 = Ifges(4,1) * t307 - Ifges(4,4) * t306 + Ifges(4,5) * qJD(3) - t262;
t389 = -mrSges(4,2) * t232 + pkin(3) * (-qJDD(3) * mrSges(5,2) - qJD(3) * t298 + t354) + qJ(4) * (-t284 * mrSges(5,1) - t306 * t272 + t349) + mrSges(4,1) * t231 + t307 * t263 - Ifges(4,6) * t284 + Ifges(4,5) * t285 + Ifges(4,3) * qJDD(3) - t352 + t378 * t306;
t271 = t306 * mrSges(4,1) + t307 * mrSges(4,2);
t295 = -qJD(3) * mrSges(4,2) - t306 * mrSges(4,3);
t171 = m(4) * t231 - t285 * mrSges(4,3) - t307 * t271 + (mrSges(4,1) - mrSges(5,2)) * qJDD(3) + (t295 - t298) * qJD(3) + t354;
t296 = qJD(3) * mrSges(4,1) - t307 * mrSges(4,3);
t187 = -qJDD(3) * mrSges(4,2) - qJD(3) * t296 + m(4) * t232 + (-t271 - t272) * t306 + (-mrSges(4,3) - mrSges(5,1)) * t284 + t349;
t164 = t384 * t171 + t337 * t187;
t286 = -t333 * t308 + t369;
t362 = Ifges(3,4) * t333 + Ifges(3,2) * t335;
t363 = Ifges(3,1) * t333 + Ifges(3,4) * t335;
t386 = qJD(1) * t335;
t387 = -mrSges(3,1) * t286 + mrSges(3,2) * t287 - pkin(2) * t164 - (t362 * t377 - t363 * t386) * qJD(1) - t389;
t382 = -Ifges(5,6) - Ifges(4,4);
t381 = mrSges(3,2) * t333;
t176 = -t332 * t179 + t334 * t180;
t264 = Ifges(5,1) * qJD(3) - Ifges(5,4) * t307 + Ifges(5,5) * t306;
t379 = -Ifges(4,5) * t307 + Ifges(4,6) * t306 - Ifges(4,3) * qJD(3) - t264;
t359 = qJDD(1) * mrSges(3,3) + t342 * (-mrSges(3,1) * t335 + t381);
t162 = m(3) * t286 - t333 * t359 + t164;
t366 = -t337 * t171 + t384 * t187;
t163 = m(3) * t287 + t335 * t359 + t366;
t367 = -t333 * t162 + t335 * t163;
t361 = Ifges(3,5) * t333 + Ifges(3,6) * t335;
t215 = t284 * pkin(3) + t345;
t172 = m(5) * t215 - t284 * mrSges(5,2) - t285 * mrSges(5,3) - t306 * t298 - t307 * t299 + t176;
t351 = -mrSges(5,1) * t216 + mrSges(5,2) * t215 - pkin(4) * t195 - qJ(5) * t176 - t334 * t165 - t332 * t167;
t156 = -mrSges(4,1) * t283 + mrSges(4,3) * t232 - pkin(3) * t172 + t379 * t307 - t382 * t285 + (-Ifges(4,2) - Ifges(5,3)) * t284 + (Ifges(4,6) - Ifges(5,5)) * qJDD(3) + t378 * qJD(3) + t351;
t353 = -mrSges(7,1) * t197 + mrSges(7,2) * t198 - Ifges(7,5) * t223 - Ifges(7,6) * t222 - Ifges(7,3) * t282 - t240 * t225 + t239 * t226;
t348 = -mrSges(6,1) * t202 + mrSges(6,2) * t203 - Ifges(6,5) * t258 - Ifges(6,6) * t257 - Ifges(6,3) * t285 - pkin(5) * t182 - t292 * t236 + t291 * t237 + t353;
t344 = -mrSges(5,1) * t221 + mrSges(5,3) * t215 - pkin(4) * t175 + t348;
t157 = mrSges(4,2) * t283 - mrSges(4,3) * t231 - qJ(4) * t172 - t344 + t379 * t306 + (Ifges(4,1) + Ifges(5,2)) * t285 + t382 * t284 + (Ifges(4,5) - Ifges(5,4)) * qJDD(3) + (t260 - t263) * qJD(3);
t304 = -qJDD(1) * pkin(1) - t342 * qJ(2) + t364;
t310 = t361 * qJD(1);
t350 = m(4) * t283 + t284 * mrSges(4,1) + t285 * mrSges(4,2) + t306 * t295 + t307 * t296 + t172;
t152 = -mrSges(3,1) * t304 + mrSges(3,3) * t287 - pkin(2) * t350 + pkin(7) * t366 + t362 * qJDD(1) + t384 * t156 + t337 * t157 - t310 * t377;
t154 = mrSges(3,2) * t304 - mrSges(3,3) * t286 - pkin(7) * t164 + t363 * qJDD(1) - t337 * t156 + t384 * t157 + t310 * t386;
t347 = -m(3) * t304 + mrSges(3,1) * t371 - t350 + (t325 * t342 + t380) * mrSges(3,3);
t355 = -mrSges(2,2) * t315 + qJ(2) * t367 + t335 * t152 + t333 * t154 + pkin(1) * (-mrSges(3,2) * t372 + t347) + mrSges(2,1) * t314 + Ifges(2,3) * qJDD(1);
t168 = -t342 * mrSges(2,2) + m(2) * t314 + (mrSges(2,1) - t381) * qJDD(1) + t347;
t160 = t335 * t162 + t333 * t163;
t158 = m(2) * t315 - t342 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t367;
t155 = mrSges(2,1) * g(3) + t342 * Ifges(2,5) + mrSges(2,3) * t315 - pkin(1) * t160 + (Ifges(2,6) - t361) * qJDD(1) + t387;
t150 = -mrSges(2,2) * g(3) - mrSges(2,3) * t314 + Ifges(2,5) * qJDD(1) - t342 * Ifges(2,6) - qJ(2) * t160 - t333 * t152 + t335 * t154;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t340 * t150 - t338 * t155 - pkin(6) * (t338 * t158 + t340 * t168), t150, t154, t157, -t306 * t262 - t352, t167, t184; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t338 * t150 + t340 * t155 + pkin(6) * (t340 * t158 - t338 * t168), t155, t152, t156, Ifges(5,4) * qJDD(3) - Ifges(5,2) * t285 + Ifges(5,6) * t284 - qJD(3) * t260 + t306 * t264 + t344, t165, t183; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t355, t355, qJDD(1) * t361 - t387, t389, Ifges(5,5) * qJDD(3) - Ifges(5,6) * t285 + Ifges(5,3) * t284 + qJD(3) * t262 + t307 * t264 - t351, -t348, -t353;];
m_new  = t1;
