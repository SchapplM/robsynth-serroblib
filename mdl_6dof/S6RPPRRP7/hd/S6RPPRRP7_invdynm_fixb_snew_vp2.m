% Calculate vector of cutting torques with Newton-Euler for
% S6RPPRRP7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
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
% Datum: 2019-05-05 15:05
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPPRRP7_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP7_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP7_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRP7_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP7_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP7_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP7_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP7_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRP7_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 15:03:18
% EndTime: 2019-05-05 15:03:32
% DurationCPUTime: 7.00s
% Computational Cost: add. (88821->315), mult. (199468->371), div. (0->0), fcn. (135652->8), ass. (0->129)
t301 = qJD(1) ^ 2;
t292 = sin(pkin(9));
t283 = t292 ^ 2;
t293 = cos(pkin(9));
t339 = t293 ^ 2 + t283;
t328 = t339 * mrSges(4,3);
t351 = t301 * t328;
t296 = sin(qJ(1));
t299 = cos(qJ(1));
t272 = t296 * g(1) - t299 * g(2);
t318 = -t301 * qJ(2) + qJDD(2) - t272;
t343 = -pkin(1) - qJ(3);
t350 = -(2 * qJD(1) * qJD(3)) + t343 * qJDD(1) + t318;
t273 = -t299 * g(1) - t296 * g(2);
t349 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t273;
t295 = sin(qJ(4));
t298 = cos(qJ(4));
t322 = t292 * t298 + t293 * t295;
t265 = t322 * qJD(1);
t321 = -t292 * t295 + t293 * t298;
t337 = t265 * qJD(4);
t248 = t321 * qJDD(1) - t337;
t266 = t321 * qJD(1);
t294 = sin(qJ(5));
t297 = cos(qJ(5));
t251 = t297 * qJD(4) - t294 * t266;
t212 = t251 * qJD(5) + t294 * qJDD(4) + t297 * t248;
t252 = t294 * qJD(4) + t297 * t266;
t216 = -t251 * mrSges(7,1) + t252 * mrSges(7,2);
t244 = t292 * g(3) + t350 * t293;
t347 = pkin(3) * t301;
t221 = (-pkin(7) * qJDD(1) - t292 * t347) * t293 + t244;
t245 = -t293 * g(3) + t350 * t292;
t334 = qJDD(1) * t292;
t222 = -pkin(7) * t334 - t283 * t347 + t245;
t189 = t295 * t221 + t298 * t222;
t246 = t265 * pkin(4) - t266 * pkin(8);
t300 = qJD(4) ^ 2;
t183 = -t300 * pkin(4) + qJDD(4) * pkin(8) - t265 * t246 + t189;
t314 = qJDD(3) + t349;
t229 = pkin(3) * t334 + (-t339 * pkin(7) + t343) * t301 + t314;
t336 = t266 * qJD(4);
t247 = -t322 * qJDD(1) - t336;
t186 = (-t248 + t337) * pkin(8) + (-t247 + t336) * pkin(4) + t229;
t179 = -t294 * t183 + t297 * t186;
t243 = qJDD(5) - t247;
t263 = qJD(5) + t265;
t173 = -0.2e1 * qJD(6) * t252 + (t251 * t263 - t212) * qJ(6) + (t251 * t252 + t243) * pkin(5) + t179;
t223 = -t263 * mrSges(7,2) + t251 * mrSges(7,3);
t331 = m(7) * t173 + t243 * mrSges(7,1) + t263 * t223;
t170 = -t212 * mrSges(7,3) - t252 * t216 + t331;
t180 = t297 * t183 + t294 * t186;
t197 = Ifges(6,4) * t252 + Ifges(6,2) * t251 + Ifges(6,6) * t263;
t198 = Ifges(7,1) * t252 + Ifges(7,4) * t251 + Ifges(7,5) * t263;
t199 = Ifges(6,1) * t252 + Ifges(6,4) * t251 + Ifges(6,5) * t263;
t211 = -t252 * qJD(5) + t297 * qJDD(4) - t294 * t248;
t225 = t263 * pkin(5) - t252 * qJ(6);
t250 = t251 ^ 2;
t176 = -t250 * pkin(5) + t211 * qJ(6) + 0.2e1 * qJD(6) * t251 - t263 * t225 + t180;
t196 = Ifges(7,4) * t252 + Ifges(7,2) * t251 + Ifges(7,6) * t263;
t316 = -mrSges(7,1) * t173 + mrSges(7,2) * t176 - Ifges(7,5) * t212 - Ifges(7,6) * t211 - Ifges(7,3) * t243 - t252 * t196;
t348 = mrSges(6,1) * t179 - mrSges(6,2) * t180 + Ifges(6,5) * t212 + Ifges(6,6) * t211 + Ifges(6,3) * t243 + pkin(5) * t170 + t252 * t197 - (t199 + t198) * t251 - t316;
t346 = mrSges(2,1) - mrSges(3,2);
t345 = -mrSges(6,2) - mrSges(7,2);
t344 = -Ifges(2,6) + Ifges(3,5);
t342 = Ifges(4,6) * t292;
t239 = t265 * mrSges(5,1) + t266 * mrSges(5,2);
t257 = qJD(4) * mrSges(5,1) - t266 * mrSges(5,3);
t217 = -t251 * mrSges(6,1) + t252 * mrSges(6,2);
t224 = -t263 * mrSges(6,2) + t251 * mrSges(6,3);
t162 = m(6) * t179 + t243 * mrSges(6,1) + t263 * t224 + (-t216 - t217) * t252 + (-mrSges(6,3) - mrSges(7,3)) * t212 + t331;
t330 = m(7) * t176 + t211 * mrSges(7,3) + t251 * t216;
t226 = t263 * mrSges(7,1) - t252 * mrSges(7,3);
t340 = -t263 * mrSges(6,1) + t252 * mrSges(6,3) - t226;
t165 = m(6) * t180 + t211 * mrSges(6,3) + t251 * t217 + t345 * t243 + t340 * t263 + t330;
t326 = -t294 * t162 + t297 * t165;
t155 = m(5) * t189 - qJDD(4) * mrSges(5,2) + t247 * mrSges(5,3) - qJD(4) * t257 - t265 * t239 + t326;
t188 = t298 * t221 - t295 * t222;
t256 = -qJD(4) * mrSges(5,2) - t265 * mrSges(5,3);
t182 = -qJDD(4) * pkin(4) - t300 * pkin(8) + t266 * t246 - t188;
t178 = -t211 * pkin(5) - t250 * qJ(6) + t252 * t225 + qJDD(6) + t182;
t325 = -m(7) * t178 + t211 * mrSges(7,1) + t251 * t223;
t307 = -m(6) * t182 + t211 * mrSges(6,1) + t345 * t212 + t251 * t224 + t340 * t252 + t325;
t167 = m(5) * t188 + qJDD(4) * mrSges(5,1) - t248 * mrSges(5,3) + qJD(4) * t256 - t266 * t239 + t307;
t147 = t295 * t155 + t298 * t167;
t159 = t297 * t162 + t294 * t165;
t338 = t301 * (Ifges(4,5) * t293 - t342);
t333 = qJDD(1) * t293;
t329 = Ifges(3,4) + t342;
t320 = -qJDD(1) * mrSges(4,3) - t301 * (mrSges(4,1) * t292 + mrSges(4,2) * t293);
t144 = m(4) * t244 + t320 * t293 + t147;
t327 = t298 * t155 - t295 * t167;
t145 = m(4) * t245 + t320 * t292 + t327;
t141 = -t292 * t144 + t293 * t145;
t324 = Ifges(4,1) * t293 - Ifges(4,4) * t292;
t323 = Ifges(4,4) * t293 - Ifges(4,2) * t292;
t140 = t293 * t144 + t292 * t145;
t317 = -mrSges(7,1) * t178 + mrSges(7,3) * t176 + Ifges(7,4) * t212 + Ifges(7,2) * t211 + Ifges(7,6) * t243 + t263 * t198;
t194 = Ifges(7,5) * t252 + Ifges(7,6) * t251 + Ifges(7,3) * t263;
t315 = mrSges(7,2) * t178 - mrSges(7,3) * t173 + Ifges(7,1) * t212 + Ifges(7,4) * t211 + Ifges(7,5) * t243 + t251 * t194;
t264 = -qJDD(1) * pkin(1) + t318;
t313 = -m(3) * t264 + t301 * mrSges(3,3) - t140;
t312 = -m(5) * t229 + t247 * mrSges(5,1) - t248 * mrSges(5,2) - t265 * t256 - t266 * t257 - t159;
t195 = Ifges(6,5) * t252 + Ifges(6,6) * t251 + Ifges(6,3) * t263;
t149 = Ifges(6,4) * t212 + Ifges(6,2) * t211 + Ifges(6,6) * t243 + t263 * t199 - mrSges(6,1) * t182 + mrSges(6,3) * t180 - pkin(5) * (t212 * mrSges(7,2) - t325) + qJ(6) * (-t243 * mrSges(7,2) - t263 * t226 + t330) + (-pkin(5) * t226 - t194 - t195) * t252 + t317;
t157 = mrSges(6,2) * t182 - mrSges(6,3) * t179 + Ifges(6,1) * t212 + Ifges(6,4) * t211 + Ifges(6,5) * t243 - qJ(6) * t170 + t251 * t195 + (-t196 - t197) * t263 + t315;
t230 = Ifges(5,5) * t266 - Ifges(5,6) * t265 + Ifges(5,3) * qJD(4);
t231 = Ifges(5,4) * t266 - Ifges(5,2) * t265 + Ifges(5,6) * qJD(4);
t136 = mrSges(5,2) * t229 - mrSges(5,3) * t188 + Ifges(5,1) * t248 + Ifges(5,4) * t247 + Ifges(5,5) * qJDD(4) - pkin(8) * t159 - qJD(4) * t231 - t294 * t149 + t297 * t157 - t265 * t230;
t232 = Ifges(5,1) * t266 - Ifges(5,4) * t265 + Ifges(5,5) * qJD(4);
t142 = -mrSges(5,1) * t229 + mrSges(5,3) * t189 + Ifges(5,4) * t248 + Ifges(5,2) * t247 + Ifges(5,6) * qJDD(4) - pkin(4) * t159 + qJD(4) * t232 - t266 * t230 - t348;
t255 = t343 * t301 + t314;
t133 = -mrSges(4,1) * t255 + mrSges(4,3) * t245 + pkin(3) * t312 + pkin(7) * t327 + t323 * qJDD(1) + t295 * t136 + t298 * t142 - t293 * t338;
t135 = mrSges(4,2) * t255 - mrSges(4,3) * t244 - pkin(7) * t147 + t324 * qJDD(1) + t298 * t136 - t295 * t142 - t292 * t338;
t259 = t301 * pkin(1) - t349;
t311 = mrSges(3,2) * t264 - mrSges(3,3) * t259 + Ifges(3,1) * qJDD(1) - qJ(3) * t140 - t292 * t133 + t293 * t135;
t309 = -m(4) * t255 - mrSges(4,1) * t334 - mrSges(4,2) * t333 + t312;
t310 = -mrSges(3,1) * t259 - pkin(2) * (t309 + t351) - qJ(3) * t141 - t293 * t133 - t292 * t135;
t308 = -mrSges(5,1) * t188 + mrSges(5,2) * t189 - Ifges(5,5) * t248 - Ifges(5,6) * t247 - Ifges(5,3) * qJDD(4) - pkin(4) * t307 - pkin(8) * t326 - t297 * t149 - t294 * t157 - t266 * t231 - t265 * t232;
t305 = -m(3) * t259 + t301 * mrSges(3,2) + qJDD(1) * mrSges(3,3) - t309;
t306 = -mrSges(2,2) * t273 + pkin(1) * (-qJDD(1) * mrSges(3,2) + t313) + qJ(2) * (t305 - t351) + mrSges(2,1) * t272 + Ifges(2,3) * qJDD(1) + t311;
t304 = -mrSges(4,1) * t244 + mrSges(4,2) * t245 - Ifges(4,5) * t333 - pkin(3) * t147 + t308 + (-t292 * t324 - t293 * t323) * t301;
t302 = -mrSges(3,1) * t264 - pkin(2) * t140 + t304;
t150 = t305 - qJDD(1) * mrSges(2,2) + (-mrSges(2,1) - t328) * t301 + m(2) * t273;
t139 = -m(3) * g(3) + t141;
t137 = m(2) * t272 - t301 * mrSges(2,2) + t346 * qJDD(1) + t313;
t132 = -t302 - qJ(2) * t139 + t344 * t301 + (Ifges(2,5) - t329) * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3) - mrSges(2,3) * t272;
t131 = mrSges(2,3) * t273 - pkin(1) * t139 + (-Ifges(3,4) + Ifges(2,5)) * t301 - t344 * qJDD(1) + t346 * g(3) + t310;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t299 * t132 - t296 * t131 - pkin(6) * (t299 * t137 + t296 * t150), t132, t311, t135, t136, t157, -t263 * t196 + t315; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t296 * t132 + t299 * t131 + pkin(6) * (-t296 * t137 + t299 * t150), t131, -mrSges(3,3) * g(3) - t301 * Ifges(3,5) + t329 * qJDD(1) + t302, t133, t142, t149, -t252 * t194 + t317; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t306, t306, mrSges(3,2) * g(3) + t301 * Ifges(3,4) + Ifges(3,5) * qJDD(1) - t310, -Ifges(4,6) * t334 - t304, -t308, t348, -t251 * t198 - t316;];
m_new  = t1;
