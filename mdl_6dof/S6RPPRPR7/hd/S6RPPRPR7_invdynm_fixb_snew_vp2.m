% Calculate vector of cutting torques with Newton-Euler for
% S6RPPRPR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
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
% Datum: 2019-05-05 14:33
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPPRPR7_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR7_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR7_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRPR7_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR7_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR7_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR7_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR7_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR7_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:31:18
% EndTime: 2019-05-05 14:31:37
% DurationCPUTime: 11.12s
% Computational Cost: add. (190634->320), mult. (448244->391), div. (0->0), fcn. (315836->10), ass. (0->137)
t303 = qJD(1) ^ 2;
t293 = sin(pkin(9));
t283 = t293 ^ 2;
t295 = cos(pkin(9));
t338 = t295 ^ 2 + t283;
t329 = t338 * mrSges(4,3);
t346 = t303 * t329;
t298 = sin(qJ(1));
t301 = cos(qJ(1));
t272 = t298 * g(1) - t301 * g(2);
t319 = -t303 * qJ(2) + qJDD(2) - t272;
t340 = -pkin(1) - qJ(3);
t345 = -(2 * qJD(1) * qJD(3)) + t340 * qJDD(1) + t319;
t273 = -t301 * g(1) - t298 * g(2);
t344 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t273;
t343 = pkin(3) * t303;
t342 = mrSges(2,1) - mrSges(3,2);
t341 = -Ifges(2,6) + Ifges(3,5);
t339 = Ifges(4,6) * t293;
t243 = t293 * g(3) + t345 * t295;
t221 = (-pkin(7) * qJDD(1) - t293 * t343) * t295 + t243;
t244 = -t295 * g(3) + t345 * t293;
t333 = qJDD(1) * t293;
t222 = -pkin(7) * t333 - t283 * t343 + t244;
t297 = sin(qJ(4));
t300 = cos(qJ(4));
t206 = t297 * t221 + t300 * t222;
t323 = t293 * t300 + t295 * t297;
t265 = t323 * qJD(1);
t322 = -t293 * t297 + t295 * t300;
t266 = t322 * qJD(1);
t238 = t265 * mrSges(5,1) + t266 * mrSges(5,2);
t335 = t266 * qJD(4);
t245 = t323 * qJDD(1) + t335;
t257 = qJD(4) * mrSges(5,1) - t266 * mrSges(5,3);
t237 = t265 * pkin(4) - t266 * qJ(5);
t302 = qJD(4) ^ 2;
t190 = -t302 * pkin(4) + qJDD(4) * qJ(5) - t265 * t237 + t206;
t318 = qJDD(3) + t344;
t230 = pkin(3) * t333 + (-t338 * pkin(7) + t340) * t303 + t318;
t336 = t265 * qJD(4);
t246 = t322 * qJDD(1) - t336;
t198 = (-t246 + t336) * qJ(5) + (t245 + t335) * pkin(4) + t230;
t292 = sin(pkin(10));
t294 = cos(pkin(10));
t251 = t292 * qJD(4) + t294 * t266;
t184 = -0.2e1 * qJD(5) * t251 - t292 * t190 + t294 * t198;
t228 = t292 * qJDD(4) + t294 * t246;
t250 = t294 * qJD(4) - t292 * t266;
t182 = (t250 * t265 - t228) * pkin(8) + (t250 * t251 + t245) * pkin(5) + t184;
t185 = 0.2e1 * qJD(5) * t250 + t294 * t190 + t292 * t198;
t226 = t265 * pkin(5) - t251 * pkin(8);
t227 = t294 * qJDD(4) - t292 * t246;
t249 = t250 ^ 2;
t183 = -t249 * pkin(5) + t227 * pkin(8) - t265 * t226 + t185;
t296 = sin(qJ(6));
t299 = cos(qJ(6));
t180 = t299 * t182 - t296 * t183;
t213 = t299 * t250 - t296 * t251;
t197 = t213 * qJD(6) + t296 * t227 + t299 * t228;
t214 = t296 * t250 + t299 * t251;
t203 = -t213 * mrSges(7,1) + t214 * mrSges(7,2);
t263 = qJD(6) + t265;
t207 = -t263 * mrSges(7,2) + t213 * mrSges(7,3);
t242 = qJDD(6) + t245;
t175 = m(7) * t180 + t242 * mrSges(7,1) - t197 * mrSges(7,3) - t214 * t203 + t263 * t207;
t181 = t296 * t182 + t299 * t183;
t196 = -t214 * qJD(6) + t299 * t227 - t296 * t228;
t208 = t263 * mrSges(7,1) - t214 * mrSges(7,3);
t176 = m(7) * t181 - t242 * mrSges(7,2) + t196 * mrSges(7,3) + t213 * t203 - t263 * t208;
t167 = t299 * t175 + t296 * t176;
t216 = -t250 * mrSges(6,1) + t251 * mrSges(6,2);
t224 = -t265 * mrSges(6,2) + t250 * mrSges(6,3);
t165 = m(6) * t184 + t245 * mrSges(6,1) - t228 * mrSges(6,3) - t251 * t216 + t265 * t224 + t167;
t225 = t265 * mrSges(6,1) - t251 * mrSges(6,3);
t326 = -t296 * t175 + t299 * t176;
t166 = m(6) * t185 - t245 * mrSges(6,2) + t227 * mrSges(6,3) + t250 * t216 - t265 * t225 + t326;
t327 = -t292 * t165 + t294 * t166;
t158 = m(5) * t206 - qJDD(4) * mrSges(5,2) - t245 * mrSges(5,3) - qJD(4) * t257 - t265 * t238 + t327;
t205 = t300 * t221 - t297 * t222;
t256 = -qJD(4) * mrSges(5,2) - t265 * mrSges(5,3);
t189 = -qJDD(4) * pkin(4) - t302 * qJ(5) + t266 * t237 + qJDD(5) - t205;
t186 = -t227 * pkin(5) - t249 * pkin(8) + t251 * t226 + t189;
t316 = m(7) * t186 - t196 * mrSges(7,1) + t197 * mrSges(7,2) - t213 * t207 + t214 * t208;
t307 = -m(6) * t189 + t227 * mrSges(6,1) - t228 * mrSges(6,2) + t250 * t224 - t251 * t225 - t316;
t171 = m(5) * t205 + qJDD(4) * mrSges(5,1) - t246 * mrSges(5,3) + qJD(4) * t256 - t266 * t238 + t307;
t148 = t297 * t158 + t300 * t171;
t160 = t294 * t165 + t292 * t166;
t337 = t303 * (Ifges(4,5) * t295 - t339);
t332 = qJDD(1) * t295;
t330 = Ifges(3,4) + t339;
t321 = -qJDD(1) * mrSges(4,3) - t303 * (mrSges(4,1) * t293 + mrSges(4,2) * t295);
t145 = m(4) * t243 + t321 * t295 + t148;
t328 = t300 * t158 - t297 * t171;
t146 = m(4) * t244 + t321 * t293 + t328;
t142 = -t293 * t145 + t295 * t146;
t325 = Ifges(4,1) * t295 - Ifges(4,4) * t293;
t324 = Ifges(4,4) * t295 - Ifges(4,2) * t293;
t141 = t295 * t145 + t293 * t146;
t264 = -qJDD(1) * pkin(1) + t319;
t317 = -m(3) * t264 + t303 * mrSges(3,3) - t141;
t315 = m(5) * t230 + t245 * mrSges(5,1) + t246 * mrSges(5,2) + t265 * t256 + t266 * t257 + t160;
t200 = Ifges(7,4) * t214 + Ifges(7,2) * t213 + Ifges(7,6) * t263;
t201 = Ifges(7,1) * t214 + Ifges(7,4) * t213 + Ifges(7,5) * t263;
t314 = -mrSges(7,1) * t180 + mrSges(7,2) * t181 - Ifges(7,5) * t197 - Ifges(7,6) * t196 - Ifges(7,3) * t242 - t214 * t200 + t213 * t201;
t199 = Ifges(7,5) * t214 + Ifges(7,6) * t213 + Ifges(7,3) * t263;
t168 = -mrSges(7,1) * t186 + mrSges(7,3) * t181 + Ifges(7,4) * t197 + Ifges(7,2) * t196 + Ifges(7,6) * t242 - t214 * t199 + t263 * t201;
t169 = mrSges(7,2) * t186 - mrSges(7,3) * t180 + Ifges(7,1) * t197 + Ifges(7,4) * t196 + Ifges(7,5) * t242 + t213 * t199 - t263 * t200;
t209 = Ifges(6,5) * t251 + Ifges(6,6) * t250 + Ifges(6,3) * t265;
t211 = Ifges(6,1) * t251 + Ifges(6,4) * t250 + Ifges(6,5) * t265;
t150 = -mrSges(6,1) * t189 + mrSges(6,3) * t185 + Ifges(6,4) * t228 + Ifges(6,2) * t227 + Ifges(6,6) * t245 - pkin(5) * t316 + pkin(8) * t326 + t299 * t168 + t296 * t169 - t251 * t209 + t265 * t211;
t210 = Ifges(6,4) * t251 + Ifges(6,2) * t250 + Ifges(6,6) * t265;
t152 = mrSges(6,2) * t189 - mrSges(6,3) * t184 + Ifges(6,1) * t228 + Ifges(6,4) * t227 + Ifges(6,5) * t245 - pkin(8) * t167 - t296 * t168 + t299 * t169 + t250 * t209 - t265 * t210;
t231 = Ifges(5,5) * t266 - Ifges(5,6) * t265 + Ifges(5,3) * qJD(4);
t232 = Ifges(5,4) * t266 - Ifges(5,2) * t265 + Ifges(5,6) * qJD(4);
t137 = mrSges(5,2) * t230 - mrSges(5,3) * t205 + Ifges(5,1) * t246 - Ifges(5,4) * t245 + Ifges(5,5) * qJDD(4) - qJ(5) * t160 - qJD(4) * t232 - t292 * t150 + t294 * t152 - t265 * t231;
t233 = Ifges(5,1) * t266 - Ifges(5,4) * t265 + Ifges(5,5) * qJD(4);
t304 = -mrSges(6,1) * t184 + mrSges(6,2) * t185 - Ifges(6,5) * t228 - Ifges(6,6) * t227 - pkin(5) * t167 - t251 * t210 + t250 * t211 + t314;
t143 = t304 + (-Ifges(5,2) - Ifges(6,3)) * t245 - t266 * t231 + Ifges(5,4) * t246 - mrSges(5,1) * t230 + qJD(4) * t233 + mrSges(5,3) * t206 + Ifges(5,6) * qJDD(4) - pkin(4) * t160;
t255 = t340 * t303 + t318;
t134 = -mrSges(4,1) * t255 + mrSges(4,3) * t244 - pkin(3) * t315 + pkin(7) * t328 + t324 * qJDD(1) + t297 * t137 + t300 * t143 - t295 * t337;
t136 = mrSges(4,2) * t255 - mrSges(4,3) * t243 - pkin(7) * t148 + t325 * qJDD(1) + t300 * t137 - t297 * t143 - t293 * t337;
t259 = t303 * pkin(1) - t344;
t313 = mrSges(3,2) * t264 - mrSges(3,3) * t259 + Ifges(3,1) * qJDD(1) - qJ(3) * t141 - t293 * t134 + t295 * t136;
t311 = -m(4) * t255 - mrSges(4,1) * t333 - mrSges(4,2) * t332 - t315;
t312 = -mrSges(3,1) * t259 - pkin(2) * (t311 + t346) - qJ(3) * t142 - t295 * t134 - t293 * t136;
t310 = -mrSges(5,1) * t205 + mrSges(5,2) * t206 - Ifges(5,5) * t246 + Ifges(5,6) * t245 - Ifges(5,3) * qJDD(4) - pkin(4) * t307 - qJ(5) * t327 - t294 * t150 - t292 * t152 - t266 * t232 - t265 * t233;
t308 = -m(3) * t259 + t303 * mrSges(3,2) + qJDD(1) * mrSges(3,3) - t311;
t309 = -mrSges(2,2) * t273 + pkin(1) * (-qJDD(1) * mrSges(3,2) + t317) + qJ(2) * (t308 - t346) + mrSges(2,1) * t272 + Ifges(2,3) * qJDD(1) + t313;
t306 = -mrSges(4,1) * t243 + mrSges(4,2) * t244 - Ifges(4,5) * t332 - pkin(3) * t148 + t310 + (-t293 * t325 - t295 * t324) * t303;
t305 = -mrSges(3,1) * t264 - pkin(2) * t141 + t306;
t153 = t308 + m(2) * t273 + (-mrSges(2,1) - t329) * t303 - qJDD(1) * mrSges(2,2);
t140 = -m(3) * g(3) + t142;
t138 = m(2) * t272 - t303 * mrSges(2,2) + t342 * qJDD(1) + t317;
t133 = -t305 + t341 * t303 + (Ifges(2,5) - t330) * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3) - mrSges(2,3) * t272 - qJ(2) * t140;
t132 = mrSges(2,3) * t273 - pkin(1) * t140 + (-Ifges(3,4) + Ifges(2,5)) * t303 - t341 * qJDD(1) + t342 * g(3) + t312;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t301 * t133 - t298 * t132 - pkin(6) * (t301 * t138 + t298 * t153), t133, t313, t136, t137, t152, t169; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t298 * t133 + t301 * t132 + pkin(6) * (-t298 * t138 + t301 * t153), t132, -mrSges(3,3) * g(3) - t303 * Ifges(3,5) + t330 * qJDD(1) + t305, t134, t143, t150, t168; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t309, t309, mrSges(3,2) * g(3) + t303 * Ifges(3,4) + Ifges(3,5) * qJDD(1) - t312, -Ifges(4,6) * t333 - t306, -t310, Ifges(6,3) * t245 - t304, -t314;];
m_new  = t1;
