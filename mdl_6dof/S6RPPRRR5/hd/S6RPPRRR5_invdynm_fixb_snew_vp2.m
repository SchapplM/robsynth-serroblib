% Calculate vector of cutting torques with Newton-Euler for
% S6RPPRRR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6]';
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
% Datum: 2019-05-05 15:51
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPPRRR5_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR5_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR5_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRR5_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR5_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR5_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR5_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR5_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR5_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 15:50:00
% EndTime: 2019-05-05 15:50:10
% DurationCPUTime: 5.31s
% Computational Cost: add. (72142->299), mult. (139844->355), div. (0->0), fcn. (82275->8), ass. (0->116)
t327 = 2 * qJD(1);
t286 = sin(qJ(1));
t290 = cos(qJ(1));
t255 = -t290 * g(1) - t286 * g(2);
t326 = qJDD(1) * qJ(2) + (qJD(2) * t327) + t255;
t254 = t286 * g(1) - t290 * g(2);
t292 = qJD(1) ^ 2;
t233 = -qJDD(1) * pkin(1) - t292 * qJ(2) + qJDD(2) - t254;
t307 = qJDD(1) * qJ(3) + (qJD(3) * t327) - t233;
t284 = sin(qJ(5));
t285 = sin(qJ(4));
t288 = cos(qJ(5));
t289 = cos(qJ(4));
t238 = (t284 * t289 + t285 * t288) * qJD(1);
t325 = mrSges(3,2) - mrSges(4,3);
t324 = (Ifges(3,4) - Ifges(4,5));
t323 = Ifges(2,6) - Ifges(3,5);
t322 = mrSges(4,3) * t292;
t226 = qJDD(3) + (-pkin(1) - qJ(3)) * t292 + t326;
t219 = -qJDD(1) * pkin(7) + t226;
t211 = t285 * g(3) + t289 * t219;
t319 = qJD(1) * qJD(4);
t316 = t285 * t319;
t249 = qJDD(1) * t289 - t316;
t188 = (-t249 - t316) * pkin(8) + (-t285 * t289 * t292 + qJDD(4)) * pkin(4) + t211;
t212 = -g(3) * t289 + t285 * t219;
t248 = -qJDD(1) * t285 - t289 * t319;
t320 = qJD(1) * t289;
t253 = qJD(4) * pkin(4) - pkin(8) * t320;
t277 = t285 ^ 2;
t189 = -pkin(4) * t277 * t292 + pkin(8) * t248 - qJD(4) * t253 + t212;
t178 = t284 * t188 + t288 * t189;
t239 = (-t284 * t285 + t288 * t289) * qJD(1);
t200 = -qJD(5) * t239 + t248 * t288 - t249 * t284;
t213 = mrSges(6,1) * t238 + mrSges(6,2) * t239;
t264 = qJD(4) + qJD(5);
t228 = mrSges(6,1) * t264 - mrSges(6,3) * t239;
t263 = qJDD(4) + qJDD(5);
t214 = pkin(5) * t238 - pkin(9) * t239;
t262 = t264 ^ 2;
t173 = -pkin(5) * t262 + pkin(9) * t263 - t214 * t238 + t178;
t191 = -t248 * pkin(4) + t253 * t320 + (-pkin(8) * t277 - pkin(7)) * t292 + t307;
t201 = -qJD(5) * t238 + t248 * t284 + t249 * t288;
t174 = (t238 * t264 - t201) * pkin(9) + (t239 * t264 - t200) * pkin(5) + t191;
t283 = sin(qJ(6));
t287 = cos(qJ(6));
t170 = -t173 * t283 + t174 * t287;
t220 = -t239 * t283 + t264 * t287;
t181 = qJD(6) * t220 + t201 * t287 + t263 * t283;
t221 = t239 * t287 + t264 * t283;
t192 = -mrSges(7,1) * t220 + mrSges(7,2) * t221;
t199 = qJDD(6) - t200;
t234 = qJD(6) + t238;
t202 = -mrSges(7,2) * t234 + mrSges(7,3) * t220;
t166 = m(7) * t170 + mrSges(7,1) * t199 - t181 * mrSges(7,3) - t192 * t221 + t202 * t234;
t171 = t173 * t287 + t174 * t283;
t180 = -qJD(6) * t221 - t201 * t283 + t263 * t287;
t203 = mrSges(7,1) * t234 - mrSges(7,3) * t221;
t167 = m(7) * t171 - mrSges(7,2) * t199 + t180 * mrSges(7,3) + t192 * t220 - t203 * t234;
t313 = -t166 * t283 + t287 * t167;
t153 = m(6) * t178 - mrSges(6,2) * t263 + mrSges(6,3) * t200 - t213 * t238 - t228 * t264 + t313;
t177 = t188 * t288 - t189 * t284;
t227 = -mrSges(6,2) * t264 - mrSges(6,3) * t238;
t172 = -pkin(5) * t263 - pkin(9) * t262 + t214 * t239 - t177;
t304 = -m(7) * t172 + t180 * mrSges(7,1) - t181 * mrSges(7,2) + t220 * t202 - t203 * t221;
t162 = m(6) * t177 + mrSges(6,1) * t263 - mrSges(6,3) * t201 - t213 * t239 + t227 * t264 + t304;
t145 = t284 * t153 + t288 * t162;
t247 = (mrSges(5,1) * t285 + mrSges(5,2) * t289) * qJD(1);
t321 = qJD(1) * t285;
t251 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t321;
t142 = m(5) * t211 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t249 + qJD(4) * t251 - t247 * t320 + t145;
t252 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t320;
t314 = t288 * t153 - t162 * t284;
t143 = m(5) * t212 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t248 - qJD(4) * t252 - t247 * t321 + t314;
t135 = t289 * t142 + t285 * t143;
t155 = t287 * t166 + t283 * t167;
t315 = -t285 * t142 + t289 * t143;
t312 = m(4) * t226 + qJDD(1) * mrSges(4,2) + t135;
t308 = m(6) * t191 - t200 * mrSges(6,1) + t201 * mrSges(6,2) + t238 * t227 + t239 * t228 + t155;
t182 = Ifges(7,5) * t221 + Ifges(7,6) * t220 + Ifges(7,3) * t234;
t184 = Ifges(7,1) * t221 + Ifges(7,4) * t220 + Ifges(7,5) * t234;
t159 = -mrSges(7,1) * t172 + mrSges(7,3) * t171 + Ifges(7,4) * t181 + Ifges(7,2) * t180 + Ifges(7,6) * t199 - t182 * t221 + t184 * t234;
t183 = Ifges(7,4) * t221 + Ifges(7,2) * t220 + Ifges(7,6) * t234;
t160 = mrSges(7,2) * t172 - mrSges(7,3) * t170 + Ifges(7,1) * t181 + Ifges(7,4) * t180 + Ifges(7,5) * t199 + t182 * t220 - t183 * t234;
t204 = Ifges(6,5) * t239 - Ifges(6,6) * t238 + Ifges(6,3) * t264;
t205 = Ifges(6,4) * t239 - Ifges(6,2) * t238 + Ifges(6,6) * t264;
t137 = mrSges(6,2) * t191 - mrSges(6,3) * t177 + Ifges(6,1) * t201 + Ifges(6,4) * t200 + Ifges(6,5) * t263 - pkin(9) * t155 - t159 * t283 + t160 * t287 - t204 * t238 - t205 * t264;
t206 = Ifges(6,1) * t239 - Ifges(6,4) * t238 + Ifges(6,5) * t264;
t299 = mrSges(7,1) * t170 - mrSges(7,2) * t171 + Ifges(7,5) * t181 + Ifges(7,6) * t180 + Ifges(7,3) * t199 + t183 * t221 - t184 * t220;
t138 = -mrSges(6,1) * t191 + mrSges(6,3) * t178 + Ifges(6,4) * t201 + Ifges(6,2) * t200 + Ifges(6,6) * t263 - pkin(5) * t155 - t204 * t239 + t206 * t264 - t299;
t218 = -t292 * pkin(7) + t307;
t235 = (Ifges(5,3) * qJD(4)) + (Ifges(5,5) * t289 - Ifges(5,6) * t285) * qJD(1);
t237 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t289 - Ifges(5,4) * t285) * qJD(1);
t125 = -mrSges(5,1) * t218 + mrSges(5,3) * t212 + Ifges(5,4) * t249 + Ifges(5,2) * t248 + Ifges(5,6) * qJDD(4) - pkin(4) * t308 + pkin(8) * t314 + qJD(4) * t237 + t284 * t137 + t288 * t138 - t235 * t320;
t236 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t289 - Ifges(5,2) * t285) * qJD(1);
t128 = mrSges(5,2) * t218 - mrSges(5,3) * t211 + Ifges(5,1) * t249 + Ifges(5,4) * t248 + Ifges(5,5) * qJDD(4) - pkin(8) * t145 - qJD(4) * t236 + t137 * t288 - t138 * t284 - t235 * t321;
t301 = -m(5) * t218 + t248 * mrSges(5,1) - t249 * mrSges(5,2) - t251 * t321 - t252 * t320 - t308;
t306 = -mrSges(4,1) * t307 + mrSges(4,2) * g(3) + (t292 * Ifges(4,4)) + Ifges(4,5) * qJDD(1) + pkin(3) * t301 + pkin(7) * t315 + t289 * t125 + t285 * t128;
t231 = pkin(1) * t292 - t326;
t305 = -m(3) * t231 + (t292 * mrSges(3,2)) + qJDD(1) * mrSges(3,3) + t312;
t303 = mrSges(4,2) * t226 + mrSges(4,3) * t307 + Ifges(4,1) * qJDD(1) - pkin(7) * t135 - t125 * t285 + t289 * t128;
t302 = mrSges(6,1) * t177 - mrSges(6,2) * t178 + Ifges(6,5) * t201 + Ifges(6,6) * t200 + Ifges(6,3) * t263 + pkin(5) * t304 + pkin(9) * t313 + t287 * t159 + t283 * t160 + t239 * t205 + t238 * t206;
t148 = -m(4) * t307 - t292 * mrSges(4,2) - qJDD(1) * mrSges(4,3) + t301;
t300 = mrSges(3,1) * t233 + pkin(2) * t148 + t306;
t298 = mrSges(3,2) * t233 - mrSges(3,3) * t231 + Ifges(3,1) * qJDD(1) - qJ(3) * t148 + t303;
t297 = mrSges(5,1) * t211 - mrSges(5,2) * t212 + Ifges(5,5) * t249 + Ifges(5,6) * t248 + Ifges(5,3) * qJDD(4) + pkin(4) * t145 + t236 * t320 + t237 * t321 + t302;
t296 = -m(3) * t233 + t292 * mrSges(3,3) - t148;
t295 = -mrSges(2,2) * t255 + qJ(2) * (t305 - t322) + pkin(1) * (-qJDD(1) * mrSges(3,2) + t296) + mrSges(2,1) * t254 + Ifges(2,3) * qJDD(1) + t298;
t294 = mrSges(4,1) * t226 - Ifges(4,4) * qJDD(1) + pkin(3) * t135 + t297;
t293 = mrSges(3,1) * t231 + pkin(2) * (-t312 + t322) + qJ(3) * (-m(4) * g(3) + t315) - t294;
t146 = t296 + (mrSges(2,1) - mrSges(3,2)) * qJDD(1) + m(2) * t254 - mrSges(2,2) * t292;
t132 = (-m(3) - m(4)) * g(3) + t315;
t129 = m(2) * t255 - qJDD(1) * mrSges(2,2) + (-mrSges(2,1) - mrSges(4,3)) * t292 + t305;
t123 = -t293 + ((Ifges(2,5) - t324) * t292) + t323 * qJDD(1) + (mrSges(2,1) - t325) * g(3) + mrSges(2,3) * t255 - pkin(1) * t132;
t122 = t300 + (-Ifges(3,4) + Ifges(2,5)) * qJDD(1) - t323 * t292 + (-mrSges(2,2) + mrSges(3,3)) * g(3) - mrSges(2,3) * t254 - qJ(2) * t132;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t290 * t122 - t286 * t123 - pkin(6) * (t129 * t286 + t146 * t290), t122, t298, t303, t128, t137, t160; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t286 * t122 + t290 * t123 + pkin(6) * (t129 * t290 - t146 * t286), t123, -mrSges(3,3) * g(3) + Ifges(3,4) * qJDD(1) - (t292 * Ifges(3,5)) - t300, -mrSges(4,3) * g(3) - t292 * Ifges(4,5) - t294, t125, t138, t159; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t295, t295, Ifges(3,5) * qJDD(1) + g(3) * t325 + (t324 * t292) + t293, t306, t297, t302, t299;];
m_new  = t1;
