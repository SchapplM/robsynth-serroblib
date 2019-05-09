% Calculate vector of cutting torques with Newton-Euler for
% S6RPRPRP8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
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
% Datum: 2019-05-05 18:03
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPRPRP8_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP8_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP8_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP8_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP8_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP8_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP8_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP8_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP8_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 18:00:31
% EndTime: 2019-05-05 18:00:44
% DurationCPUTime: 6.12s
% Computational Cost: add. (95505->338), mult. (206453->402), div. (0->0), fcn. (131950->8), ass. (0->126)
t347 = -2 * qJD(4);
t303 = sin(qJ(1));
t305 = cos(qJ(1));
t282 = t303 * g(1) - t305 * g(2);
t307 = qJD(1) ^ 2;
t321 = -t307 * qJ(2) + qJDD(2) - t282;
t345 = -pkin(1) - pkin(7);
t254 = t345 * qJDD(1) + t321;
t302 = sin(qJ(3));
t304 = cos(qJ(3));
t244 = t302 * g(3) + t304 * t254;
t333 = qJD(1) * qJD(3);
t330 = t302 * t333;
t277 = qJDD(1) * t304 - t330;
t214 = (-t277 - t330) * qJ(4) + (-t302 * t304 * t307 + qJDD(3)) * pkin(3) + t244;
t245 = -g(3) * t304 + t302 * t254;
t276 = -qJDD(1) * t302 - t304 * t333;
t335 = qJD(1) * t304;
t280 = qJD(3) * pkin(3) - qJ(4) * t335;
t296 = t302 ^ 2;
t215 = -pkin(3) * t296 * t307 + qJ(4) * t276 - qJD(3) * t280 + t245;
t299 = sin(pkin(9));
t300 = cos(pkin(9));
t336 = qJD(1) * t302;
t265 = -t299 * t336 + t300 * t335;
t189 = t300 * t214 - t299 * t215 + t265 * t347;
t283 = -t305 * g(1) - t303 * g(2);
t322 = qJDD(1) * qJ(2) + 0.2e1 * qJD(2) * qJD(1) + t283;
t264 = (t299 * t304 + t300 * t302) * qJD(1);
t190 = t299 * t214 + t300 * t215 + t264 * t347;
t232 = pkin(4) * t264 - pkin(8) * t265;
t306 = qJD(3) ^ 2;
t185 = -pkin(4) * t306 + qJDD(3) * pkin(8) - t232 * t264 + t190;
t217 = -t276 * pkin(3) + qJDD(4) + t280 * t335 + (-qJ(4) * t296 + t345) * t307 + t322;
t241 = t276 * t300 - t277 * t299;
t242 = t276 * t299 + t277 * t300;
t187 = (qJD(3) * t264 - t242) * pkin(8) + (qJD(3) * t265 - t241) * pkin(4) + t217;
t301 = sin(qJ(5));
t344 = cos(qJ(5));
t181 = -t301 * t185 + t344 * t187;
t182 = t344 * t185 + t301 * t187;
t246 = -t344 * qJD(3) + t301 * t265;
t247 = t301 * qJD(3) + t344 * t265;
t262 = qJD(5) + t264;
t193 = Ifges(7,5) * t247 + Ifges(7,6) * t262 + Ifges(7,3) * t246;
t196 = Ifges(6,4) * t247 - Ifges(6,2) * t246 + Ifges(6,6) * t262;
t198 = Ifges(6,1) * t247 - Ifges(6,4) * t246 + Ifges(6,5) * t262;
t207 = t247 * qJD(5) - t344 * qJDD(3) + t301 * t242;
t208 = -t246 * qJD(5) + t301 * qJDD(3) + t344 * t242;
t220 = mrSges(7,1) * t246 - mrSges(7,3) * t247;
t240 = qJDD(5) - t241;
t219 = pkin(5) * t246 - qJ(6) * t247;
t261 = t262 ^ 2;
t177 = -pkin(5) * t261 + qJ(6) * t240 + 0.2e1 * qJD(6) * t262 - t219 * t246 + t182;
t179 = -t240 * pkin(5) - t261 * qJ(6) + t247 * t219 + qJDD(6) - t181;
t197 = Ifges(7,1) * t247 + Ifges(7,4) * t262 + Ifges(7,5) * t246;
t320 = mrSges(7,1) * t179 - mrSges(7,3) * t177 - Ifges(7,4) * t208 - Ifges(7,2) * t240 - Ifges(7,6) * t207 - t246 * t197;
t223 = -mrSges(7,2) * t246 + mrSges(7,3) * t262;
t326 = -m(7) * t179 + t240 * mrSges(7,1) + t262 * t223;
t226 = -mrSges(7,1) * t262 + mrSges(7,2) * t247;
t331 = m(7) * t177 + t240 * mrSges(7,3) + t262 * t226;
t346 = -(-t196 + t193) * t247 + mrSges(6,1) * t181 - mrSges(6,2) * t182 + Ifges(6,5) * t208 - Ifges(6,6) * t207 + Ifges(6,3) * t240 + pkin(5) * (-t208 * mrSges(7,2) - t247 * t220 + t326) + qJ(6) * (-t207 * mrSges(7,2) - t246 * t220 + t331) + t246 * t198 - t320;
t343 = mrSges(2,1) - mrSges(3,2);
t342 = -mrSges(6,3) - mrSges(7,2);
t341 = Ifges(2,5) - Ifges(3,4);
t340 = -Ifges(2,6) + Ifges(3,5);
t231 = mrSges(5,1) * t264 + mrSges(5,2) * t265;
t253 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t265;
t225 = mrSges(6,1) * t262 - mrSges(6,3) * t247;
t337 = -mrSges(6,1) * t246 - mrSges(6,2) * t247 - t220;
t167 = m(6) * t182 - t240 * mrSges(6,2) + t342 * t207 - t262 * t225 + t337 * t246 + t331;
t224 = -mrSges(6,2) * t262 - mrSges(6,3) * t246;
t169 = m(6) * t181 + t240 * mrSges(6,1) + t342 * t208 + t262 * t224 + t337 * t247 + t326;
t328 = t344 * t167 - t169 * t301;
t155 = m(5) * t190 - qJDD(3) * mrSges(5,2) + mrSges(5,3) * t241 - qJD(3) * t253 - t231 * t264 + t328;
t252 = -qJD(3) * mrSges(5,2) - mrSges(5,3) * t264;
t184 = -qJDD(3) * pkin(4) - t306 * pkin(8) + t265 * t232 - t189;
t180 = -0.2e1 * qJD(6) * t247 + (t246 * t262 - t208) * qJ(6) + (t247 * t262 + t207) * pkin(5) + t184;
t174 = m(7) * t180 + t207 * mrSges(7,1) - t208 * mrSges(7,3) + t223 * t246 - t247 * t226;
t312 = -m(6) * t184 - t207 * mrSges(6,1) - mrSges(6,2) * t208 - t246 * t224 - t225 * t247 - t174;
t164 = m(5) * t189 + qJDD(3) * mrSges(5,1) - mrSges(5,3) * t242 + qJD(3) * t252 - t231 * t265 + t312;
t149 = t299 * t155 + t300 * t164;
t161 = t301 * t167 + t344 * t169;
t195 = Ifges(7,4) * t247 + Ifges(7,2) * t262 + Ifges(7,6) * t246;
t339 = -Ifges(6,5) * t247 + Ifges(6,6) * t246 - Ifges(6,3) * t262 - t195;
t275 = (mrSges(4,1) * t302 + mrSges(4,2) * t304) * qJD(1);
t279 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t336;
t146 = m(4) * t244 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t277 + qJD(3) * t279 - t275 * t335 + t149;
t281 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t335;
t329 = t300 * t155 - t164 * t299;
t147 = m(4) * t245 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t276 - qJD(3) * t281 - t275 * t336 + t329;
t142 = -t146 * t302 + t304 * t147;
t325 = -mrSges(7,1) * t180 + mrSges(7,2) * t177;
t141 = t304 * t146 + t302 * t147;
t319 = mrSges(7,2) * t179 - mrSges(7,3) * t180 + Ifges(7,1) * t208 + Ifges(7,4) * t240 + Ifges(7,5) * t207 + t262 * t193;
t263 = -qJDD(1) * pkin(1) + t321;
t318 = -m(3) * t263 + t307 * mrSges(3,3) - t141;
t317 = m(5) * t217 - mrSges(5,1) * t241 + t242 * mrSges(5,2) + t252 * t264 + t265 * t253 + t161;
t157 = -mrSges(6,1) * t184 + mrSges(6,3) * t182 - pkin(5) * t174 + (t197 + t198) * t262 + t339 * t247 + (Ifges(6,6) - Ifges(7,6)) * t240 + (Ifges(6,4) - Ifges(7,5)) * t208 + (-Ifges(6,2) - Ifges(7,3)) * t207 + t325;
t159 = mrSges(6,2) * t184 - mrSges(6,3) * t181 + Ifges(6,1) * t208 - Ifges(6,4) * t207 + Ifges(6,5) * t240 - qJ(6) * t174 - t262 * t196 + t339 * t246 + t319;
t227 = Ifges(5,5) * t265 - Ifges(5,6) * t264 + Ifges(5,3) * qJD(3);
t228 = Ifges(5,4) * t265 - Ifges(5,2) * t264 + Ifges(5,6) * qJD(3);
t143 = mrSges(5,2) * t217 - mrSges(5,3) * t189 + Ifges(5,1) * t242 + Ifges(5,4) * t241 + Ifges(5,5) * qJDD(3) - pkin(8) * t161 - qJD(3) * t228 - t301 * t157 + t344 * t159 - t264 * t227;
t229 = Ifges(5,1) * t265 - Ifges(5,4) * t264 + Ifges(5,5) * qJD(3);
t144 = -mrSges(5,1) * t217 + mrSges(5,3) * t190 + Ifges(5,4) * t242 + Ifges(5,2) * t241 + Ifges(5,6) * qJDD(3) - pkin(4) * t161 + qJD(3) * t229 - t265 * t227 - t346;
t251 = t345 * t307 + t322;
t266 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t304 - Ifges(4,6) * t302) * qJD(1);
t268 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t304 - Ifges(4,4) * t302) * qJD(1);
t135 = -mrSges(4,1) * t251 + mrSges(4,3) * t245 + Ifges(4,4) * t277 + Ifges(4,2) * t276 + Ifges(4,6) * qJDD(3) - pkin(3) * t317 + qJ(4) * t329 + qJD(3) * t268 + t299 * t143 + t300 * t144 - t266 * t335;
t267 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t304 - Ifges(4,2) * t302) * qJD(1);
t137 = mrSges(4,2) * t251 - mrSges(4,3) * t244 + Ifges(4,1) * t277 + Ifges(4,4) * t276 + Ifges(4,5) * qJDD(3) - qJ(4) * t149 - qJD(3) * t267 + t143 * t300 - t144 * t299 - t266 * t336;
t257 = t307 * pkin(1) - t322;
t316 = mrSges(3,2) * t263 - mrSges(3,3) * t257 + Ifges(3,1) * qJDD(1) - pkin(7) * t141 - t135 * t302 + t304 * t137;
t152 = -m(4) * t251 + mrSges(4,1) * t276 - t277 * mrSges(4,2) - t279 * t336 - t281 * t335 - t317;
t315 = -mrSges(3,1) * t257 - pkin(2) * t152 - pkin(7) * t142 - t304 * t135 - t302 * t137;
t314 = -mrSges(5,1) * t189 + mrSges(5,2) * t190 - Ifges(5,5) * t242 - Ifges(5,6) * t241 - Ifges(5,3) * qJDD(3) - pkin(4) * t312 - pkin(8) * t328 - t344 * t157 - t301 * t159 - t265 * t228 - t264 * t229;
t311 = -m(3) * t257 + t307 * mrSges(3,2) + qJDD(1) * mrSges(3,3) - t152;
t313 = -mrSges(2,2) * t283 + pkin(1) * (-qJDD(1) * mrSges(3,2) + t318) + qJ(2) * t311 + mrSges(2,1) * t282 + Ifges(2,3) * qJDD(1) + t316;
t310 = -mrSges(4,1) * t244 + mrSges(4,2) * t245 - Ifges(4,5) * t277 - Ifges(4,6) * t276 - Ifges(4,3) * qJDD(3) - pkin(3) * t149 - t267 * t335 - t268 * t336 + t314;
t308 = -mrSges(3,1) * t263 - pkin(2) * t141 + t310;
t150 = m(2) * t283 - mrSges(2,1) * t307 - qJDD(1) * mrSges(2,2) + t311;
t140 = -m(3) * g(3) + t142;
t138 = m(2) * t282 - t307 * mrSges(2,2) + t343 * qJDD(1) + t318;
t134 = -t308 + t340 * t307 + t341 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3) - mrSges(2,3) * t282 - qJ(2) * t140;
t133 = mrSges(2,3) * t283 - pkin(1) * t140 + t343 * g(3) - t340 * qJDD(1) + t341 * t307 + t315;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t305 * t134 - t303 * t133 - pkin(6) * (t138 * t305 + t150 * t303), t134, t316, t137, t143, t159, -t195 * t246 + t319; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t303 * t134 + t305 * t133 + pkin(6) * (-t138 * t303 + t150 * t305), t133, -mrSges(3,3) * g(3) + Ifges(3,4) * qJDD(1) - t307 * Ifges(3,5) + t308, t135, t144, t157, -t247 * t193 - t320; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t313, t313, mrSges(3,2) * g(3) + t307 * Ifges(3,4) + Ifges(3,5) * qJDD(1) - t315, -t310, -t314, t346, Ifges(7,5) * t208 + Ifges(7,6) * t240 + Ifges(7,3) * t207 + t247 * t195 - t262 * t197 - t325;];
m_new  = t1;
