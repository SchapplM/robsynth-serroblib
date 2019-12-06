% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5PPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d4,d5,theta1,theta2,theta3]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 14:59
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5PPPRR2_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR2_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPPRR2_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPPRR2_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPPRR2_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR2_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPPRR2_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPPRR2_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPPRR2_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 14:59:32
% EndTime: 2019-12-05 14:59:33
% DurationCPUTime: 0.95s
% Computational Cost: add. (10831->151), mult. (16416->200), div. (0->0), fcn. (11772->10), ass. (0->71)
t258 = sin(pkin(7));
t261 = cos(pkin(7));
t252 = -g(1) * t261 - g(2) * t258;
t255 = -g(3) + qJDD(1);
t257 = sin(pkin(8));
t260 = cos(pkin(8));
t240 = t252 * t260 + t255 * t257;
t251 = g(1) * t258 - g(2) * t261;
t250 = qJDD(2) - t251;
t256 = sin(pkin(9));
t259 = cos(pkin(9));
t236 = t240 * t259 + t250 * t256;
t239 = -t252 * t257 + t255 * t260;
t238 = qJDD(3) - t239;
t263 = sin(qJ(4));
t265 = cos(qJ(4));
t233 = t265 * t236 + t263 * t238;
t266 = qJD(4) ^ 2;
t231 = -pkin(4) * t266 + qJDD(4) * pkin(6) + t233;
t235 = t240 * t256 - t259 * t250;
t262 = sin(qJ(5));
t264 = cos(qJ(5));
t228 = -t231 * t262 + t235 * t264;
t247 = (-mrSges(6,1) * t264 + mrSges(6,2) * t262) * qJD(4);
t275 = qJD(4) * qJD(5);
t248 = qJDD(4) * t262 + t264 * t275;
t276 = qJD(4) * t264;
t254 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t276;
t277 = qJD(4) * t262;
t226 = m(6) * t228 + qJDD(5) * mrSges(6,1) - mrSges(6,3) * t248 + qJD(5) * t254 - t247 * t277;
t229 = t231 * t264 + t235 * t262;
t249 = qJDD(4) * t264 - t262 * t275;
t253 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t277;
t227 = m(6) * t229 - qJDD(5) * mrSges(6,2) + mrSges(6,3) * t249 - qJD(5) * t253 + t247 * t276;
t270 = -t226 * t262 + t264 * t227;
t220 = m(5) * t233 - mrSges(5,1) * t266 - qJDD(4) * mrSges(5,2) + t270;
t232 = -t236 * t263 + t238 * t265;
t230 = -qJDD(4) * pkin(4) - pkin(6) * t266 - t232;
t267 = -m(6) * t230 + t249 * mrSges(6,1) - mrSges(6,2) * t248 - t253 * t277 + t254 * t276;
t224 = m(5) * t232 + qJDD(4) * mrSges(5,1) - mrSges(5,2) * t266 + t267;
t271 = t265 * t220 - t224 * t263;
t216 = m(4) * t236 + t271;
t221 = t264 * t226 + t262 * t227;
t218 = (-m(4) - m(5)) * t235 - t221;
t272 = t259 * t216 - t218 * t256;
t209 = m(3) * t240 + t272;
t217 = t220 * t263 + t224 * t265;
t268 = -m(4) * t238 - t217;
t214 = m(3) * t239 + t268;
t273 = t260 * t209 - t214 * t257;
t203 = m(2) * t252 + t273;
t210 = t216 * t256 + t218 * t259;
t269 = -m(3) * t250 - t210;
t208 = m(2) * t251 + t269;
t278 = t258 * t203 + t261 * t208;
t204 = t257 * t209 + t260 * t214;
t274 = t261 * t203 - t208 * t258;
t243 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t262 + Ifges(6,4) * t264) * qJD(4);
t242 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t262 + Ifges(6,2) * t264) * qJD(4);
t241 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t262 + Ifges(6,6) * t264) * qJD(4);
t223 = mrSges(6,2) * t230 - mrSges(6,3) * t228 + Ifges(6,1) * t248 + Ifges(6,4) * t249 + Ifges(6,5) * qJDD(5) - qJD(5) * t242 + t241 * t276;
t222 = -mrSges(6,1) * t230 + mrSges(6,3) * t229 + Ifges(6,4) * t248 + Ifges(6,2) * t249 + Ifges(6,6) * qJDD(5) + qJD(5) * t243 - t241 * t277;
t212 = -mrSges(5,1) * t235 - mrSges(6,1) * t228 + mrSges(6,2) * t229 + mrSges(5,3) * t233 + t266 * Ifges(5,5) - Ifges(6,5) * t248 + Ifges(5,6) * qJDD(4) - Ifges(6,6) * t249 - Ifges(6,3) * qJDD(5) - pkin(4) * t221 + (-t242 * t262 + t243 * t264) * qJD(4);
t211 = mrSges(5,2) * t235 - mrSges(5,3) * t232 + Ifges(5,5) * qJDD(4) - Ifges(5,6) * t266 - pkin(6) * t221 - t222 * t262 + t223 * t264;
t200 = -mrSges(4,1) * t238 - mrSges(5,1) * t232 + mrSges(5,2) * t233 + mrSges(4,3) * t236 - Ifges(5,3) * qJDD(4) - pkin(3) * t217 - pkin(4) * t267 - pkin(6) * t270 - t264 * t222 - t262 * t223;
t199 = mrSges(4,2) * t238 + mrSges(4,3) * t235 - pkin(5) * t217 + t211 * t265 - t212 * t263;
t198 = -mrSges(3,1) * t250 + mrSges(3,3) * t240 + mrSges(4,1) * t235 + mrSges(4,2) * t236 - t263 * t211 - t265 * t212 - pkin(3) * (-m(5) * t235 - t221) - pkin(5) * t271 - pkin(2) * t210;
t197 = mrSges(3,2) * t250 - mrSges(3,3) * t239 - qJ(3) * t210 + t199 * t259 - t200 * t256;
t196 = -mrSges(2,1) * t255 - mrSges(3,1) * t239 + mrSges(3,2) * t240 + mrSges(2,3) * t252 - pkin(1) * t204 - pkin(2) * t268 - qJ(3) * t272 - t256 * t199 - t259 * t200;
t195 = mrSges(2,2) * t255 - mrSges(2,3) * t251 - qJ(2) * t204 + t197 * t260 - t198 * t257;
t1 = [-m(1) * g(1) + t274; -m(1) * g(2) + t278; -m(1) * g(3) + m(2) * t255 + t204; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t278 + t261 * t195 - t258 * t196; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t274 + t258 * t195 + t261 * t196; -mrSges(1,1) * g(2) + mrSges(2,1) * t251 + mrSges(1,2) * g(1) - mrSges(2,2) * t252 + pkin(1) * t269 + qJ(2) * t273 + t257 * t197 + t260 * t198;];
tauB = t1;
