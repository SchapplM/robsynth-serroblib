% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S4PRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S4PRPR5_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR5_invdynB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR5_invdynB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPR5_invdynB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR5_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR5_invdynB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR5_invdynB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPR5_invdynB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRPR5_invdynB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:23:01
% EndTime: 2019-12-31 16:23:03
% DurationCPUTime: 0.59s
% Computational Cost: add. (4809->142), mult. (7850->185), div. (0->0), fcn. (4448->8), ass. (0->61)
t259 = sin(pkin(6));
t261 = cos(pkin(6));
t251 = -t261 * g(1) - t259 * g(2);
t257 = -g(3) + qJDD(1);
t263 = sin(qJ(2));
t265 = cos(qJ(2));
t238 = -t263 * t251 + t265 * t257;
t236 = qJDD(2) * pkin(2) + t238;
t239 = t265 * t251 + t263 * t257;
t266 = qJD(2) ^ 2;
t237 = -t266 * pkin(2) + t239;
t258 = sin(pkin(7));
t260 = cos(pkin(7));
t233 = t258 * t236 + t260 * t237;
t231 = -t266 * pkin(3) + qJDD(2) * pkin(5) + t233;
t250 = t259 * g(1) - t261 * g(2);
t249 = qJDD(3) - t250;
t262 = sin(qJ(4));
t264 = cos(qJ(4));
t228 = -t262 * t231 + t264 * t249;
t246 = (-mrSges(5,1) * t264 + mrSges(5,2) * t262) * qJD(2);
t273 = qJD(2) * qJD(4);
t247 = t262 * qJDD(2) + t264 * t273;
t274 = qJD(2) * t264;
t253 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t274;
t275 = qJD(2) * t262;
t226 = m(5) * t228 + qJDD(4) * mrSges(5,1) - t247 * mrSges(5,3) + qJD(4) * t253 - t246 * t275;
t229 = t264 * t231 + t262 * t249;
t248 = t264 * qJDD(2) - t262 * t273;
t252 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t275;
t227 = m(5) * t229 - qJDD(4) * mrSges(5,2) + t248 * mrSges(5,3) - qJD(4) * t252 + t246 * t274;
t269 = -t262 * t226 + t264 * t227;
t215 = m(4) * t233 - t266 * mrSges(4,1) - qJDD(2) * mrSges(4,2) + t269;
t232 = t260 * t236 - t258 * t237;
t230 = -qJDD(2) * pkin(3) - t266 * pkin(5) - t232;
t267 = -m(5) * t230 + t248 * mrSges(5,1) - t247 * mrSges(5,2) - t252 * t275 + t253 * t274;
t222 = m(4) * t232 + qJDD(2) * mrSges(4,1) - t266 * mrSges(4,2) + t267;
t212 = t258 * t215 + t260 * t222;
t210 = m(3) * t238 + qJDD(2) * mrSges(3,1) - t266 * mrSges(3,2) + t212;
t270 = t260 * t215 - t258 * t222;
t211 = m(3) * t239 - t266 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t270;
t271 = -t263 * t210 + t265 * t211;
t203 = m(2) * t251 + t271;
t218 = t264 * t226 + t262 * t227;
t268 = m(4) * t249 + t218;
t217 = (m(2) + m(3)) * t250 - t268;
t276 = t259 * t203 + t261 * t217;
t204 = t265 * t210 + t263 * t211;
t272 = t261 * t203 - t259 * t217;
t242 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t262 + Ifges(5,4) * t264) * qJD(2);
t241 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t262 + Ifges(5,2) * t264) * qJD(2);
t240 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t262 + Ifges(5,6) * t264) * qJD(2);
t220 = mrSges(5,2) * t230 - mrSges(5,3) * t228 + Ifges(5,1) * t247 + Ifges(5,4) * t248 + Ifges(5,5) * qJDD(4) - qJD(4) * t241 + t240 * t274;
t219 = -mrSges(5,1) * t230 + mrSges(5,3) * t229 + Ifges(5,4) * t247 + Ifges(5,2) * t248 + Ifges(5,6) * qJDD(4) + qJD(4) * t242 - t240 * t275;
t206 = -mrSges(4,1) * t249 - mrSges(5,1) * t228 + mrSges(5,2) * t229 + mrSges(4,3) * t233 + t266 * Ifges(4,5) - Ifges(5,5) * t247 + Ifges(4,6) * qJDD(2) - Ifges(5,6) * t248 - Ifges(5,3) * qJDD(4) - pkin(3) * t218 + (-t241 * t262 + t242 * t264) * qJD(2);
t205 = mrSges(4,2) * t249 - mrSges(4,3) * t232 + Ifges(4,5) * qJDD(2) - t266 * Ifges(4,6) - pkin(5) * t218 - t262 * t219 + t264 * t220;
t200 = -mrSges(3,2) * t250 - mrSges(3,3) * t238 + Ifges(3,5) * qJDD(2) - t266 * Ifges(3,6) - qJ(3) * t212 + t260 * t205 - t258 * t206;
t199 = mrSges(3,1) * t250 + mrSges(3,3) * t239 + t266 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t268 + qJ(3) * t270 + t258 * t205 + t260 * t206;
t198 = -pkin(1) * t204 + mrSges(2,3) * t251 - pkin(2) * t212 - mrSges(3,1) * t238 + mrSges(3,2) * t239 - pkin(5) * t269 - mrSges(4,1) * t232 + mrSges(4,2) * t233 - t262 * t220 - t264 * t219 - pkin(3) * t267 - mrSges(2,1) * t257 + (-Ifges(3,3) - Ifges(4,3)) * qJDD(2);
t197 = mrSges(2,2) * t257 - mrSges(2,3) * t250 - pkin(4) * t204 - t263 * t199 + t265 * t200;
t1 = [-m(1) * g(1) + t272; -m(1) * g(2) + t276; -m(1) * g(3) + m(2) * t257 + t204; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t276 + t261 * t197 - t259 * t198; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t272 + t259 * t197 + t261 * t198; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + mrSges(2,1) * t250 - mrSges(2,2) * t251 + t263 * t200 + t265 * t199 + pkin(1) * (m(3) * t250 - t268) + pkin(4) * t271;];
tauB = t1;
