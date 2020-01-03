% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S4RRPR4
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
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
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
% tauJ [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S4RRPR4_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR4_invdynJ_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR4_invdynJ_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR4_invdynJ_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR4_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR4_invdynJ_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR4_invdynJ_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR4_invdynJ_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR4_invdynJ_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:02:30
% EndTime: 2019-12-31 17:02:30
% DurationCPUTime: 0.38s
% Computational Cost: add. (2672->117), mult. (3695->157), div. (0->0), fcn. (2154->8), ass. (0->60)
t249 = qJD(1) + qJD(2);
t245 = t249 ^ 2;
t251 = cos(pkin(7));
t276 = pkin(3) * t251;
t250 = sin(pkin(7));
t275 = mrSges(4,2) * t250;
t248 = t251 ^ 2;
t273 = t245 * t248;
t246 = qJDD(1) + qJDD(2);
t272 = t246 * t251;
t254 = sin(qJ(1));
t257 = cos(qJ(1));
t269 = t254 * g(1) - g(2) * t257;
t236 = qJDD(1) * pkin(1) + t269;
t265 = -g(1) * t257 - t254 * g(2);
t237 = -qJD(1) ^ 2 * pkin(1) + t265;
t253 = sin(qJ(2));
t256 = cos(qJ(2));
t226 = t253 * t236 + t256 * t237;
t223 = -pkin(2) * t245 + qJ(3) * t246 + t226;
t270 = qJD(3) * t249;
t268 = -t251 * g(3) - 0.2e1 * t250 * t270;
t208 = (-pkin(6) * t246 + t245 * t276 - t223) * t250 + t268;
t212 = -g(3) * t250 + (t223 + 0.2e1 * t270) * t251;
t209 = -pkin(3) * t273 + pkin(6) * t272 + t212;
t252 = sin(qJ(4));
t255 = cos(qJ(4));
t206 = t208 * t255 - t209 * t252;
t261 = -t250 * t252 + t251 * t255;
t229 = t261 * t249;
t262 = t250 * t255 + t251 * t252;
t230 = t262 * t249;
t219 = -mrSges(5,1) * t229 + mrSges(5,2) * t230;
t222 = t229 * qJD(4) + t262 * t246;
t227 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t229;
t204 = m(5) * t206 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t222 + qJD(4) * t227 - t219 * t230;
t207 = t208 * t252 + t209 * t255;
t221 = -t230 * qJD(4) + t261 * t246;
t228 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t230;
t205 = m(5) * t207 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t221 - qJD(4) * t228 + t219 * t229;
t271 = t255 * t204 + t252 * t205;
t211 = -t250 * t223 + t268;
t263 = mrSges(4,3) * t246 + (-mrSges(4,1) * t251 + t275) * t245;
t266 = -t252 * t204 + t255 * t205;
t267 = -(m(4) * t211 - t263 * t250 + t271) * t250 + t251 * (m(4) * t212 + t263 * t251 + t266);
t225 = t256 * t236 - t253 * t237;
t264 = qJDD(3) - t225;
t247 = t250 ^ 2;
t210 = (-pkin(2) - t276) * t246 + (-qJ(3) + (-t247 - t248) * pkin(6)) * t245 + t264;
t213 = Ifges(5,5) * t230 + Ifges(5,6) * t229 + Ifges(5,3) * qJD(4);
t215 = Ifges(5,1) * t230 + Ifges(5,4) * t229 + Ifges(5,5) * qJD(4);
t197 = -mrSges(5,1) * t210 + mrSges(5,3) * t207 + Ifges(5,4) * t222 + Ifges(5,2) * t221 + Ifges(5,6) * qJDD(4) + qJD(4) * t215 - t213 * t230;
t214 = Ifges(5,4) * t230 + Ifges(5,2) * t229 + Ifges(5,6) * qJD(4);
t198 = mrSges(5,2) * t210 - mrSges(5,3) * t206 + Ifges(5,1) * t222 + Ifges(5,4) * t221 + Ifges(5,5) * qJDD(4) - qJD(4) * t214 + t213 * t229;
t220 = -t246 * pkin(2) - t245 * qJ(3) + t264;
t259 = m(5) * t210 - t221 * mrSges(5,1) + mrSges(5,2) * t222 - t229 * t227 + t228 * t230;
t258 = -m(4) * t220 + mrSges(4,1) * t272 - t259 + (t245 * t247 + t273) * mrSges(4,3);
t200 = t246 * t275 - t258;
t260 = -mrSges(3,2) * t226 + t251 * (-mrSges(4,1) * t220 + mrSges(4,3) * t212 + t252 * t198 + t255 * t197 - pkin(3) * t259 + pkin(6) * t266 + (Ifges(4,4) * t250 + Ifges(4,2) * t251) * t246) + t250 * (mrSges(4,2) * t220 - mrSges(4,3) * t211 + t255 * t198 - t252 * t197 - pkin(6) * t271 + (Ifges(4,1) * t250 + Ifges(4,4) * t251) * t246) + qJ(3) * t267 - pkin(2) * t200 + mrSges(3,1) * t225 + Ifges(3,3) * t246;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t269 - mrSges(2,2) * t265 + pkin(1) * (t253 * (m(3) * t226 - mrSges(3,1) * t245 - mrSges(3,2) * t246 + t267) + t256 * (m(3) * t225 - mrSges(3,2) * t245 + (mrSges(3,1) - t275) * t246 + t258)) + t260; t260; t200; mrSges(5,1) * t206 - mrSges(5,2) * t207 + Ifges(5,5) * t222 + Ifges(5,6) * t221 + Ifges(5,3) * qJDD(4) + t214 * t230 - t215 * t229;];
tauJ = t1;
