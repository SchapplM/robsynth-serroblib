% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5PPRRP4
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1]';
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
% tauJ [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5PPRRP4_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP4_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP4_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRP4_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRP4_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRRP4_invdynJ_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRP4_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRP4_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRP4_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:34:29
% EndTime: 2019-12-31 17:34:30
% DurationCPUTime: 0.46s
% Computational Cost: add. (639->125), mult. (1173->150), div. (0->0), fcn. (576->6), ass. (0->57)
t252 = cos(qJ(4));
t269 = Ifges(5,4) + Ifges(6,4);
t276 = t252 * t269;
t275 = Ifges(5,1) + Ifges(6,1);
t268 = Ifges(5,5) + Ifges(6,5);
t274 = Ifges(5,2) + Ifges(6,2);
t273 = Ifges(5,6) + Ifges(6,6);
t250 = sin(qJ(4));
t272 = -t268 * qJD(4) + (-t275 * t250 - t276) * qJD(3);
t254 = qJD(3) ^ 2;
t271 = pkin(4) * t254;
t270 = -mrSges(5,2) - mrSges(6,2);
t248 = sin(pkin(7));
t249 = cos(pkin(7));
t236 = -g(1) * t248 + g(2) * t249 + qJDD(2);
t237 = -g(1) * t249 - g(2) * t248;
t251 = sin(qJ(3));
t253 = cos(qJ(3));
t215 = t251 * t236 + t253 * t237;
t213 = -pkin(3) * t254 + qJDD(3) * pkin(6) + t215;
t247 = g(3) - qJDD(1);
t210 = t252 * t213 + t250 * t247;
t267 = t273 * qJD(4) + (t269 * t250 + t252 * t274) * qJD(3);
t265 = qJD(3) * t250;
t239 = qJD(4) * mrSges(6,1) - mrSges(6,3) * t265;
t266 = -qJD(4) * mrSges(5,1) + mrSges(5,3) * t265 - t239;
t264 = qJD(3) * t252;
t263 = qJD(3) * qJD(4);
t262 = qJD(3) * qJD(5);
t259 = t252 * t263;
t233 = qJDD(3) * t250 + t259;
t244 = t252 * t247;
t206 = qJDD(4) * pkin(4) + t244 + (-t233 + t259) * qJ(5) + (t252 * t271 - t213 - 0.2e1 * t262) * t250;
t241 = -qJD(4) * mrSges(6,2) + mrSges(6,3) * t264;
t261 = m(6) * t206 + qJDD(4) * mrSges(6,1) + qJD(4) * t241;
t234 = qJDD(3) * t252 - t250 * t263;
t238 = qJD(4) * pkin(4) - qJ(5) * t265;
t246 = t252 ^ 2;
t207 = qJ(5) * t234 - qJD(4) * t238 - t246 * t271 + 0.2e1 * t252 * t262 + t210;
t231 = (-t252 * mrSges(6,1) + mrSges(6,2) * t250) * qJD(3);
t260 = m(6) * t207 + t234 * mrSges(6,3) + t231 * t264;
t214 = t236 * t253 - t251 * t237;
t258 = qJD(3) * t266;
t256 = -qJDD(3) * pkin(3) - t214;
t208 = t238 * t265 - pkin(4) * t234 + qJDD(5) + (-qJ(5) * t246 - pkin(6)) * t254 + t256;
t257 = m(6) * t208 - t234 * mrSges(6,1) - t241 * t264;
t212 = -pkin(6) * t254 + t256;
t242 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t264;
t255 = -m(5) * t212 + t234 * mrSges(5,1) + t242 * t264 - t257;
t232 = (-t252 * mrSges(5,1) + mrSges(5,2) * t250) * qJD(3);
t209 = -t213 * t250 + t244;
t203 = mrSges(6,2) * t233 + t239 * t265 + t257;
t202 = -mrSges(6,3) * t233 - t231 * t265 + t261;
t201 = m(5) * t210 + mrSges(5,3) * t234 + t266 * qJD(4) + t270 * qJDD(4) + t232 * t264 + t260;
t200 = m(5) * t209 + qJDD(4) * mrSges(5,1) + qJD(4) * t242 + (-mrSges(5,3) - mrSges(6,3)) * t233 + (-t231 - t232) * t265 + t261;
t199 = t252 * t201;
t1 = [-t200 * t252 - t201 * t250 + (-m(2) - m(3) - m(4)) * t247; m(3) * t236 + t251 * (m(4) * t215 - mrSges(4,1) * t254 - qJDD(3) * mrSges(4,2) - t200 * t250 + t199) + t253 * (m(4) * t214 + qJDD(3) * mrSges(4,1) - mrSges(4,2) * t254 + t270 * t233 + t250 * t258 + t255); Ifges(4,3) * qJDD(3) + mrSges(4,1) * t214 - mrSges(4,2) * t215 + t252 * (-mrSges(5,1) * t212 - mrSges(6,1) * t208 + mrSges(5,3) * t210 + mrSges(6,3) * t207 - pkin(4) * t203 + qJ(5) * t260 + t274 * t234 + (-qJ(5) * mrSges(6,2) + t273) * qJDD(4) + (-qJ(5) * t239 - t272) * qJD(4)) + pkin(3) * t255 + pkin(6) * t199 + (pkin(3) * t270 + t276) * t233 + (mrSges(5,2) * t212 + mrSges(6,2) * t208 - mrSges(5,3) * t209 - mrSges(6,3) * t206 + pkin(3) * t258 - pkin(6) * t200 - qJ(5) * t202 - t267 * qJD(4) + t268 * qJDD(4) + t275 * t233 + t269 * t234) * t250; mrSges(5,1) * t209 + mrSges(6,1) * t206 - mrSges(5,2) * t210 - mrSges(6,2) * t207 + pkin(4) * t202 + t273 * t234 + t268 * t233 + (Ifges(5,3) + Ifges(6,3)) * qJDD(4) + (t267 * t250 + t272 * t252) * qJD(3); t203;];
tauJ = t1;
