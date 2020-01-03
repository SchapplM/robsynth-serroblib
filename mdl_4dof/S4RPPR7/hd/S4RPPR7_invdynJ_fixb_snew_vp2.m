% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S4RPPR7
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta3]';
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
% Datum: 2019-12-31 16:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S4RPPR7_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR7_invdynJ_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR7_invdynJ_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPR7_invdynJ_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR7_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR7_invdynJ_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR7_invdynJ_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPR7_invdynJ_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPPR7_invdynJ_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:41:38
% EndTime: 2019-12-31 16:41:38
% DurationCPUTime: 0.39s
% Computational Cost: add. (1026->110), mult. (2170->147), div. (0->0), fcn. (1205->6), ass. (0->55)
t255 = qJD(1) ^ 2;
t252 = sin(qJ(1));
t254 = cos(qJ(1));
t268 = t252 * g(1) - t254 * g(2);
t259 = -t255 * qJ(2) + qJDD(2) - t268;
t276 = -pkin(1) - qJ(3);
t280 = -(2 * qJD(1) * qJD(3)) + t276 * qJDD(1) + t259;
t250 = cos(pkin(6));
t279 = t250 ^ 2;
t265 = -t254 * g(1) - t252 * g(2);
t278 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t265;
t277 = pkin(3) * t255;
t275 = t250 * mrSges(4,2);
t249 = sin(pkin(6));
t222 = t249 * g(3) + t280 * t250;
t212 = (-pkin(5) * qJDD(1) - t249 * t277) * t250 + t222;
t223 = -t250 * g(3) + t280 * t249;
t246 = t249 ^ 2;
t270 = t249 * qJDD(1);
t213 = -pkin(5) * t270 - t246 * t277 + t223;
t251 = sin(qJ(4));
t253 = cos(qJ(4));
t210 = t253 * t212 - t251 * t213;
t263 = -t249 * t253 - t250 * t251;
t235 = t263 * qJD(1);
t262 = -t249 * t251 + t250 * t253;
t236 = t262 * qJD(1);
t220 = -t235 * mrSges(5,1) + t236 * mrSges(5,2);
t225 = t235 * qJD(4) + t262 * qJDD(1);
t230 = -qJD(4) * mrSges(5,2) + t235 * mrSges(5,3);
t208 = m(5) * t210 + qJDD(4) * mrSges(5,1) - t225 * mrSges(5,3) + qJD(4) * t230 - t236 * t220;
t211 = t251 * t212 + t253 * t213;
t224 = -t236 * qJD(4) + t263 * qJDD(1);
t231 = qJD(4) * mrSges(5,1) - t236 * mrSges(5,3);
t209 = m(5) * t211 - qJDD(4) * mrSges(5,2) + t224 * mrSges(5,3) - qJD(4) * t231 + t235 * t220;
t274 = t253 * t208 + t251 * t209;
t273 = -t246 - t279;
t267 = t273 * mrSges(4,3);
t266 = -t251 * t208 + t253 * t209;
t261 = -qJDD(1) * mrSges(4,3) - t255 * (t249 * mrSges(4,1) + t275);
t264 = t250 * (m(4) * t222 + t261 * t250 + t274) + t249 * (m(4) * t223 + t261 * t249 + t266);
t258 = qJDD(3) + t278;
t215 = pkin(3) * t270 + (t273 * pkin(5) + t276) * t255 + t258;
t257 = m(5) * t215 - t224 * mrSges(5,1) + t225 * mrSges(5,2) - t235 * t230 + t236 * t231;
t229 = t276 * t255 + t258;
t256 = m(4) * t229 + mrSges(4,1) * t270 + qJDD(1) * t275 + t257;
t234 = -qJDD(1) * pkin(1) + t259;
t233 = t255 * pkin(1) - t278;
t218 = Ifges(5,1) * t236 + Ifges(5,4) * t235 + Ifges(5,5) * qJD(4);
t217 = Ifges(5,4) * t236 + Ifges(5,2) * t235 + Ifges(5,6) * qJD(4);
t216 = Ifges(5,5) * t236 + Ifges(5,6) * t235 + Ifges(5,3) * qJD(4);
t204 = mrSges(5,2) * t215 - mrSges(5,3) * t210 + Ifges(5,1) * t225 + Ifges(5,4) * t224 + Ifges(5,5) * qJDD(4) - qJD(4) * t217 + t235 * t216;
t203 = -mrSges(5,1) * t215 + mrSges(5,3) * t211 + Ifges(5,4) * t225 + Ifges(5,2) * t224 + Ifges(5,6) * qJDD(4) + qJD(4) * t218 - t236 * t216;
t200 = m(3) * t234 + qJDD(1) * mrSges(3,2) - t255 * mrSges(3,3) + t264;
t1 = [mrSges(2,1) * t268 - mrSges(2,2) * t265 + mrSges(3,2) * t234 - mrSges(3,3) * t233 + t250 * (mrSges(4,2) * t229 - mrSges(4,3) * t222 - pkin(5) * t274 - t251 * t203 + t253 * t204) - t249 * (-mrSges(4,1) * t229 + mrSges(4,3) * t223 - pkin(3) * t257 + pkin(5) * t266 + t253 * t203 + t251 * t204) - qJ(3) * t264 - pkin(1) * t200 + qJ(2) * (-m(3) * t233 + (mrSges(3,2) + t267) * t255 + t256) + (Ifges(4,1) * t279 + qJ(2) * mrSges(3,3) + Ifges(3,1) + Ifges(2,3) + (-0.2e1 * Ifges(4,4) * t250 + Ifges(4,2) * t249) * t249) * qJDD(1); t200; t255 * t267 + t256; mrSges(5,1) * t210 - mrSges(5,2) * t211 + Ifges(5,5) * t225 + Ifges(5,6) * t224 + Ifges(5,3) * qJDD(4) + t236 * t217 - t235 * t218;];
tauJ = t1;
