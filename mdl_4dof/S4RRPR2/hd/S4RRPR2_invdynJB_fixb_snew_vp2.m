% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S4RRPR2
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
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
% tauJB [(6+4)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 18:16
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S4RRPR2_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR2_invdynJB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR2_invdynJB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR2_invdynJB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR2_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPR2_invdynJB_fixb_snew_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR2_invdynJB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR2_invdynJB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR2_invdynJB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 18:16:33
% EndTime: 2019-07-18 18:16:33
% DurationCPUTime: 0.39s
% Computational Cost: add. (4530->118), mult. (5160->132), div. (0->0), fcn. (2118->6), ass. (0->57)
t334 = -m(4) - m(5);
t333 = -pkin(2) - pkin(3);
t332 = -mrSges(3,1) - mrSges(4,1);
t331 = Ifges(4,4) + Ifges(3,5);
t330 = Ifges(3,6) - Ifges(4,6);
t315 = sin(qJ(1));
t318 = cos(qJ(1));
t295 = t315 * g(1) - t318 * g(2);
t292 = qJDD(1) * pkin(1) + t295;
t296 = -t318 * g(1) - t315 * g(2);
t319 = qJD(1) ^ 2;
t293 = -t319 * pkin(1) + t296;
t314 = sin(qJ(2));
t317 = cos(qJ(2));
t288 = t314 * t292 + t317 * t293;
t311 = (qJD(1) + qJD(2));
t309 = t311 ^ 2;
t310 = qJDD(1) + qJDD(2);
t325 = t310 * qJ(3) + (2 * qJD(3) * t311) + t288;
t280 = t333 * t309 + t325;
t287 = t317 * t292 - t314 * t293;
t322 = -t309 * qJ(3) + qJDD(3) - t287;
t283 = t333 * t310 + t322;
t313 = sin(qJ(4));
t316 = cos(qJ(4));
t278 = -t313 * t280 + t316 * t283;
t307 = qJD(4) - t311;
t304 = t307 ^ 2;
t305 = qJDD(4) - t310;
t275 = m(5) * t278 + t305 * mrSges(5,1) - (t304 * mrSges(5,2));
t279 = t316 * t280 + t313 * t283;
t276 = m(5) * t279 - t304 * mrSges(5,1) - t305 * mrSges(5,2);
t284 = -t309 * pkin(2) + t325;
t324 = m(4) * t284 + t310 * mrSges(4,3) - t313 * t275 + t316 * t276;
t265 = m(3) * t288 - t310 * mrSges(3,2) + t332 * t309 + t324;
t271 = t316 * t275 + t313 * t276;
t285 = -t310 * pkin(2) + t322;
t270 = m(4) * t285 - t310 * mrSges(4,1) - t309 * mrSges(4,3) + t271;
t267 = m(3) * t287 + t310 * mrSges(3,1) - t309 * mrSges(3,2) - t270;
t260 = t314 * t265 + t317 * t267;
t257 = m(2) * t295 + qJDD(1) * mrSges(2,1) - t319 * mrSges(2,2) + t260;
t326 = t317 * t265 - t314 * t267;
t258 = m(2) * t296 - t319 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t326;
t329 = t318 * t257 + t315 * t258;
t328 = -m(3) + t334;
t327 = -t315 * t257 + t318 * t258;
t323 = mrSges(5,1) * t278 - mrSges(5,2) * t279 + Ifges(5,3) * t305;
t321 = -mrSges(4,1) * t285 - mrSges(3,2) * t288 - pkin(3) * t271 + qJ(3) * (-t309 * mrSges(4,1) + t324) - pkin(2) * t270 + mrSges(4,3) * t284 + mrSges(3,1) * t287 - t323 + (Ifges(3,3) + Ifges(4,2)) * t310;
t320 = mrSges(2,1) * t295 - mrSges(2,2) * t296 + Ifges(2,3) * qJDD(1) + pkin(1) * t260 + t321;
t299 = t334 * g(3);
t274 = mrSges(5,2) * g(3) - mrSges(5,3) * t278 + Ifges(5,5) * t305 - (t304 * Ifges(5,6));
t273 = -mrSges(5,1) * g(3) + mrSges(5,3) * t279 + t304 * Ifges(5,5) + Ifges(5,6) * t305;
t262 = mrSges(4,2) * t285 - mrSges(3,3) * t287 - qJ(3) * t299 - t313 * t273 + t316 * t274 + t331 * t310 - t330 * t309 + (-mrSges(3,2) + mrSges(4,3)) * g(3);
t261 = mrSges(4,2) * t284 + mrSges(3,3) * t288 - pkin(2) * t299 - t316 * t273 - t313 * t274 + t330 * t310 + t331 * t309 + (m(5) * pkin(3) - t332) * g(3);
t253 = -mrSges(2,2) * g(3) - mrSges(2,3) * t295 + Ifges(2,5) * qJDD(1) - t319 * Ifges(2,6) - pkin(5) * t260 - t314 * t261 + t317 * t262;
t252 = Ifges(2,6) * qJDD(1) + t319 * Ifges(2,5) + mrSges(2,3) * t296 + t314 * t262 + t317 * t261 + pkin(5) * t326 + (-pkin(1) * t328 + mrSges(2,1)) * g(3);
t1 = [-m(1) * g(1) + t327; -m(1) * g(2) + t329; (-m(1) - m(2) + t328) * g(3); -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(4) * t329 - t315 * t252 + t318 * t253; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(4) * t327 + t318 * t252 + t315 * t253; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t320; t320; t321; t270; t323;];
tauJB  = t1;
