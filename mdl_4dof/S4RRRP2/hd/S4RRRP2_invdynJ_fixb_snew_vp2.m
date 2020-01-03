% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S4RRRP2
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
%   pkin=[a2,a3,a4,d1,d2,d3]';
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
% Datum: 2019-12-31 17:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S4RRRP2_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP2_invdynJ_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP2_invdynJ_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP2_invdynJ_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP2_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP2_invdynJ_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP2_invdynJ_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRP2_invdynJ_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRP2_invdynJ_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:12:56
% EndTime: 2019-12-31 17:12:57
% DurationCPUTime: 0.38s
% Computational Cost: add. (1402->127), mult. (1779->155), div. (0->0), fcn. (796->6), ass. (0->59)
t285 = Ifges(4,1) + Ifges(5,1);
t279 = Ifges(4,4) + Ifges(5,4);
t278 = Ifges(4,5) + Ifges(5,5);
t284 = Ifges(4,2) + Ifges(5,2);
t277 = Ifges(4,6) + Ifges(5,6);
t283 = Ifges(4,3) + Ifges(5,3);
t250 = qJD(1) + qJD(2);
t248 = t250 ^ 2;
t282 = pkin(3) * t248;
t256 = cos(qJ(3));
t281 = t256 * g(3);
t280 = -mrSges(4,2) - mrSges(5,2);
t253 = sin(qJ(3));
t276 = t250 * t253;
t275 = t250 * t256;
t274 = (t253 * t278 + t277 * t256) * t250 + t283 * qJD(3);
t273 = (t253 * t279 + t284 * t256) * t250 + t277 * qJD(3);
t272 = (-t285 * t253 - t256 * t279) * t250 - t278 * qJD(3);
t255 = sin(qJ(1));
t258 = cos(qJ(1));
t265 = t255 * g(1) - t258 * g(2);
t239 = qJDD(1) * pkin(1) + t265;
t262 = -t258 * g(1) - t255 * g(2);
t240 = -qJD(1) ^ 2 * pkin(1) + t262;
t254 = sin(qJ(2));
t257 = cos(qJ(2));
t218 = t254 * t239 + t257 * t240;
t242 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t276;
t271 = -qJD(3) * mrSges(4,1) + mrSges(4,3) * t276 - t242;
t270 = qJD(3) * t250;
t269 = qJD(4) * t250;
t249 = qJDD(1) + qJDD(2);
t215 = -t248 * pkin(2) + t249 * pkin(6) + t218;
t266 = t256 * t270;
t233 = t253 * t249 + t266;
t208 = qJDD(3) * pkin(3) - t281 + (-t233 + t266) * qJ(4) + (t256 * t282 - t215 - 0.2e1 * t269) * t253;
t244 = -qJD(3) * mrSges(5,2) + mrSges(5,3) * t275;
t268 = m(5) * t208 + qJDD(3) * mrSges(5,1) + qJD(3) * t244;
t212 = -t253 * g(3) + t256 * t215;
t234 = t256 * t249 - t253 * t270;
t241 = qJD(3) * pkin(3) - qJ(4) * t276;
t252 = t256 ^ 2;
t209 = t234 * qJ(4) - qJD(3) * t241 - t252 * t282 + 0.2e1 * t256 * t269 + t212;
t231 = (-mrSges(5,1) * t256 + mrSges(5,2) * t253) * t250;
t267 = m(5) * t209 + t234 * mrSges(5,3) + t231 * t275;
t211 = -t253 * t215 - t281;
t232 = (-mrSges(4,1) * t256 + mrSges(4,2) * t253) * t250;
t245 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t275;
t264 = -t253 * (m(4) * t211 + qJDD(3) * mrSges(4,1) + qJD(3) * t245 + (-t231 - t232) * t276 + (-mrSges(4,3) - mrSges(5,3)) * t233 + t268) + t256 * (m(4) * t212 + t234 * mrSges(4,3) + t271 * qJD(3) + t280 * qJDD(3) + t232 * t275 + t267);
t217 = t257 * t239 - t254 * t240;
t261 = -t249 * pkin(2) - t217;
t210 = t241 * t276 - t234 * pkin(3) + qJDD(4) + (-qJ(4) * t252 - pkin(6)) * t248 + t261;
t263 = -m(5) * t210 + t234 * mrSges(5,1) + t244 * t275;
t204 = t233 * mrSges(5,2) + t242 * t276 - t263;
t205 = -t233 * mrSges(5,3) - t231 * t276 + t268;
t214 = -t248 * pkin(6) + t261;
t259 = -m(4) * t214 + t234 * mrSges(4,1) + t280 * t233 + t245 * t275 + t271 * t276 + t263;
t260 = -mrSges(3,2) * t218 + t256 * (-mrSges(4,1) * t214 + mrSges(4,3) * t212 - mrSges(5,1) * t210 + mrSges(5,3) * t209 - pkin(3) * t204 + qJ(4) * t267 - t274 * t276 + t284 * t234 + t279 * t233 + (-qJ(4) * mrSges(5,2) + t277) * qJDD(3) + (-qJ(4) * t242 - t272) * qJD(3)) + t253 * (mrSges(4,2) * t214 + mrSges(5,2) * t210 - mrSges(4,3) * t211 - mrSges(5,3) * t208 - qJ(4) * t205 - t273 * qJD(3) + t278 * qJDD(3) + t285 * t233 + t279 * t234 + t274 * t275) + pkin(6) * t264 + pkin(2) * t259 + mrSges(3,1) * t217 + Ifges(3,3) * t249;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t265 - mrSges(2,2) * t262 + pkin(1) * (t254 * (m(3) * t218 - t248 * mrSges(3,1) - t249 * mrSges(3,2) + t264) + t257 * (m(3) * t217 + t249 * mrSges(3,1) - t248 * mrSges(3,2) + t259)) + t260; t260; mrSges(4,1) * t211 + mrSges(5,1) * t208 - mrSges(4,2) * t212 - mrSges(5,2) * t209 + pkin(3) * t205 + t277 * t234 + t278 * t233 + t283 * qJDD(3) + (t273 * t253 + t272 * t256) * t250; t204;];
tauJ = t1;
