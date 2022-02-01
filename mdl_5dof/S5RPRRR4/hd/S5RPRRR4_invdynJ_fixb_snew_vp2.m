% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RPRRR4
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% m [6x1]
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
% Datum: 2022-01-23 09:35
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RPRRR4_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR4_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR4_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR4_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR4_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR4_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR4_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR4_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR4_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:34:36
% EndTime: 2022-01-23 09:34:36
% DurationCPUTime: 0.38s
% Computational Cost: add. (3863->112), mult. (5048->148), div. (0->0), fcn. (2712->10), ass. (0->60)
t258 = qJD(1) + qJD(3);
t254 = qJD(4) + t258;
t262 = sin(qJ(5));
t281 = t254 * t262;
t266 = cos(qJ(5));
t280 = t254 * t266;
t265 = sin(qJ(1));
t269 = cos(qJ(1));
t276 = t265 * g(1) - t269 * g(2);
t248 = qJDD(1) * pkin(1) + t276;
t270 = qJD(1) ^ 2;
t274 = -t269 * g(1) - t265 * g(2);
t249 = -t270 * pkin(1) + t274;
t260 = sin(pkin(9));
t261 = cos(pkin(9));
t233 = t261 * t248 - t260 * t249;
t231 = qJDD(1) * pkin(2) + t233;
t234 = t260 * t248 + t261 * t249;
t232 = -t270 * pkin(2) + t234;
t264 = sin(qJ(3));
t268 = cos(qJ(3));
t226 = t268 * t231 - t264 * t232;
t256 = t258 ^ 2;
t257 = qJDD(1) + qJDD(3);
t223 = t257 * pkin(3) + t226;
t227 = t264 * t231 + t268 * t232;
t224 = -t256 * pkin(3) + t227;
t263 = sin(qJ(4));
t267 = cos(qJ(4));
t220 = t263 * t223 + t267 * t224;
t252 = t254 ^ 2;
t253 = qJDD(4) + t257;
t217 = -t252 * pkin(4) + t253 * pkin(8) + t220;
t259 = -g(3) + qJDD(2);
t214 = -t262 * t217 + t266 * t259;
t240 = (-mrSges(6,1) * t266 + mrSges(6,2) * t262) * t254;
t277 = qJD(5) * t254;
t241 = t262 * t253 + t266 * t277;
t247 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t280;
t212 = m(6) * t214 + qJDD(5) * mrSges(6,1) - t241 * mrSges(6,3) + qJD(5) * t247 - t240 * t281;
t215 = t266 * t217 + t262 * t259;
t242 = t266 * t253 - t262 * t277;
t246 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t281;
t213 = m(6) * t215 - qJDD(5) * mrSges(6,2) + t242 * mrSges(6,3) - qJD(5) * t246 + t240 * t280;
t275 = -t262 * t212 + t266 * t213;
t204 = m(5) * t220 - t252 * mrSges(5,1) - t253 * mrSges(5,2) + t275;
t219 = t267 * t223 - t263 * t224;
t216 = -t253 * pkin(4) - t252 * pkin(8) - t219;
t272 = -m(6) * t216 + t242 * mrSges(6,1) - t241 * mrSges(6,2) - t246 * t281 + t247 * t280;
t209 = m(5) * t219 + t253 * mrSges(5,1) - t252 * mrSges(5,2) + t272;
t278 = t263 * t204 + t267 * t209;
t200 = m(4) * t226 + t257 * mrSges(4,1) - t256 * mrSges(4,2) + t278;
t201 = m(4) * t227 - t256 * mrSges(4,1) - t257 * mrSges(4,2) + t267 * t204 - t263 * t209;
t279 = t268 * t200 + t264 * t201;
t235 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t262 + Ifges(6,6) * t266) * t254;
t236 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t262 + Ifges(6,2) * t266) * t254;
t237 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t262 + Ifges(6,4) * t266) * t254;
t273 = -mrSges(5,2) * t220 + pkin(8) * t275 + t262 * (mrSges(6,2) * t216 - mrSges(6,3) * t214 + Ifges(6,1) * t241 + Ifges(6,4) * t242 + Ifges(6,5) * qJDD(5) - qJD(5) * t236 + t235 * t280) + t266 * (-mrSges(6,1) * t216 + mrSges(6,3) * t215 + Ifges(6,4) * t241 + Ifges(6,2) * t242 + Ifges(6,6) * qJDD(5) + qJD(5) * t237 - t235 * t281) + pkin(4) * t272 + mrSges(5,1) * t219 + Ifges(5,3) * t253;
t271 = mrSges(4,1) * t226 - mrSges(4,2) * t227 + Ifges(4,3) * t257 + pkin(3) * t278 + t273;
t1 = [t271 + mrSges(3,1) * t233 - mrSges(3,2) * t234 + Ifges(2,3) * qJDD(1) + Ifges(3,3) * qJDD(1) + pkin(2) * t279 + mrSges(2,1) * t276 - mrSges(2,2) * t274 + pkin(1) * (t260 * (m(3) * t234 - t270 * mrSges(3,1) - qJDD(1) * mrSges(3,2) - t264 * t200 + t268 * t201) + t261 * (m(3) * t233 + qJDD(1) * mrSges(3,1) - t270 * mrSges(3,2) + t279)); t266 * t212 + t262 * t213 + (m(3) + m(4) + m(5)) * t259; t271; t273; mrSges(6,1) * t214 - mrSges(6,2) * t215 + Ifges(6,5) * t241 + Ifges(6,6) * t242 + Ifges(6,3) * qJDD(5) + (t236 * t262 - t237 * t266) * t254;];
tauJ = t1;
