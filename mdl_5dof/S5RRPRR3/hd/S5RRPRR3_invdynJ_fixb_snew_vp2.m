% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRPRR3
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2022-01-20 10:34
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRPRR3_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR3_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR3_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR3_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR3_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR3_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR3_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR3_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR3_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:34:09
% EndTime: 2022-01-20 10:34:10
% DurationCPUTime: 0.48s
% Computational Cost: add. (4776->111), mult. (5787->147), div. (0->0), fcn. (3110->10), ass. (0->59)
t260 = qJD(1) + qJD(2);
t256 = qJD(4) + t260;
t264 = sin(qJ(5));
t282 = t256 * t264;
t268 = cos(qJ(5));
t281 = t256 * t268;
t267 = sin(qJ(1));
t271 = cos(qJ(1));
t277 = t267 * g(1) - t271 * g(2);
t249 = qJDD(1) * pkin(1) + t277;
t275 = -t271 * g(1) - t267 * g(2);
t250 = -qJD(1) ^ 2 * pkin(1) + t275;
t266 = sin(qJ(2));
t270 = cos(qJ(2));
t234 = t270 * t249 - t266 * t250;
t259 = qJDD(1) + qJDD(2);
t231 = t259 * pkin(2) + t234;
t235 = t266 * t249 + t270 * t250;
t258 = t260 ^ 2;
t232 = -t258 * pkin(2) + t235;
t262 = sin(pkin(9));
t263 = cos(pkin(9));
t226 = t263 * t231 - t262 * t232;
t223 = t259 * pkin(3) + t226;
t227 = t262 * t231 + t263 * t232;
t224 = -t258 * pkin(3) + t227;
t265 = sin(qJ(4));
t269 = cos(qJ(4));
t220 = t265 * t223 + t269 * t224;
t254 = t256 ^ 2;
t255 = qJDD(4) + t259;
t217 = -t254 * pkin(4) + t255 * pkin(8) + t220;
t261 = -g(3) + qJDD(3);
t214 = -t264 * t217 + t268 * t261;
t241 = (-mrSges(6,1) * t268 + mrSges(6,2) * t264) * t256;
t278 = qJD(5) * t256;
t242 = t264 * t255 + t268 * t278;
t248 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t281;
t212 = m(6) * t214 + qJDD(5) * mrSges(6,1) - t242 * mrSges(6,3) + qJD(5) * t248 - t241 * t282;
t215 = t268 * t217 + t264 * t261;
t243 = t268 * t255 - t264 * t278;
t247 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t282;
t213 = m(6) * t215 - qJDD(5) * mrSges(6,2) + t243 * mrSges(6,3) - qJD(5) * t247 + t241 * t281;
t276 = -t264 * t212 + t268 * t213;
t204 = m(5) * t220 - t254 * mrSges(5,1) - t255 * mrSges(5,2) + t276;
t219 = t269 * t223 - t265 * t224;
t216 = -t255 * pkin(4) - t254 * pkin(8) - t219;
t273 = -m(6) * t216 + t243 * mrSges(6,1) - t242 * mrSges(6,2) - t247 * t282 + t248 * t281;
t209 = m(5) * t219 + t255 * mrSges(5,1) - t254 * mrSges(5,2) + t273;
t279 = t265 * t204 + t269 * t209;
t200 = m(4) * t226 + t259 * mrSges(4,1) - t258 * mrSges(4,2) + t279;
t201 = m(4) * t227 - t258 * mrSges(4,1) - t259 * mrSges(4,2) + t269 * t204 - t265 * t209;
t280 = t263 * t200 + t262 * t201;
t236 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t264 + Ifges(6,6) * t268) * t256;
t237 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t264 + Ifges(6,2) * t268) * t256;
t238 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t264 + Ifges(6,4) * t268) * t256;
t274 = -mrSges(5,2) * t220 + pkin(8) * t276 + t264 * (mrSges(6,2) * t216 - mrSges(6,3) * t214 + Ifges(6,1) * t242 + Ifges(6,4) * t243 + Ifges(6,5) * qJDD(5) - qJD(5) * t237 + t236 * t281) + t268 * (-mrSges(6,1) * t216 + mrSges(6,3) * t215 + Ifges(6,4) * t242 + Ifges(6,2) * t243 + Ifges(6,6) * qJDD(5) + qJD(5) * t238 - t236 * t282) + pkin(4) * t273 + mrSges(5,1) * t219 + Ifges(5,3) * t255;
t272 = mrSges(3,1) * t234 + mrSges(4,1) * t226 - mrSges(3,2) * t235 - mrSges(4,2) * t227 + pkin(2) * t280 + pkin(3) * t279 + t274 + (Ifges(3,3) + Ifges(4,3)) * t259;
t1 = [Ifges(2,3) * qJDD(1) + pkin(1) * (t266 * (m(3) * t235 - t258 * mrSges(3,1) - t259 * mrSges(3,2) - t262 * t200 + t263 * t201) + t270 * (m(3) * t234 + t259 * mrSges(3,1) - t258 * mrSges(3,2) + t280)) + mrSges(2,1) * t277 - mrSges(2,2) * t275 + t272; t272; t268 * t212 + t264 * t213 + (m(4) + m(5)) * t261; t274; mrSges(6,1) * t214 - mrSges(6,2) * t215 + Ifges(6,5) * t242 + Ifges(6,6) * t243 + Ifges(6,3) * qJDD(5) + (t237 * t264 - t238 * t268) * t256;];
tauJ = t1;
