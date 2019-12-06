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
% Datum: 2019-12-05 18:31
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 18:30:26
% EndTime: 2019-12-05 18:30:27
% DurationCPUTime: 0.44s
% Computational Cost: add. (4776->111), mult. (5787->147), div. (0->0), fcn. (3110->10), ass. (0->59)
t266 = qJD(1) + qJD(2);
t260 = qJD(4) + t266;
t270 = sin(qJ(5));
t288 = t260 * t270;
t274 = cos(qJ(5));
t287 = t260 * t274;
t273 = sin(qJ(1));
t277 = cos(qJ(1));
t284 = t277 * g(2) + t273 * g(3);
t253 = qJDD(1) * pkin(1) + t284;
t282 = t273 * g(2) - t277 * g(3);
t254 = -qJD(1) ^ 2 * pkin(1) + t282;
t272 = sin(qJ(2));
t276 = cos(qJ(2));
t238 = t276 * t253 - t272 * t254;
t265 = qJDD(1) + qJDD(2);
t235 = t265 * pkin(2) + t238;
t239 = t272 * t253 + t276 * t254;
t264 = t266 ^ 2;
t236 = -t264 * pkin(2) + t239;
t268 = sin(pkin(9));
t269 = cos(pkin(9));
t230 = t269 * t235 - t268 * t236;
t227 = t265 * pkin(3) + t230;
t231 = t268 * t235 + t269 * t236;
t228 = -t264 * pkin(3) + t231;
t271 = sin(qJ(4));
t275 = cos(qJ(4));
t224 = t271 * t227 + t275 * t228;
t258 = t260 ^ 2;
t259 = qJDD(4) + t265;
t221 = -t258 * pkin(4) + t259 * pkin(8) + t224;
t267 = -g(1) + qJDD(3);
t218 = -t270 * t221 + t274 * t267;
t245 = (-mrSges(6,1) * t274 + mrSges(6,2) * t270) * t260;
t283 = qJD(5) * t260;
t246 = t270 * t259 + t274 * t283;
t252 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t287;
t216 = m(6) * t218 + qJDD(5) * mrSges(6,1) - t246 * mrSges(6,3) + qJD(5) * t252 - t245 * t288;
t219 = t274 * t221 + t270 * t267;
t247 = t274 * t259 - t270 * t283;
t251 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t288;
t217 = m(6) * t219 - qJDD(5) * mrSges(6,2) + t247 * mrSges(6,3) - qJD(5) * t251 + t245 * t287;
t281 = -t270 * t216 + t274 * t217;
t208 = m(5) * t224 - t258 * mrSges(5,1) - t259 * mrSges(5,2) + t281;
t223 = t275 * t227 - t271 * t228;
t220 = -t259 * pkin(4) - t258 * pkin(8) - t223;
t279 = -m(6) * t220 + t247 * mrSges(6,1) - t246 * mrSges(6,2) - t251 * t288 + t252 * t287;
t213 = m(5) * t223 + t259 * mrSges(5,1) - t258 * mrSges(5,2) + t279;
t285 = t271 * t208 + t275 * t213;
t204 = m(4) * t230 + t265 * mrSges(4,1) - t264 * mrSges(4,2) + t285;
t205 = m(4) * t231 - t264 * mrSges(4,1) - t265 * mrSges(4,2) + t275 * t208 - t271 * t213;
t286 = t269 * t204 + t268 * t205;
t240 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t270 + Ifges(6,6) * t274) * t260;
t241 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t270 + Ifges(6,2) * t274) * t260;
t242 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t270 + Ifges(6,4) * t274) * t260;
t280 = -mrSges(5,2) * t224 + pkin(8) * t281 + t270 * (mrSges(6,2) * t220 - mrSges(6,3) * t218 + Ifges(6,1) * t246 + Ifges(6,4) * t247 + Ifges(6,5) * qJDD(5) - qJD(5) * t241 + t240 * t287) + t274 * (-mrSges(6,1) * t220 + mrSges(6,3) * t219 + Ifges(6,4) * t246 + Ifges(6,2) * t247 + Ifges(6,6) * qJDD(5) + qJD(5) * t242 - t240 * t288) + pkin(4) * t279 + mrSges(5,1) * t223 + Ifges(5,3) * t259;
t278 = mrSges(3,1) * t238 + mrSges(4,1) * t230 - mrSges(3,2) * t239 - mrSges(4,2) * t231 + pkin(2) * t286 + pkin(3) * t285 + t280 + (Ifges(3,3) + Ifges(4,3)) * t265;
t1 = [mrSges(2,1) * t284 + pkin(1) * (t272 * (m(3) * t239 - t264 * mrSges(3,1) - t265 * mrSges(3,2) - t268 * t204 + t269 * t205) + t276 * (m(3) * t238 + t265 * mrSges(3,1) - t264 * mrSges(3,2) + t286)) - mrSges(2,2) * t282 + t278 + Ifges(2,3) * qJDD(1); t278; t274 * t216 + t270 * t217 + (m(4) + m(5)) * t267; t280; mrSges(6,1) * t218 - mrSges(6,2) * t219 + Ifges(6,5) * t246 + Ifges(6,6) * t247 + Ifges(6,3) * qJDD(5) + (t241 * t270 - t242 * t274) * t260;];
tauJ = t1;
