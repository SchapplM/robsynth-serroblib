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
% Datum: 2019-12-05 18:15
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 18:14:29
% EndTime: 2019-12-05 18:14:29
% DurationCPUTime: 0.36s
% Computational Cost: add. (3863->112), mult. (5048->148), div. (0->0), fcn. (2712->10), ass. (0->60)
t264 = qJD(1) + qJD(3);
t258 = qJD(4) + t264;
t268 = sin(qJ(5));
t287 = t258 * t268;
t272 = cos(qJ(5));
t286 = t258 * t272;
t271 = sin(qJ(1));
t275 = cos(qJ(1));
t283 = t275 * g(2) + t271 * g(3);
t252 = qJDD(1) * pkin(1) + t283;
t276 = qJD(1) ^ 2;
t281 = t271 * g(2) - t275 * g(3);
t253 = -t276 * pkin(1) + t281;
t266 = sin(pkin(9));
t267 = cos(pkin(9));
t237 = t267 * t252 - t266 * t253;
t235 = qJDD(1) * pkin(2) + t237;
t238 = t266 * t252 + t267 * t253;
t236 = -t276 * pkin(2) + t238;
t270 = sin(qJ(3));
t274 = cos(qJ(3));
t230 = t274 * t235 - t270 * t236;
t262 = t264 ^ 2;
t263 = qJDD(1) + qJDD(3);
t227 = t263 * pkin(3) + t230;
t231 = t270 * t235 + t274 * t236;
t228 = -t262 * pkin(3) + t231;
t269 = sin(qJ(4));
t273 = cos(qJ(4));
t224 = t269 * t227 + t273 * t228;
t256 = t258 ^ 2;
t257 = qJDD(4) + t263;
t221 = -t256 * pkin(4) + t257 * pkin(8) + t224;
t265 = -g(1) + qJDD(2);
t218 = -t268 * t221 + t272 * t265;
t244 = (-mrSges(6,1) * t272 + mrSges(6,2) * t268) * t258;
t282 = qJD(5) * t258;
t245 = t268 * t257 + t272 * t282;
t251 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t286;
t216 = m(6) * t218 + qJDD(5) * mrSges(6,1) - t245 * mrSges(6,3) + qJD(5) * t251 - t244 * t287;
t219 = t272 * t221 + t268 * t265;
t246 = t272 * t257 - t268 * t282;
t250 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t287;
t217 = m(6) * t219 - qJDD(5) * mrSges(6,2) + t246 * mrSges(6,3) - qJD(5) * t250 + t244 * t286;
t280 = -t268 * t216 + t272 * t217;
t208 = m(5) * t224 - t256 * mrSges(5,1) - t257 * mrSges(5,2) + t280;
t223 = t273 * t227 - t269 * t228;
t220 = -t257 * pkin(4) - t256 * pkin(8) - t223;
t278 = -m(6) * t220 + t246 * mrSges(6,1) - t245 * mrSges(6,2) - t250 * t287 + t251 * t286;
t213 = m(5) * t223 + t257 * mrSges(5,1) - t256 * mrSges(5,2) + t278;
t284 = t269 * t208 + t273 * t213;
t204 = m(4) * t230 + t263 * mrSges(4,1) - t262 * mrSges(4,2) + t284;
t205 = m(4) * t231 - t262 * mrSges(4,1) - t263 * mrSges(4,2) + t273 * t208 - t269 * t213;
t285 = t274 * t204 + t270 * t205;
t239 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t268 + Ifges(6,6) * t272) * t258;
t240 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t268 + Ifges(6,2) * t272) * t258;
t241 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t268 + Ifges(6,4) * t272) * t258;
t279 = -mrSges(5,2) * t224 + pkin(8) * t280 + t268 * (mrSges(6,2) * t220 - mrSges(6,3) * t218 + Ifges(6,1) * t245 + Ifges(6,4) * t246 + Ifges(6,5) * qJDD(5) - qJD(5) * t240 + t239 * t286) + t272 * (-mrSges(6,1) * t220 + mrSges(6,3) * t219 + Ifges(6,4) * t245 + Ifges(6,2) * t246 + Ifges(6,6) * qJDD(5) + qJD(5) * t241 - t239 * t287) + pkin(4) * t278 + mrSges(5,1) * t223 + Ifges(5,3) * t257;
t277 = mrSges(4,1) * t230 - mrSges(4,2) * t231 + Ifges(4,3) * t263 + pkin(3) * t284 + t279;
t1 = [pkin(2) * t285 + t277 + Ifges(2,3) * qJDD(1) + Ifges(3,3) * qJDD(1) + mrSges(2,1) * t283 + mrSges(3,1) * t237 - mrSges(3,2) * t238 - mrSges(2,2) * t281 + pkin(1) * (t266 * (m(3) * t238 - t276 * mrSges(3,1) - qJDD(1) * mrSges(3,2) - t270 * t204 + t274 * t205) + t267 * (m(3) * t237 + qJDD(1) * mrSges(3,1) - t276 * mrSges(3,2) + t285)); t272 * t216 + t268 * t217 + (m(3) + m(4) + m(5)) * t265; t277; t279; mrSges(6,1) * t218 - mrSges(6,2) * t219 + Ifges(6,5) * t245 + Ifges(6,6) * t246 + Ifges(6,3) * qJDD(5) + (t240 * t268 - t241 * t272) * t258;];
tauJ = t1;
