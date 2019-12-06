% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRRPR2
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-05 18:41
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRRPR2_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR2_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR2_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR2_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR2_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR2_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR2_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR2_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR2_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:40:52
% EndTime: 2019-12-05 18:40:53
% DurationCPUTime: 0.46s
% Computational Cost: add. (5275->110), mult. (6167->147), div. (0->0), fcn. (3316->10), ass. (0->59)
t263 = qJD(1) + qJD(2);
t257 = qJD(3) + t263;
t267 = sin(qJ(5));
t285 = t257 * t267;
t271 = cos(qJ(5));
t284 = t257 * t271;
t270 = sin(qJ(1));
t274 = cos(qJ(1));
t281 = t274 * g(2) + t270 * g(3);
t250 = qJDD(1) * pkin(1) + t281;
t279 = t270 * g(2) - t274 * g(3);
t251 = -qJD(1) ^ 2 * pkin(1) + t279;
t269 = sin(qJ(2));
t273 = cos(qJ(2));
t235 = t273 * t250 - t269 * t251;
t262 = qJDD(1) + qJDD(2);
t232 = t262 * pkin(2) + t235;
t236 = t269 * t250 + t273 * t251;
t261 = t263 ^ 2;
t233 = -t261 * pkin(2) + t236;
t268 = sin(qJ(3));
t272 = cos(qJ(3));
t227 = t272 * t232 - t268 * t233;
t255 = t257 ^ 2;
t256 = qJDD(3) + t262;
t224 = t256 * pkin(3) + t227;
t228 = t268 * t232 + t272 * t233;
t225 = -t255 * pkin(3) + t228;
t265 = sin(pkin(9));
t266 = cos(pkin(9));
t221 = t265 * t224 + t266 * t225;
t218 = -t255 * pkin(4) + t256 * pkin(8) + t221;
t264 = -g(1) + qJDD(4);
t215 = -t267 * t218 + t271 * t264;
t242 = (-mrSges(6,1) * t271 + mrSges(6,2) * t267) * t257;
t280 = qJD(5) * t257;
t243 = t267 * t256 + t271 * t280;
t249 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t284;
t213 = m(6) * t215 + qJDD(5) * mrSges(6,1) - t243 * mrSges(6,3) + qJD(5) * t249 - t242 * t285;
t216 = t271 * t218 + t267 * t264;
t244 = t271 * t256 - t267 * t280;
t248 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t285;
t214 = m(6) * t216 - qJDD(5) * mrSges(6,2) + t244 * mrSges(6,3) - qJD(5) * t248 + t242 * t284;
t278 = -t267 * t213 + t271 * t214;
t205 = m(5) * t221 - t255 * mrSges(5,1) - t256 * mrSges(5,2) + t278;
t220 = t266 * t224 - t265 * t225;
t217 = -t256 * pkin(4) - t255 * pkin(8) - t220;
t277 = -m(6) * t217 + t244 * mrSges(6,1) - t243 * mrSges(6,2) - t248 * t285 + t249 * t284;
t210 = m(5) * t220 + t256 * mrSges(5,1) - t255 * mrSges(5,2) + t277;
t282 = t265 * t205 + t266 * t210;
t201 = m(4) * t227 + t256 * mrSges(4,1) - t255 * mrSges(4,2) + t282;
t202 = m(4) * t228 - t255 * mrSges(4,1) - t256 * mrSges(4,2) + t266 * t205 - t265 * t210;
t283 = t272 * t201 + t268 * t202;
t237 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t267 + Ifges(6,6) * t271) * t257;
t238 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t267 + Ifges(6,2) * t271) * t257;
t239 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t267 + Ifges(6,4) * t271) * t257;
t276 = -mrSges(4,2) * t228 - mrSges(5,2) * t221 + pkin(3) * t282 + pkin(8) * t278 + t267 * (mrSges(6,2) * t217 - mrSges(6,3) * t215 + Ifges(6,1) * t243 + Ifges(6,4) * t244 + Ifges(6,5) * qJDD(5) - qJD(5) * t238 + t237 * t284) + t271 * (-mrSges(6,1) * t217 + mrSges(6,3) * t216 + Ifges(6,4) * t243 + Ifges(6,2) * t244 + Ifges(6,6) * qJDD(5) + qJD(5) * t239 - t237 * t285) + pkin(4) * t277 + mrSges(5,1) * t220 + mrSges(4,1) * t227 + (Ifges(5,3) + Ifges(4,3)) * t256;
t275 = mrSges(3,1) * t235 - mrSges(3,2) * t236 + Ifges(3,3) * t262 + pkin(2) * t283 + t276;
t1 = [t275 + Ifges(2,3) * qJDD(1) + mrSges(2,1) * t281 + pkin(1) * (t269 * (m(3) * t236 - t261 * mrSges(3,1) - t262 * mrSges(3,2) - t268 * t201 + t272 * t202) + t273 * (m(3) * t235 + t262 * mrSges(3,1) - t261 * mrSges(3,2) + t283)) - mrSges(2,2) * t279; t275; t276; m(5) * t264 + t271 * t213 + t267 * t214; mrSges(6,1) * t215 - mrSges(6,2) * t216 + Ifges(6,5) * t243 + Ifges(6,6) * t244 + Ifges(6,3) * qJDD(5) + (t238 * t267 - t239 * t271) * t257;];
tauJ = t1;
