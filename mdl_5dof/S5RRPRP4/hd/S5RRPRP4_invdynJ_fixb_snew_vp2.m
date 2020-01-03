% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRPRP4
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
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
% Datum: 2019-12-31 19:53
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRPRP4_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP4_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP4_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP4_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP4_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP4_invdynJ_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP4_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP4_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP4_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:52:41
% EndTime: 2019-12-31 19:52:42
% DurationCPUTime: 0.48s
% Computational Cost: add. (2039->136), mult. (2415->165), div. (0->0), fcn. (1018->6), ass. (0->59)
t324 = Ifges(5,1) + Ifges(6,1);
t318 = Ifges(5,4) - Ifges(6,5);
t317 = Ifges(5,5) + Ifges(6,4);
t323 = Ifges(5,2) + Ifges(6,3);
t316 = Ifges(5,6) - Ifges(6,6);
t294 = sin(qJ(1));
t297 = cos(qJ(1));
t307 = t294 * g(1) - g(2) * t297;
t273 = qJDD(1) * pkin(1) + t307;
t303 = -g(1) * t297 - g(2) * t294;
t274 = -qJD(1) ^ 2 * pkin(1) + t303;
t293 = sin(qJ(2));
t296 = cos(qJ(2));
t251 = t293 * t273 + t296 * t274;
t288 = qJDD(1) + qJDD(2);
t289 = qJD(1) + qJD(2);
t322 = -t288 * qJ(3) - 0.2e1 * qJD(3) * t289 - t251;
t321 = -pkin(2) - pkin(7);
t292 = sin(qJ(4));
t320 = g(3) * t292;
t319 = -mrSges(5,3) - mrSges(6,2);
t315 = t289 * t292;
t295 = cos(qJ(4));
t314 = t289 * t295;
t313 = (t323 * t292 - t295 * t318) * t289 - t316 * qJD(4);
t312 = (-t292 * t318 + t324 * t295) * t289 + t317 * qJD(4);
t310 = qJD(4) * t289;
t250 = t273 * t296 - t293 * t274;
t287 = t289 ^ 2;
t302 = -qJ(3) * t287 + qJDD(3) - t250;
t245 = t321 * t288 + t302;
t241 = -g(3) * t295 + t292 * t245;
t263 = (pkin(4) * t292 - qJ(5) * t295) * t289;
t298 = qJD(4) ^ 2;
t238 = -pkin(4) * t298 + qJDD(4) * qJ(5) + 0.2e1 * qJD(5) * qJD(4) - t263 * t315 + t241;
t277 = -qJD(4) * mrSges(6,1) + mrSges(6,2) * t314;
t309 = m(6) * t238 + qJDD(4) * mrSges(6,3) + qJD(4) * t277;
t264 = (mrSges(6,1) * t292 - mrSges(6,3) * t295) * t289;
t305 = t289 * (-t264 - (mrSges(5,1) * t292 + mrSges(5,2) * t295) * t289);
t239 = -qJDD(4) * pkin(4) - t320 - qJ(5) * t298 + qJDD(5) + (t263 * t289 - t245) * t295;
t278 = -mrSges(6,2) * t315 + qJD(4) * mrSges(6,3);
t304 = -m(6) * t239 + qJDD(4) * mrSges(6,1) + qJD(4) * t278;
t240 = t245 * t295 + t320;
t266 = t288 * t292 + t295 * t310;
t267 = t288 * t295 - t292 * t310;
t275 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t315;
t276 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t314;
t229 = (m(5) * t241 - qJDD(4) * mrSges(5,2) - qJD(4) * t276 + t319 * t266 + t292 * t305 + t309) * t292 + (m(5) * t240 + qJDD(4) * mrSges(5,1) + qJD(4) * t275 + t319 * t267 + t295 * t305 + t304) * t295;
t243 = t321 * t287 - t322;
t248 = -pkin(2) * t288 + t302;
t301 = -m(4) * t248 + t287 * mrSges(4,3) - t229;
t236 = pkin(4) * t266 - qJ(5) * t267 + (-0.2e1 * qJD(5) * t295 + (pkin(4) * t295 + qJ(5) * t292) * qJD(4)) * t289 + t243;
t233 = m(6) * t236 + t266 * mrSges(6,1) - mrSges(6,3) * t267 - t277 * t314 + t278 * t315;
t228 = mrSges(4,2) * t288 - t301;
t246 = pkin(2) * t287 + t322;
t299 = -m(4) * t246 + m(5) * t243 + mrSges(5,1) * t266 + t287 * mrSges(4,2) + t267 * mrSges(5,2) + t288 * mrSges(4,3) + t275 * t315 + t276 * t314 + t233;
t300 = -mrSges(3,2) * t251 - mrSges(4,3) * t246 - pkin(7) * t229 + t295 * (mrSges(5,2) * t243 + mrSges(6,2) * t239 - mrSges(5,3) * t240 - mrSges(6,3) * t236 - qJ(5) * t233 + t313 * qJD(4) + t317 * qJDD(4) - t318 * t266 + t324 * t267) - pkin(2) * t228 + qJ(3) * t299 + mrSges(4,2) * t248 + mrSges(3,1) * t250 + (Ifges(3,3) + Ifges(4,1)) * t288 + (mrSges(5,1) * t243 + mrSges(6,1) * t236 - mrSges(6,2) * t238 - mrSges(5,3) * t241 + pkin(4) * t233 - t312 * qJD(4) - t316 * qJDD(4) + t323 * t266 - t318 * t267) * t292;
t234 = mrSges(6,2) * t267 + t264 * t314 - t304;
t1 = [mrSges(2,1) * t307 - mrSges(2,2) * t303 + Ifges(2,3) * qJDD(1) + t300 + pkin(1) * (t293 * (m(3) * t251 - mrSges(3,1) * t287 - mrSges(3,2) * t288 + t299) + t296 * (m(3) * t250 - mrSges(3,2) * t287 + (mrSges(3,1) - mrSges(4,2)) * t288 + t301)); t300; t228; mrSges(5,1) * t240 - mrSges(5,2) * t241 - mrSges(6,1) * t239 + mrSges(6,3) * t238 - pkin(4) * t234 + qJ(5) * t309 + t317 * t267 + (-mrSges(6,2) * qJ(5) - t316) * t266 + (Ifges(5,3) + Ifges(6,2)) * qJDD(4) + (-t313 * t295 + (-qJ(5) * t264 + t312) * t292) * t289; t234;];
tauJ = t1;
