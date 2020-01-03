% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S4PPRR3
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
%   pkin=[a2,a3,a4,d3,d4,theta1]';
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
% Datum: 2019-12-31 16:17
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S4PPRR3_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR3_invdynJB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR3_invdynJB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRR3_invdynJB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR3_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR3_invdynJB_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR3_invdynJB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPRR3_invdynJB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PPRR3_invdynJB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:17:24
% EndTime: 2019-12-31 16:17:25
% DurationCPUTime: 0.34s
% Computational Cost: add. (2263->124), mult. (3826->159), div. (0->0), fcn. (2068->6), ass. (0->55)
t262 = sin(pkin(6));
t263 = cos(pkin(6));
t255 = t262 * g(1) - t263 * g(2);
t253 = qJDD(2) - t255;
t256 = -t263 * g(1) - t262 * g(2);
t265 = sin(qJ(3));
t267 = cos(qJ(3));
t242 = t265 * t253 + t267 * t256;
t268 = qJD(3) ^ 2;
t240 = -t268 * pkin(3) + qJDD(3) * pkin(5) + t242;
t261 = g(3) - qJDD(1);
t264 = sin(qJ(4));
t266 = cos(qJ(4));
t237 = -t264 * t240 + t266 * t261;
t238 = t266 * t240 + t264 * t261;
t244 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t264 + Ifges(5,2) * t266) * qJD(3);
t245 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t264 + Ifges(5,4) * t266) * qJD(3);
t276 = qJD(3) * qJD(4);
t251 = t264 * qJDD(3) + t266 * t276;
t252 = t266 * qJDD(3) - t264 * t276;
t281 = mrSges(5,1) * t237 - mrSges(5,2) * t238 + Ifges(5,5) * t251 + Ifges(5,6) * t252 + Ifges(5,3) * qJDD(4) + (t244 * t264 - t245 * t266) * qJD(3);
t280 = -mrSges(2,2) + mrSges(3,3);
t250 = (-mrSges(5,1) * t266 + mrSges(5,2) * t264) * qJD(3);
t277 = qJD(3) * t266;
t258 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t277;
t278 = qJD(3) * t264;
t234 = m(5) * t237 + qJDD(4) * mrSges(5,1) - t251 * mrSges(5,3) + qJD(4) * t258 - t250 * t278;
t257 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t278;
t235 = m(5) * t238 - qJDD(4) * mrSges(5,2) + t252 * mrSges(5,3) - qJD(4) * t257 + t250 * t277;
t229 = -t264 * t234 + t266 * t235;
t226 = m(4) * t242 - t268 * mrSges(4,1) - qJDD(3) * mrSges(4,2) + t229;
t241 = t267 * t253 - t265 * t256;
t239 = -qJDD(3) * pkin(3) - t268 * pkin(5) - t241;
t236 = -m(5) * t239 + t252 * mrSges(5,1) - t251 * mrSges(5,2) - t257 * t278 + t258 * t277;
t232 = m(4) * t241 + qJDD(3) * mrSges(4,1) - t268 * mrSges(4,2) + t236;
t224 = t265 * t226 + t267 * t232;
t223 = m(3) * t253 + t224;
t221 = m(2) * t255 - t223;
t274 = t267 * t226 - t265 * t232;
t273 = m(3) * t256 + t274;
t222 = m(2) * t256 + t273;
t279 = t263 * t221 + t262 * t222;
t275 = -t262 * t221 + t263 * t222;
t228 = t266 * t234 + t264 * t235;
t227 = -t228 + (-m(3) - m(4)) * t261;
t271 = -m(2) * t261 + t227;
t243 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t264 + Ifges(5,6) * t266) * qJD(3);
t230 = -mrSges(5,1) * t239 + mrSges(5,3) * t238 + Ifges(5,4) * t251 + Ifges(5,2) * t252 + Ifges(5,6) * qJDD(4) + qJD(4) * t245 - t243 * t278;
t231 = mrSges(5,2) * t239 - mrSges(5,3) * t237 + Ifges(5,1) * t251 + Ifges(5,4) * t252 + Ifges(5,5) * qJDD(4) - qJD(4) * t244 + t243 * t277;
t269 = mrSges(4,1) * t241 - mrSges(4,2) * t242 + Ifges(4,3) * qJDD(3) + pkin(3) * t236 + pkin(5) * t229 + t266 * t230 + t264 * t231;
t217 = -mrSges(4,1) * t261 + mrSges(4,3) * t242 + t268 * Ifges(4,5) + Ifges(4,6) * qJDD(3) - pkin(3) * t228 - t281;
t216 = mrSges(4,2) * t261 - mrSges(4,3) * t241 + Ifges(4,5) * qJDD(3) - t268 * Ifges(4,6) - pkin(5) * t228 - t264 * t230 + t266 * t231;
t215 = mrSges(3,2) * t253 - mrSges(2,3) * t255 - pkin(4) * t224 - qJ(2) * t227 + t267 * t216 - t265 * t217 + t280 * t261;
t214 = -t265 * t216 - t267 * t217 + pkin(2) * t228 - pkin(4) * t274 - pkin(1) * t227 + (m(4) * pkin(2) + mrSges(2,1) + mrSges(3,1)) * t261 + (mrSges(3,2) + mrSges(2,3)) * t256;
t1 = [-m(1) * g(1) + t275; -m(1) * g(2) + t279; -m(1) * g(3) + t271; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t279 - t262 * t214 + t263 * t215; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t275 + t263 * t214 + t262 * t215; -mrSges(1,1) * g(2) + mrSges(2,1) * t255 - mrSges(3,1) * t253 + mrSges(1,2) * g(1) - pkin(1) * t223 - pkin(2) * t224 + qJ(2) * t273 + t280 * t256 - t269; t271; t223; t269; t281;];
tauJB = t1;
