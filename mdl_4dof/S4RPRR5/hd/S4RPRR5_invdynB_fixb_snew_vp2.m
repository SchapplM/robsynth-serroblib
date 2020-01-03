% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S4RPRR5
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
%   pkin=[a2,a3,a4,d1,d3,d4]';
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
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S4RPRR5_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR5_invdynB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR5_invdynB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR5_invdynB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR5_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR5_invdynB_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR5_invdynB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR5_invdynB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR5_invdynB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:51:35
% EndTime: 2019-12-31 16:51:35
% DurationCPUTime: 0.38s
% Computational Cost: add. (3950->151), mult. (5170->182), div. (0->0), fcn. (1848->6), ass. (0->63)
t276 = -m(3) - m(4);
t275 = -pkin(1) - pkin(2);
t274 = -mrSges(2,1) - mrSges(3,1);
t273 = Ifges(3,4) + Ifges(2,5);
t272 = Ifges(2,6) - Ifges(3,6);
t248 = -qJD(1) + qJD(3);
t253 = sin(qJ(4));
t271 = t248 * t253;
t256 = cos(qJ(4));
t270 = t248 * t256;
t255 = sin(qJ(1));
t258 = cos(qJ(1));
t243 = -t258 * g(1) - t255 * g(2);
t259 = qJD(1) ^ 2;
t263 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t243;
t233 = -t259 * pkin(1) + t263;
t227 = t275 * t259 + t263;
t242 = t255 * g(1) - t258 * g(2);
t262 = -t259 * qJ(2) + qJDD(2) - t242;
t229 = t275 * qJDD(1) + t262;
t254 = sin(qJ(3));
t257 = cos(qJ(3));
t224 = t257 * t227 + t254 * t229;
t246 = t248 ^ 2;
t247 = -qJDD(1) + qJDD(3);
t222 = -t246 * pkin(3) + t247 * pkin(6) + t224;
t219 = t256 * g(3) - t253 * t222;
t237 = (-mrSges(5,1) * t256 + mrSges(5,2) * t253) * t248;
t268 = qJD(4) * t248;
t238 = t253 * t247 + t256 * t268;
t241 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t270;
t217 = m(5) * t219 + qJDD(4) * mrSges(5,1) - t238 * mrSges(5,3) + qJD(4) * t241 - t237 * t271;
t220 = t253 * g(3) + t256 * t222;
t239 = t256 * t247 - t253 * t268;
t240 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t271;
t218 = m(5) * t220 - qJDD(4) * mrSges(5,2) + t239 * mrSges(5,3) - qJD(4) * t240 + t237 * t270;
t265 = -t253 * t217 + t256 * t218;
t210 = m(4) * t224 - t246 * mrSges(4,1) - t247 * mrSges(4,2) + t265;
t223 = -t254 * t227 + t257 * t229;
t221 = -t247 * pkin(3) - t246 * pkin(6) - t223;
t260 = -m(5) * t221 + t239 * mrSges(5,1) - t238 * mrSges(5,2) - t240 * t271 + t241 * t270;
t215 = m(4) * t223 + t247 * mrSges(4,1) - t246 * mrSges(4,2) + t260;
t266 = t257 * t210 - t254 * t215;
t264 = m(3) * t233 + qJDD(1) * mrSges(3,3) + t266;
t205 = m(2) * t243 - qJDD(1) * mrSges(2,2) + t274 * t259 + t264;
t208 = t254 * t210 + t257 * t215;
t236 = -qJDD(1) * pkin(1) + t262;
t261 = -m(3) * t236 + qJDD(1) * mrSges(3,1) + t259 * mrSges(3,3) - t208;
t206 = m(2) * t242 + qJDD(1) * mrSges(2,1) - t259 * mrSges(2,2) + t261;
t269 = t255 * t205 + t258 * t206;
t267 = t258 * t205 - t255 * t206;
t212 = t256 * t217 + t253 * t218;
t232 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t253 + Ifges(5,4) * t256) * t248;
t231 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t253 + Ifges(5,2) * t256) * t248;
t230 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t253 + Ifges(5,6) * t256) * t248;
t214 = mrSges(5,2) * t221 - mrSges(5,3) * t219 + Ifges(5,1) * t238 + Ifges(5,4) * t239 + Ifges(5,5) * qJDD(4) - qJD(4) * t231 + t230 * t270;
t213 = -mrSges(5,1) * t221 + mrSges(5,3) * t220 + Ifges(5,4) * t238 + Ifges(5,2) * t239 + Ifges(5,6) * qJDD(4) + qJD(4) * t232 - t230 * t271;
t211 = t276 * g(3) - t212;
t207 = -mrSges(4,1) * g(3) - mrSges(5,1) * t219 + mrSges(5,2) * t220 + mrSges(4,3) * t224 + t246 * Ifges(4,5) - Ifges(5,5) * t238 + Ifges(4,6) * t247 - Ifges(5,6) * t239 - Ifges(5,3) * qJDD(4) - pkin(3) * t212 + (-t231 * t253 + t232 * t256) * t248;
t201 = mrSges(4,2) * g(3) - mrSges(4,3) * t223 + Ifges(4,5) * t247 - t246 * Ifges(4,6) - pkin(6) * t212 - t253 * t213 + t256 * t214;
t200 = mrSges(3,2) * t236 - mrSges(2,3) * t242 - pkin(5) * t208 - qJ(2) * t211 + t257 * t201 - t254 * t207 - t272 * t259 + t273 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3);
t199 = mrSges(2,3) * t243 + mrSges(3,2) * t233 - t254 * t201 - t257 * t207 + pkin(2) * t212 - pkin(5) * t266 - pkin(1) * t211 + t273 * t259 + t272 * qJDD(1) + (pkin(2) * m(4) - t274) * g(3);
t1 = [-m(1) * g(1) + t267; -m(1) * g(2) + t269; (-m(1) - m(2) + t276) * g(3) - t212; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(4) * t269 - t255 * t199 + t258 * t200; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(4) * t267 + t258 * t199 + t255 * t200; qJ(2) * (-t259 * mrSges(3,1) + t264) + pkin(1) * t261 + mrSges(2,1) * t242 - mrSges(2,2) * t243 - pkin(2) * t208 - mrSges(3,1) * t236 + mrSges(3,3) * t233 - t253 * t214 - t256 * t213 - pkin(3) * t260 - pkin(6) * t265 - mrSges(4,1) * t223 + mrSges(4,2) * t224 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) - Ifges(4,3) * t247 + (Ifges(3,2) + Ifges(2,3)) * qJDD(1);];
tauB = t1;
