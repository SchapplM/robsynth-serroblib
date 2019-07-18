% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S4PRRR2
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4]';
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
% Datum: 2019-07-18 13:27
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S4PRRR2_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(2,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR2_invdynJB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR2_invdynJB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRR2_invdynJB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR2_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S4PRRR2_invdynJB_fixb_snew_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR2_invdynJB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRR2_invdynJB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRR2_invdynJB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:27:23
% EndTime: 2019-07-18 13:27:23
% DurationCPUTime: 0.27s
% Computational Cost: add. (1983->100), mult. (2561->114), div. (0->0), fcn. (1318->6), ass. (0->48)
t254 = g(2) + qJDD(1);
t268 = (m(4) + m(5)) * t254;
t272 = (m(2) + m(3)) * t254 + t268;
t271 = -m(1) - m(2);
t257 = sin(qJ(2));
t260 = cos(qJ(2));
t239 = t260 * g(1) + t257 * g(3);
t237 = qJDD(2) * pkin(1) + t239;
t240 = t257 * g(1) - t260 * g(3);
t261 = qJD(2) ^ 2;
t238 = -t261 * pkin(1) + t240;
t256 = sin(qJ(3));
t259 = cos(qJ(3));
t232 = t259 * t237 - t256 * t238;
t253 = qJD(2) + qJD(3);
t251 = t253 ^ 2;
t252 = qJDD(2) + qJDD(3);
t229 = t252 * pkin(2) + t232;
t233 = t256 * t237 + t259 * t238;
t230 = -t251 * pkin(2) + t233;
t255 = sin(qJ(4));
t258 = cos(qJ(4));
t227 = t258 * t229 - t255 * t230;
t245 = qJD(4) + t253;
t243 = t245 ^ 2;
t244 = qJDD(4) + t252;
t224 = m(5) * t227 + t244 * mrSges(5,1) - t243 * mrSges(5,2);
t228 = t255 * t229 + t258 * t230;
t225 = m(5) * t228 - t243 * mrSges(5,1) - t244 * mrSges(5,2);
t269 = t258 * t224 + t255 * t225;
t216 = m(4) * t232 + t252 * mrSges(4,1) - t251 * mrSges(4,2) + t269;
t217 = m(4) * t233 - t251 * mrSges(4,1) - t252 * mrSges(4,2) - t255 * t224 + t258 * t225;
t211 = t259 * t216 + t256 * t217;
t209 = m(3) * t239 + qJDD(2) * mrSges(3,1) - t261 * mrSges(3,2) + t211;
t210 = m(3) * t240 - t261 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t256 * t216 + t259 * t217;
t267 = -t257 * t209 + t260 * t210;
t266 = qJ(1) * m(2) - mrSges(1,2) + mrSges(2,3);
t265 = -mrSges(5,1) * t227 + mrSges(5,2) * t228 - Ifges(5,3) * t244;
t264 = -t260 * t209 - t257 * t210;
t263 = -mrSges(4,1) * t232 + mrSges(4,2) * t233 - Ifges(4,3) * t252 - pkin(2) * t269 + t265;
t262 = mrSges(3,1) * t239 - mrSges(3,2) * t240 + Ifges(3,3) * qJDD(2) + pkin(1) * t211 - t263;
t221 = mrSges(5,2) * t254 - mrSges(5,3) * t227 + Ifges(5,5) * t244 - t243 * Ifges(5,6);
t220 = -mrSges(5,1) * t254 + mrSges(5,3) * t228 + t243 * Ifges(5,5) + Ifges(5,6) * t244;
t213 = mrSges(4,2) * t254 - mrSges(4,3) * t232 + Ifges(4,5) * t252 - t251 * Ifges(4,6) - t255 * t220 + t258 * t221;
t212 = mrSges(4,3) * t233 + t251 * Ifges(4,5) + Ifges(4,6) * t252 + t258 * t220 + t255 * t221 + (-m(5) * pkin(2) - mrSges(4,1)) * t254;
t207 = mrSges(3,2) * t254 - mrSges(3,3) * t239 + Ifges(3,5) * qJDD(2) - t261 * Ifges(3,6) - t256 * t212 + t259 * t213;
t206 = -mrSges(3,1) * t254 + mrSges(3,3) * t240 + t261 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(1) * t268 + t259 * t212 + t256 * t213;
t1 = [t271 * g(1) + t264; -m(1) * g(2) - t272; t271 * g(3) + t267; mrSges(2,1) * t254 + mrSges(1,3) * g(2) + t266 * g(3) - qJ(1) * t267 - t260 * t206 - t257 * t207; -t262 + (mrSges(1,1) - mrSges(2,2)) * g(3) + (-mrSges(2,1) - mrSges(1,3)) * g(1); -mrSges(1,1) * g(2) + mrSges(2,2) * t254 - t266 * g(1) + qJ(1) * t264 - t257 * t206 + t260 * t207; t272; t262; -t263; -t265;];
tauJB  = t1;
