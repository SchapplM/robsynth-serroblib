% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S4PRPR4
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
%   pkin=[a2,a3,a4,d2,d4,theta1]';
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
% Datum: 2019-12-31 16:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S4PRPR4_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR4_invdynB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR4_invdynB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPR4_invdynB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR4_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR4_invdynB_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR4_invdynB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPR4_invdynB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRPR4_invdynB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:21:58
% EndTime: 2019-12-31 16:21:58
% DurationCPUTime: 0.39s
% Computational Cost: add. (2511->135), mult. (4290->169), div. (0->0), fcn. (2134->6), ass. (0->59)
t259 = -pkin(2) - pkin(5);
t258 = mrSges(3,1) - mrSges(4,2);
t257 = -Ifges(4,4) + Ifges(3,5);
t256 = Ifges(4,5) - Ifges(3,6);
t236 = sin(pkin(6));
t237 = cos(pkin(6));
t226 = t236 * g(1) - t237 * g(2);
t227 = -t237 * g(1) - t236 * g(2);
t239 = sin(qJ(2));
t241 = cos(qJ(2));
t212 = t241 * t226 - t239 * t227;
t242 = qJD(2) ^ 2;
t246 = -t242 * qJ(3) + qJDD(3) - t212;
t209 = t259 * qJDD(2) + t246;
t233 = -g(3) + qJDD(1);
t238 = sin(qJ(4));
t240 = cos(qJ(4));
t205 = t240 * t209 - t238 * t233;
t223 = (mrSges(5,1) * t238 + mrSges(5,2) * t240) * qJD(2);
t252 = qJD(2) * qJD(4);
t225 = t240 * qJDD(2) - t238 * t252;
t254 = qJD(2) * t238;
t228 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t254;
t253 = qJD(2) * t240;
t203 = m(5) * t205 + qJDD(4) * mrSges(5,1) - t225 * mrSges(5,3) + qJD(4) * t228 - t223 * t253;
t206 = t238 * t209 + t240 * t233;
t224 = -t238 * qJDD(2) - t240 * t252;
t229 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t253;
t204 = m(5) * t206 - qJDD(4) * mrSges(5,2) + t224 * mrSges(5,3) - qJD(4) * t229 - t223 * t254;
t196 = t240 * t203 + t238 * t204;
t211 = -qJDD(2) * pkin(2) + t246;
t244 = -m(4) * t211 + t242 * mrSges(4,3) - t196;
t194 = m(3) * t212 - t242 * mrSges(3,2) + t258 * qJDD(2) + t244;
t213 = t239 * t226 + t241 * t227;
t245 = qJDD(2) * qJ(3) + 0.2e1 * qJD(3) * qJD(2) + t213;
t210 = t242 * pkin(2) - t245;
t208 = t259 * t242 + t245;
t247 = -m(5) * t208 + t224 * mrSges(5,1) - t225 * mrSges(5,2) - t228 * t254 - t229 * t253;
t243 = -m(4) * t210 + t242 * mrSges(4,2) + qJDD(2) * mrSges(4,3) - t247;
t199 = m(3) * t213 - t242 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t243;
t192 = t241 * t194 + t239 * t199;
t190 = m(2) * t226 + t192;
t250 = -t239 * t194 + t241 * t199;
t191 = m(2) * t227 + t250;
t255 = t237 * t190 + t236 * t191;
t251 = -t236 * t190 + t237 * t191;
t249 = -t238 * t203 + t240 * t204;
t195 = m(4) * t233 + t249;
t248 = m(3) * t233 + t195;
t216 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t240 - Ifges(5,4) * t238) * qJD(2);
t215 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t240 - Ifges(5,2) * t238) * qJD(2);
t214 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t240 - Ifges(5,6) * t238) * qJD(2);
t201 = mrSges(5,2) * t208 - mrSges(5,3) * t205 + Ifges(5,1) * t225 + Ifges(5,4) * t224 + Ifges(5,5) * qJDD(4) - qJD(4) * t215 - t214 * t254;
t200 = -mrSges(5,1) * t208 + mrSges(5,3) * t206 + Ifges(5,4) * t225 + Ifges(5,2) * t224 + Ifges(5,6) * qJDD(4) + qJD(4) * t216 - t214 * t253;
t186 = mrSges(4,1) * t211 + mrSges(5,1) * t205 - mrSges(5,2) * t206 - mrSges(3,3) * t212 + Ifges(5,5) * t225 + Ifges(5,6) * t224 + Ifges(5,3) * qJDD(4) + pkin(3) * t196 - qJ(3) * t195 + t256 * t242 + (mrSges(3,2) - mrSges(4,3)) * t233 + t257 * qJDD(2) + (t240 * t215 + t238 * t216) * qJD(2);
t185 = -mrSges(4,1) * t210 + mrSges(3,3) * t213 - pkin(2) * t195 - pkin(3) * t247 - pkin(5) * t249 - t256 * qJDD(2) - t240 * t200 - t238 * t201 - t258 * t233 + t257 * t242;
t184 = mrSges(2,2) * t233 - mrSges(2,3) * t226 - pkin(4) * t192 - t239 * t185 + t241 * t186;
t183 = -mrSges(2,1) * t233 + mrSges(2,3) * t227 - pkin(1) * t248 + pkin(4) * t250 + t241 * t185 + t239 * t186;
t1 = [-m(1) * g(1) + t251; -m(1) * g(2) + t255; -m(1) * g(3) + m(2) * t233 + t248; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t255 - t236 * t183 + t237 * t184; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t251 + t237 * t183 + t236 * t184; pkin(1) * t192 + mrSges(2,1) * t226 - mrSges(2,2) * t227 + pkin(2) * t244 + qJ(3) * t243 - t238 * t200 - pkin(5) * t196 + mrSges(3,1) * t212 - mrSges(3,2) * t213 + mrSges(4,2) * t211 - mrSges(4,3) * t210 + t240 * t201 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (-pkin(2) * mrSges(4,2) + Ifges(4,1) + Ifges(3,3)) * qJDD(2);];
tauB = t1;
