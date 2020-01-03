% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S4PRPR6
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
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
% tauJ [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:24
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S4PRPR6_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR6_invdynJ_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR6_invdynJ_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPR6_invdynJ_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR6_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR6_invdynJ_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR6_invdynJ_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPR6_invdynJ_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRPR6_invdynJ_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:24:23
% EndTime: 2019-12-31 16:24:24
% DurationCPUTime: 0.29s
% Computational Cost: add. (937->110), mult. (2006->150), div. (0->0), fcn. (1258->8), ass. (0->57)
t241 = qJD(2) ^ 2;
t235 = cos(pkin(7));
t231 = t235 ^ 2;
t257 = 0.2e1 * t235;
t256 = pkin(3) * t235;
t233 = sin(pkin(7));
t255 = mrSges(4,2) * t233;
t254 = t231 * t241;
t234 = sin(pkin(6));
t236 = cos(pkin(6));
t223 = -t236 * g(1) - t234 * g(2);
t232 = -g(3) + qJDD(1);
t238 = sin(qJ(2));
t240 = cos(qJ(2));
t213 = t240 * t223 + t238 * t232;
t209 = -t241 * pkin(2) + qJDD(2) * qJ(3) + t213;
t222 = -t234 * g(1) + t236 * g(2);
t251 = qJD(2) * qJD(3);
t252 = t235 * t222 - 0.2e1 * t233 * t251;
t194 = (-pkin(5) * qJDD(2) + t241 * t256 - t209) * t233 + t252;
t197 = t235 * t209 + t233 * t222 + t251 * t257;
t250 = t235 * qJDD(2);
t195 = -pkin(3) * t254 + pkin(5) * t250 + t197;
t237 = sin(qJ(4));
t239 = cos(qJ(4));
t192 = t239 * t194 - t237 * t195;
t245 = -t233 * t237 + t235 * t239;
t214 = t245 * qJD(2);
t246 = t233 * t239 + t235 * t237;
t215 = t246 * qJD(2);
t203 = -t214 * mrSges(5,1) + t215 * mrSges(5,2);
t206 = t214 * qJD(4) + qJDD(2) * t246;
t210 = -qJD(4) * mrSges(5,2) + t214 * mrSges(5,3);
t189 = m(5) * t192 + qJDD(4) * mrSges(5,1) - t206 * mrSges(5,3) + qJD(4) * t210 - t215 * t203;
t193 = t237 * t194 + t239 * t195;
t205 = -t215 * qJD(4) + qJDD(2) * t245;
t211 = qJD(4) * mrSges(5,1) - t215 * mrSges(5,3);
t190 = m(5) * t193 - qJDD(4) * mrSges(5,2) + t205 * mrSges(5,3) - qJD(4) * t211 + t214 * t203;
t253 = t239 * t189 + t237 * t190;
t196 = -t233 * t209 + t252;
t244 = mrSges(4,3) * qJDD(2) + t241 * (-t235 * mrSges(4,1) + t255);
t248 = -t237 * t189 + t239 * t190;
t249 = -t233 * (m(4) * t196 - t233 * t244 + t253) + t235 * (m(4) * t197 + t235 * t244 + t248);
t212 = -t238 * t223 + t240 * t232;
t247 = qJDD(3) - t212;
t230 = t233 ^ 2;
t198 = (-pkin(2) - t256) * qJDD(2) + (-qJ(3) + (-t230 - t231) * pkin(5)) * t241 + t247;
t243 = m(5) * t198 - t205 * mrSges(5,1) + t206 * mrSges(5,2) - t214 * t210 + t215 * t211;
t208 = -qJDD(2) * pkin(2) - t241 * qJ(3) + t247;
t242 = -m(4) * t208 + mrSges(4,1) * t250 - t243 + (t230 * t241 + t254) * mrSges(4,3);
t201 = Ifges(5,1) * t215 + Ifges(5,4) * t214 + Ifges(5,5) * qJD(4);
t200 = Ifges(5,4) * t215 + Ifges(5,2) * t214 + Ifges(5,6) * qJD(4);
t199 = Ifges(5,5) * t215 + Ifges(5,6) * t214 + Ifges(5,3) * qJD(4);
t191 = qJDD(2) * t255 - t242;
t185 = mrSges(5,2) * t198 - mrSges(5,3) * t192 + Ifges(5,1) * t206 + Ifges(5,4) * t205 + Ifges(5,5) * qJDD(4) - qJD(4) * t200 + t214 * t199;
t184 = -mrSges(5,1) * t198 + mrSges(5,3) * t193 + Ifges(5,4) * t206 + Ifges(5,2) * t205 + Ifges(5,6) * qJDD(4) + qJD(4) * t201 - t215 * t199;
t1 = [m(2) * t232 + t238 * (m(3) * t213 - t241 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t249) + t240 * (m(3) * t212 - t241 * mrSges(3,2) + (mrSges(3,1) - t255) * qJDD(2) + t242); mrSges(3,1) * t212 - mrSges(3,2) * t213 + t233 * (mrSges(4,2) * t208 - mrSges(4,3) * t196 - pkin(5) * t253 - t237 * t184 + t239 * t185) + t235 * (-mrSges(4,1) * t208 + mrSges(4,3) * t197 - pkin(3) * t243 + pkin(5) * t248 + t239 * t184 + t237 * t185) - pkin(2) * t191 + qJ(3) * t249 + (Ifges(4,2) * t231 + Ifges(3,3) + (Ifges(4,1) * t233 + Ifges(4,4) * t257) * t233) * qJDD(2); t191; mrSges(5,1) * t192 - mrSges(5,2) * t193 + Ifges(5,5) * t206 + Ifges(5,6) * t205 + Ifges(5,3) * qJDD(4) + t215 * t200 - t214 * t201;];
tauJ = t1;
