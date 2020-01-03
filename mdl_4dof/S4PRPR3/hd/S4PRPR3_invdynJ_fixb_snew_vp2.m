% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S4PRPR3
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
% Datum: 2019-12-31 16:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S4PRPR3_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR3_invdynJ_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR3_invdynJ_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPR3_invdynJ_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR3_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR3_invdynJ_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR3_invdynJ_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPR3_invdynJ_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRPR3_invdynJ_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:20:51
% EndTime: 2019-12-31 16:20:52
% DurationCPUTime: 0.27s
% Computational Cost: add. (890->102), mult. (1935->142), div. (0->0), fcn. (1227->8), ass. (0->56)
t239 = qJD(2) ^ 2;
t233 = cos(pkin(7));
t229 = t233 ^ 2;
t255 = 0.2e1 * t233;
t254 = pkin(3) * t239;
t253 = pkin(5) * qJDD(2);
t232 = sin(pkin(6));
t234 = cos(pkin(6));
t222 = g(1) * t232 - g(2) * t234;
t223 = -g(1) * t234 - g(2) * t232;
t236 = sin(qJ(2));
t238 = cos(qJ(2));
t251 = t236 * t222 + t238 * t223;
t212 = -pkin(2) * t239 + qJDD(2) * qJ(3) + t251;
t231 = sin(pkin(7));
t230 = -g(3) + qJDD(1);
t248 = qJD(2) * qJD(3);
t250 = t233 * t230 - 0.2e1 * t231 * t248;
t199 = (t233 * t254 - t212 - t253) * t231 + t250;
t203 = t233 * t212 + t231 * t230 + t248 * t255;
t200 = -t229 * t254 + t233 * t253 + t203;
t235 = sin(qJ(4));
t237 = cos(qJ(4));
t197 = t199 * t237 - t200 * t235;
t242 = -t231 * t235 + t233 * t237;
t215 = t242 * qJD(2);
t243 = t231 * t237 + t233 * t235;
t216 = t243 * qJD(2);
t208 = -mrSges(5,1) * t215 + mrSges(5,2) * t216;
t211 = t215 * qJD(4) + t243 * qJDD(2);
t213 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t215;
t195 = m(5) * t197 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t211 + qJD(4) * t213 - t208 * t216;
t198 = t199 * t235 + t200 * t237;
t210 = -t216 * qJD(4) + t242 * qJDD(2);
t214 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t216;
t196 = m(5) * t198 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t210 - qJD(4) * t214 + t208 * t215;
t252 = t237 * t195 + t235 * t196;
t249 = -t231 ^ 2 - t229;
t247 = -t235 * t195 + t237 * t196;
t246 = t238 * t222 - t236 * t223;
t245 = -mrSges(4,1) * t233 + mrSges(4,2) * t231;
t244 = qJDD(3) - t246;
t241 = mrSges(4,3) * qJDD(2) + t239 * t245;
t201 = (-pkin(3) * t233 - pkin(2)) * qJDD(2) + (t249 * pkin(5) - qJ(3)) * t239 + t244;
t240 = m(5) * t201 - t210 * mrSges(5,1) + t211 * mrSges(5,2) - t215 * t213 + t216 * t214;
t209 = -qJDD(2) * pkin(2) - t239 * qJ(3) + t244;
t206 = Ifges(5,1) * t216 + Ifges(5,4) * t215 + Ifges(5,5) * qJD(4);
t205 = Ifges(5,4) * t216 + Ifges(5,2) * t215 + Ifges(5,6) * qJD(4);
t204 = Ifges(5,5) * t216 + Ifges(5,6) * t215 + Ifges(5,3) * qJD(4);
t202 = -t212 * t231 + t250;
t191 = t249 * t239 * mrSges(4,3) + m(4) * t209 + t245 * qJDD(2) + t240;
t190 = mrSges(5,2) * t201 - mrSges(5,3) * t197 + Ifges(5,1) * t211 + Ifges(5,4) * t210 + Ifges(5,5) * qJDD(4) - qJD(4) * t205 + t204 * t215;
t189 = -mrSges(5,1) * t201 + mrSges(5,3) * t198 + Ifges(5,4) * t211 + Ifges(5,2) * t210 + Ifges(5,6) * qJDD(4) + qJD(4) * t206 - t204 * t216;
t188 = m(4) * t203 + t241 * t233 + t247;
t187 = m(4) * t202 - t241 * t231 + t252;
t1 = [t233 * t187 + t231 * t188 + (m(2) + m(3)) * t230; mrSges(3,1) * t246 - mrSges(3,2) * t251 + t231 * (mrSges(4,2) * t209 - mrSges(4,3) * t202 - pkin(5) * t252 - t235 * t189 + t237 * t190) + t233 * (-mrSges(4,1) * t209 + mrSges(4,3) * t203 - pkin(3) * t240 + pkin(5) * t247 + t237 * t189 + t235 * t190) - pkin(2) * t191 + qJ(3) * (-t187 * t231 + t188 * t233) + (Ifges(4,2) * t229 + Ifges(3,3) + (Ifges(4,1) * t231 + Ifges(4,4) * t255) * t231) * qJDD(2); t191; mrSges(5,1) * t197 - mrSges(5,2) * t198 + Ifges(5,5) * t211 + Ifges(5,6) * t210 + Ifges(5,3) * qJDD(4) + t205 * t216 - t206 * t215;];
tauJ = t1;
