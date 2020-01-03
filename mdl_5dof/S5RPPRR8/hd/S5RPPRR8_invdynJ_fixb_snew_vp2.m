% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RPPRR8
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
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
% Datum: 2019-12-31 18:01
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RPPRR8_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR8_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR8_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR8_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR8_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR8_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR8_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR8_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR8_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:01:09
% EndTime: 2019-12-31 18:01:09
% DurationCPUTime: 0.31s
% Computational Cost: add. (2820->112), mult. (3920->144), div. (0->0), fcn. (1404->8), ass. (0->56)
t261 = -pkin(1) - pkin(2);
t239 = -qJD(1) + qJD(4);
t244 = sin(qJ(5));
t260 = t239 * t244;
t247 = cos(qJ(5));
t259 = t239 * t247;
t250 = qJD(1) ^ 2;
t246 = sin(qJ(1));
t249 = cos(qJ(1));
t255 = -g(1) * t249 - g(2) * t246;
t253 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t255;
t221 = t261 * t250 + t253;
t256 = g(1) * t246 - t249 * g(2);
t252 = -qJ(2) * t250 + qJDD(2) - t256;
t222 = t261 * qJDD(1) + t252;
t242 = sin(pkin(8));
t243 = cos(pkin(8));
t216 = -t221 * t242 + t243 * t222;
t214 = -qJDD(1) * pkin(3) + t216;
t217 = t243 * t221 + t242 * t222;
t215 = -pkin(3) * t250 + t217;
t245 = sin(qJ(4));
t248 = cos(qJ(4));
t211 = t245 * t214 + t248 * t215;
t237 = t239 ^ 2;
t238 = -qJDD(1) + qJDD(4);
t209 = -pkin(4) * t237 + pkin(7) * t238 + t211;
t241 = g(3) + qJDD(3);
t206 = -t209 * t244 + t241 * t247;
t230 = (-mrSges(6,1) * t247 + mrSges(6,2) * t244) * t239;
t257 = qJD(5) * t239;
t231 = t238 * t244 + t247 * t257;
t234 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t259;
t204 = m(6) * t206 + qJDD(5) * mrSges(6,1) - mrSges(6,3) * t231 + qJD(5) * t234 - t230 * t260;
t207 = t209 * t247 + t241 * t244;
t232 = t238 * t247 - t244 * t257;
t233 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t260;
t205 = m(6) * t207 - qJDD(5) * mrSges(6,2) + mrSges(6,3) * t232 - qJD(5) * t233 + t230 * t259;
t197 = -t204 * t244 + t247 * t205;
t196 = m(5) * t211 - mrSges(5,1) * t237 - mrSges(5,2) * t238 + t197;
t210 = t214 * t248 - t215 * t245;
t208 = -pkin(4) * t238 - pkin(7) * t237 - t210;
t202 = -m(6) * t208 + t232 * mrSges(6,1) - mrSges(6,2) * t231 - t233 * t260 + t234 * t259;
t201 = m(5) * t210 + mrSges(5,1) * t238 - mrSges(5,2) * t237 + t202;
t258 = t245 * t196 + t248 * t201;
t193 = m(4) * t216 - qJDD(1) * mrSges(4,1) - mrSges(4,2) * t250 + t258;
t194 = m(4) * t217 - mrSges(4,1) * t250 + qJDD(1) * mrSges(4,2) + t196 * t248 - t201 * t245;
t254 = t193 * t243 + t194 * t242;
t223 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t244 + Ifges(6,6) * t247) * t239;
t224 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t244 + Ifges(6,2) * t247) * t239;
t225 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t244 + Ifges(6,4) * t247) * t239;
t251 = mrSges(5,1) * t210 - mrSges(5,2) * t211 + Ifges(5,3) * t238 + pkin(4) * t202 + pkin(7) * t197 + t247 * (-mrSges(6,1) * t208 + mrSges(6,3) * t207 + Ifges(6,4) * t231 + Ifges(6,2) * t232 + Ifges(6,6) * qJDD(5) + qJD(5) * t225 - t223 * t260) + t244 * (mrSges(6,2) * t208 - mrSges(6,3) * t206 + Ifges(6,1) * t231 + Ifges(6,4) * t232 + Ifges(6,5) * qJDD(5) - qJD(5) * t224 + t223 * t259);
t229 = -qJDD(1) * pkin(1) + t252;
t226 = -pkin(1) * t250 + t253;
t192 = m(3) * t229 - qJDD(1) * mrSges(3,1) - mrSges(3,3) * t250 + t254;
t1 = [-pkin(1) * t192 - pkin(3) * t258 - mrSges(4,1) * t216 + mrSges(4,2) * t217 + mrSges(3,3) * t226 - mrSges(3,1) * t229 - t251 - mrSges(2,2) * t255 + qJ(2) * (m(3) * t226 - mrSges(3,1) * t250 - t193 * t242 + t194 * t243) - pkin(2) * t254 + mrSges(2,1) * t256 + (mrSges(3,3) * qJ(2) + Ifges(3,2) + Ifges(2,3) + Ifges(4,3)) * qJDD(1); t192; t247 * t204 + t244 * t205 + (m(4) + m(5)) * t241; t251; mrSges(6,1) * t206 - mrSges(6,2) * t207 + Ifges(6,5) * t231 + Ifges(6,6) * t232 + Ifges(6,3) * qJDD(5) + (t224 * t244 - t225 * t247) * t239;];
tauJ = t1;
