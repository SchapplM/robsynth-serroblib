% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S4RPRR8
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
% tauJ [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S4RPRR8_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR8_invdynJ_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR8_invdynJ_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR8_invdynJ_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR8_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR8_invdynJ_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR8_invdynJ_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR8_invdynJ_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR8_invdynJ_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:55:09
% EndTime: 2019-12-31 16:55:10
% DurationCPUTime: 0.44s
% Computational Cost: add. (1430->142), mult. (2755->187), div. (0->0), fcn. (1486->6), ass. (0->59)
t252 = sin(qJ(1));
t255 = cos(qJ(1));
t263 = -g(1) * t255 - g(2) * t252;
t260 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t263;
t271 = -pkin(1) - pkin(5);
t256 = qJD(1) ^ 2;
t265 = g(1) * t252 - t255 * g(2);
t259 = -qJ(2) * t256 + qJDD(2) - t265;
t228 = t271 * qJDD(1) + t259;
t251 = sin(qJ(3));
t254 = cos(qJ(3));
t221 = t251 * g(3) + t254 * t228;
t268 = qJD(1) * qJD(3);
t266 = t251 * t268;
t238 = qJDD(1) * t254 - t266;
t206 = (-t238 - t266) * pkin(6) + (-t251 * t254 * t256 + qJDD(3)) * pkin(3) + t221;
t222 = -g(3) * t254 + t251 * t228;
t237 = -qJDD(1) * t251 - t254 * t268;
t269 = t254 * qJD(1);
t241 = qJD(3) * pkin(3) - pkin(6) * t269;
t249 = t251 ^ 2;
t207 = -pkin(3) * t249 * t256 + pkin(6) * t237 - qJD(3) * t241 + t222;
t250 = sin(qJ(4));
t253 = cos(qJ(4));
t204 = t206 * t253 - t207 * t250;
t234 = (-t254 * t250 - t251 * t253) * qJD(1);
t215 = qJD(4) * t234 + t237 * t250 + t238 * t253;
t235 = (-t251 * t250 + t254 * t253) * qJD(1);
t220 = -mrSges(5,1) * t234 + mrSges(5,2) * t235;
t247 = qJD(3) + qJD(4);
t225 = -mrSges(5,2) * t247 + mrSges(5,3) * t234;
t246 = qJDD(3) + qJDD(4);
t201 = m(5) * t204 + mrSges(5,1) * t246 - mrSges(5,3) * t215 - t220 * t235 + t225 * t247;
t205 = t206 * t250 + t207 * t253;
t214 = -qJD(4) * t235 + t237 * t253 - t238 * t250;
t226 = mrSges(5,1) * t247 - mrSges(5,3) * t235;
t202 = m(5) * t205 - mrSges(5,2) * t246 + mrSges(5,3) * t214 + t220 * t234 - t226 * t247;
t195 = t253 * t201 + t250 * t202;
t270 = qJD(1) * t251;
t264 = -t201 * t250 + t253 * t202;
t236 = (t251 * mrSges(4,1) + t254 * mrSges(4,2)) * qJD(1);
t239 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t270;
t240 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t269;
t262 = (m(4) * t221 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t238 + qJD(3) * t239 - t236 * t269 + t195) * t254 + (m(4) * t222 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t237 - qJD(3) * t240 - t236 * t270 + t264) * t251;
t210 = t241 * t269 - pkin(3) * t237 + (-pkin(6) * t249 + t271) * t256 + t260;
t258 = m(5) * t210 - mrSges(5,1) * t214 + t215 * mrSges(5,2) - t225 * t234 + t235 * t226;
t217 = Ifges(5,4) * t235 + Ifges(5,2) * t234 + Ifges(5,6) * t247;
t218 = Ifges(5,1) * t235 + Ifges(5,4) * t234 + Ifges(5,5) * t247;
t257 = mrSges(5,1) * t204 - mrSges(5,2) * t205 + Ifges(5,5) * t215 + Ifges(5,6) * t214 + Ifges(5,3) * t246 + t235 * t217 - t218 * t234;
t233 = (Ifges(4,5) * qJD(3)) + (t254 * Ifges(4,1) - t251 * Ifges(4,4)) * qJD(1);
t232 = (Ifges(4,6) * qJD(3)) + (t254 * Ifges(4,4) - t251 * Ifges(4,2)) * qJD(1);
t230 = -qJDD(1) * pkin(1) + t259;
t229 = pkin(1) * t256 - t260;
t227 = t271 * t256 + t260;
t216 = Ifges(5,5) * t235 + Ifges(5,6) * t234 + Ifges(5,3) * t247;
t197 = mrSges(5,2) * t210 - mrSges(5,3) * t204 + Ifges(5,1) * t215 + Ifges(5,4) * t214 + Ifges(5,5) * t246 + t216 * t234 - t217 * t247;
t196 = -mrSges(5,1) * t210 + mrSges(5,3) * t205 + Ifges(5,4) * t215 + Ifges(5,2) * t214 + Ifges(5,6) * t246 - t216 * t235 + t218 * t247;
t192 = m(3) * t230 + qJDD(1) * mrSges(3,2) - (mrSges(3,3) * t256) + t262;
t1 = [mrSges(2,1) * t265 - mrSges(2,2) * t263 + mrSges(3,2) * t230 - mrSges(3,3) * t229 + t254 * (mrSges(4,2) * t227 - mrSges(4,3) * t221 + Ifges(4,1) * t238 + Ifges(4,4) * t237 + Ifges(4,5) * qJDD(3) - pkin(6) * t195 - qJD(3) * t232 - t196 * t250 + t197 * t253) - t251 * (-mrSges(4,1) * t227 + mrSges(4,3) * t222 + Ifges(4,4) * t238 + Ifges(4,2) * t237 + Ifges(4,6) * qJDD(3) - pkin(3) * t258 + pkin(6) * t264 + qJD(3) * t233 + t253 * t196 + t250 * t197) - pkin(5) * t262 - pkin(1) * t192 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1) + (-m(3) * t229 + m(4) * t227 - mrSges(4,1) * t237 + mrSges(3,2) * t256 + mrSges(4,2) * t238 + t258 + qJDD(1) * mrSges(3,3) + (t239 * t251 + t240 * t254) * qJD(1)) * qJ(2); t192; mrSges(4,1) * t221 - mrSges(4,2) * t222 + Ifges(4,5) * t238 + Ifges(4,6) * t237 + Ifges(4,3) * qJDD(3) + pkin(3) * t195 + (t254 * t232 + t251 * t233) * qJD(1) + t257; t257;];
tauJ = t1;
