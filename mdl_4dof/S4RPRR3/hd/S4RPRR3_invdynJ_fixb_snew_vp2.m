% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S4RPRR3
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
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 16:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S4RPRR3_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR3_invdynJ_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR3_invdynJ_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR3_invdynJ_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR3_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR3_invdynJ_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR3_invdynJ_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR3_invdynJ_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR3_invdynJ_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:49:09
% EndTime: 2019-12-31 16:49:10
% DurationCPUTime: 0.39s
% Computational Cost: add. (1684->145), mult. (3267->194), div. (0->0), fcn. (1847->8), ass. (0->64)
t254 = sin(qJ(1));
t257 = cos(qJ(1));
t266 = t254 * g(1) - g(2) * t257;
t234 = qJDD(1) * pkin(1) + t266;
t258 = qJD(1) ^ 2;
t263 = -g(1) * t257 - g(2) * t254;
t236 = -pkin(1) * t258 + t263;
t250 = sin(pkin(7));
t251 = cos(pkin(7));
t221 = t250 * t234 + t251 * t236;
t218 = -pkin(2) * t258 + qJDD(1) * pkin(5) + t221;
t249 = -g(3) + qJDD(2);
t253 = sin(qJ(3));
t256 = cos(qJ(3));
t208 = -t218 * t253 + t256 * t249;
t268 = qJD(1) * qJD(3);
t267 = t256 * t268;
t237 = qJDD(1) * t253 + t267;
t201 = (-t237 + t267) * pkin(6) + (t253 * t256 * t258 + qJDD(3)) * pkin(3) + t208;
t209 = t256 * t218 + t253 * t249;
t238 = qJDD(1) * t256 - t253 * t268;
t270 = qJD(1) * t253;
t241 = qJD(3) * pkin(3) - pkin(6) * t270;
t248 = t256 ^ 2;
t202 = -pkin(3) * t248 * t258 + pkin(6) * t238 - qJD(3) * t241 + t209;
t252 = sin(qJ(4));
t255 = cos(qJ(4));
t199 = t201 * t255 - t202 * t252;
t230 = (-t253 * t252 + t256 * t255) * qJD(1);
t211 = qJD(4) * t230 + t237 * t255 + t238 * t252;
t231 = (t256 * t252 + t253 * t255) * qJD(1);
t219 = -mrSges(5,1) * t230 + mrSges(5,2) * t231;
t247 = qJD(3) + qJD(4);
t222 = -mrSges(5,2) * t247 + mrSges(5,3) * t230;
t246 = qJDD(3) + qJDD(4);
t196 = m(5) * t199 + mrSges(5,1) * t246 - mrSges(5,3) * t211 - t219 * t231 + t222 * t247;
t200 = t201 * t252 + t202 * t255;
t210 = -qJD(4) * t231 - t237 * t252 + t238 * t255;
t223 = mrSges(5,1) * t247 - mrSges(5,3) * t231;
t197 = m(5) * t200 - mrSges(5,2) * t246 + mrSges(5,3) * t210 + t219 * t230 - t223 * t247;
t190 = t255 * t196 + t252 * t197;
t269 = qJD(1) * t256;
t235 = (-t256 * mrSges(4,1) + t253 * mrSges(4,2)) * qJD(1);
t240 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t269;
t188 = m(4) * t208 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t237 + qJD(3) * t240 - t235 * t270 + t190;
t239 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t270;
t264 = -t196 * t252 + t255 * t197;
t189 = m(4) * t209 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t238 - qJD(3) * t239 + t235 * t269 + t264;
t265 = -t188 * t253 + t256 * t189;
t220 = t234 * t251 - t250 * t236;
t262 = -qJDD(1) * pkin(2) - t220;
t203 = t241 * t270 - pkin(3) * t238 + (-pkin(6) * t248 - pkin(5)) * t258 + t262;
t261 = m(5) * t203 - t210 * mrSges(5,1) + mrSges(5,2) * t211 - t230 * t222 + t223 * t231;
t214 = Ifges(5,4) * t231 + Ifges(5,2) * t230 + Ifges(5,6) * t247;
t215 = Ifges(5,1) * t231 + Ifges(5,4) * t230 + Ifges(5,5) * t247;
t260 = mrSges(5,1) * t199 - mrSges(5,2) * t200 + Ifges(5,5) * t211 + Ifges(5,6) * t210 + Ifges(5,3) * t246 + t231 * t214 - t215 * t230;
t217 = -pkin(5) * t258 + t262;
t259 = -m(4) * t217 + t238 * mrSges(4,1) - mrSges(4,2) * t237 - t239 * t270 + t240 * t269 - t261;
t229 = Ifges(4,5) * qJD(3) + (t253 * Ifges(4,1) + t256 * Ifges(4,4)) * qJD(1);
t228 = Ifges(4,6) * qJD(3) + (t253 * Ifges(4,4) + t256 * Ifges(4,2)) * qJD(1);
t213 = Ifges(5,5) * t231 + Ifges(5,6) * t230 + Ifges(5,3) * t247;
t192 = mrSges(5,2) * t203 - mrSges(5,3) * t199 + Ifges(5,1) * t211 + Ifges(5,4) * t210 + Ifges(5,5) * t246 + t213 * t230 - t214 * t247;
t191 = -mrSges(5,1) * t203 + mrSges(5,3) * t200 + Ifges(5,4) * t211 + Ifges(5,2) * t210 + Ifges(5,6) * t246 - t213 * t231 + t215 * t247;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t266 - mrSges(2,2) * t263 + Ifges(3,3) * qJDD(1) + mrSges(3,1) * t220 - mrSges(3,2) * t221 + t253 * (mrSges(4,2) * t217 - mrSges(4,3) * t208 + Ifges(4,1) * t237 + Ifges(4,4) * t238 + Ifges(4,5) * qJDD(3) - pkin(6) * t190 - qJD(3) * t228 - t252 * t191 + t255 * t192) + t256 * (-mrSges(4,1) * t217 + mrSges(4,3) * t209 + Ifges(4,4) * t237 + Ifges(4,2) * t238 + Ifges(4,6) * qJDD(3) - pkin(3) * t261 + pkin(6) * t264 + qJD(3) * t229 + t255 * t191 + t252 * t192) + pkin(2) * t259 + pkin(5) * t265 + pkin(1) * (t250 * (m(3) * t221 - mrSges(3,1) * t258 - qJDD(1) * mrSges(3,2) + t265) + t251 * (m(3) * t220 + qJDD(1) * mrSges(3,1) - mrSges(3,2) * t258 + t259)); m(3) * t249 + t188 * t256 + t189 * t253; mrSges(4,1) * t208 - mrSges(4,2) * t209 + Ifges(4,5) * t237 + Ifges(4,6) * t238 + Ifges(4,3) * qJDD(3) + pkin(3) * t190 + (t253 * t228 - t256 * t229) * qJD(1) + t260; t260;];
tauJ = t1;
