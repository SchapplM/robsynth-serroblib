% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRPPR4
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
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
% Datum: 2019-12-31 19:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRPPR4_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR4_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR4_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR4_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR4_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR4_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR4_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR4_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR4_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:27:45
% EndTime: 2019-12-31 19:27:45
% DurationCPUTime: 0.37s
% Computational Cost: add. (2922->111), mult. (3375->142), div. (0->0), fcn. (1430->8), ass. (0->53)
t244 = qJD(1) + qJD(2);
t265 = t244 ^ 2;
t264 = -pkin(2) - pkin(3);
t249 = sin(qJ(5));
t263 = t244 * t249;
t252 = cos(qJ(5));
t262 = t244 * t252;
t251 = sin(qJ(1));
t254 = cos(qJ(1));
t260 = t251 * g(1) - t254 * g(2);
t229 = qJDD(1) * pkin(1) + t260;
t258 = -t254 * g(1) - t251 * g(2);
t230 = -qJD(1) ^ 2 * pkin(1) + t258;
t250 = sin(qJ(2));
t253 = cos(qJ(2));
t217 = t250 * t229 + t253 * t230;
t243 = qJDD(1) + qJDD(2);
t259 = t243 * qJ(3) + 0.2e1 * qJD(3) * t244 + t217;
t208 = t264 * t265 + t259;
t216 = t253 * t229 - t250 * t230;
t256 = -qJ(3) * t265 + qJDD(3) - t216;
t212 = t264 * t243 + t256;
t247 = sin(pkin(8));
t248 = cos(pkin(8));
t206 = t248 * t208 + t247 * t212;
t261 = qJD(5) * t244;
t203 = -pkin(4) * t265 - t243 * pkin(7) + t206;
t246 = g(3) + qJDD(4);
t200 = -t249 * t203 + t252 * t246;
t223 = (mrSges(6,1) * t252 - mrSges(6,2) * t249) * t244;
t224 = -t249 * t243 - t252 * t261;
t232 = -qJD(5) * mrSges(6,2) - mrSges(6,3) * t262;
t198 = m(6) * t200 + qJDD(5) * mrSges(6,1) - t224 * mrSges(6,3) + qJD(5) * t232 + t223 * t263;
t201 = t252 * t203 + t249 * t246;
t225 = -t252 * t243 + t249 * t261;
t231 = qJD(5) * mrSges(6,1) + mrSges(6,3) * t263;
t199 = m(6) * t201 - qJDD(5) * mrSges(6,2) + t225 * mrSges(6,3) - qJD(5) * t231 - t223 * t262;
t192 = -t249 * t198 + t252 * t199;
t191 = m(5) * t206 - mrSges(5,1) * t265 + t243 * mrSges(5,2) + t192;
t205 = -t247 * t208 + t248 * t212;
t202 = t243 * pkin(4) - pkin(7) * t265 - t205;
t196 = -m(6) * t202 + t225 * mrSges(6,1) - t224 * mrSges(6,2) + t231 * t263 - t232 * t262;
t195 = m(5) * t205 - t243 * mrSges(5,1) - mrSges(5,2) * t265 + t196;
t189 = t247 * t191 + t248 * t195;
t213 = -pkin(2) * t265 + t259;
t257 = m(4) * t213 + t243 * mrSges(4,3) + t248 * t191 - t247 * t195;
t214 = -t243 * pkin(2) + t256;
t188 = m(4) * t214 - t243 * mrSges(4,1) - mrSges(4,3) * t265 + t189;
t218 = Ifges(6,3) * qJD(5) + (-Ifges(6,5) * t249 - Ifges(6,6) * t252) * t244;
t219 = Ifges(6,6) * qJD(5) + (-Ifges(6,4) * t249 - Ifges(6,2) * t252) * t244;
t220 = Ifges(6,5) * qJD(5) + (-Ifges(6,1) * t249 - Ifges(6,4) * t252) * t244;
t255 = -mrSges(4,1) * t214 - mrSges(5,1) * t205 - mrSges(3,2) * t217 - pkin(3) * t189 - pkin(4) * t196 - pkin(7) * t192 - t252 * (-mrSges(6,1) * t202 + mrSges(6,3) * t201 + Ifges(6,4) * t224 + Ifges(6,2) * t225 + Ifges(6,6) * qJDD(5) + qJD(5) * t220 + t218 * t263) - t249 * (mrSges(6,2) * t202 - mrSges(6,3) * t200 + Ifges(6,1) * t224 + Ifges(6,4) * t225 + Ifges(6,5) * qJDD(5) - qJD(5) * t219 - t218 * t262) + qJ(3) * (-mrSges(4,1) * t265 + t257) - pkin(2) * t188 + mrSges(5,2) * t206 + mrSges(4,3) * t213 + mrSges(3,1) * t216 + (Ifges(5,3) + Ifges(3,3) + Ifges(4,2)) * t243;
t1 = [Ifges(2,3) * qJDD(1) + pkin(1) * (t250 * (m(3) * t217 - t243 * mrSges(3,2) + (-mrSges(3,1) - mrSges(4,1)) * t265 + t257) + t253 * (m(3) * t216 + t243 * mrSges(3,1) - mrSges(3,2) * t265 - t188)) + t255 - mrSges(2,2) * t258 + mrSges(2,1) * t260; t255; t188; m(5) * t246 + t252 * t198 + t249 * t199; mrSges(6,1) * t200 - mrSges(6,2) * t201 + Ifges(6,5) * t224 + Ifges(6,6) * t225 + Ifges(6,3) * qJDD(5) + (-t219 * t249 + t220 * t252) * t244;];
tauJ = t1;
