% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S4RPRR9
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
% Datum: 2019-12-31 16:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S4RPRR9_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR9_invdynJ_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR9_invdynJ_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR9_invdynJ_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR9_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR9_invdynJ_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR9_invdynJ_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR9_invdynJ_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR9_invdynJ_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:56:12
% EndTime: 2019-12-31 16:56:12
% DurationCPUTime: 0.39s
% Computational Cost: add. (1292->142), mult. (2409->183), div. (0->0), fcn. (1230->6), ass. (0->61)
t241 = sin(qJ(1));
t244 = cos(qJ(1));
t252 = -t244 * g(1) - t241 * g(2);
t263 = -qJDD(1) * qJ(2) - (2 * qJD(2) * qJD(1)) - t252;
t262 = (-pkin(1) - pkin(5));
t240 = sin(qJ(3));
t261 = t240 * g(3);
t246 = qJD(1) ^ 2;
t219 = (t262 * t246) - t263;
t243 = cos(qJ(3));
t258 = qJD(1) * qJD(3);
t255 = t243 * t258;
t231 = -t240 * qJDD(1) - t255;
t256 = t240 * t258;
t232 = t243 * qJDD(1) - t256;
t203 = (-t232 + t256) * pkin(6) + (-t231 + t255) * pkin(3) + t219;
t254 = t241 * g(1) - t244 * g(2);
t249 = -t246 * qJ(2) + qJDD(2) - t254;
t220 = t262 * qJDD(1) + t249;
t215 = -t243 * g(3) + t240 * t220;
t230 = (t240 * pkin(3) - t243 * pkin(6)) * qJD(1);
t245 = qJD(3) ^ 2;
t259 = t240 * qJD(1);
t205 = -t245 * pkin(3) + qJDD(3) * pkin(6) - t230 * t259 + t215;
t239 = sin(qJ(4));
t242 = cos(qJ(4));
t201 = t242 * t203 - t239 * t205;
t260 = qJD(1) * t243;
t227 = t242 * qJD(3) - t239 * t260;
t212 = t227 * qJD(4) + t239 * qJDD(3) + t242 * t232;
t228 = t239 * qJD(3) + t242 * t260;
t213 = -t227 * mrSges(5,1) + t228 * mrSges(5,2);
t235 = qJD(4) + t259;
t216 = -t235 * mrSges(5,2) + t227 * mrSges(5,3);
t226 = qJDD(4) - t231;
t199 = m(5) * t201 + t226 * mrSges(5,1) - t212 * mrSges(5,3) - t228 * t213 + t235 * t216;
t202 = t239 * t203 + t242 * t205;
t211 = -t228 * qJD(4) + t242 * qJDD(3) - t239 * t232;
t217 = t235 * mrSges(5,1) - t228 * mrSges(5,3);
t200 = m(5) * t202 - t226 * mrSges(5,2) + t211 * mrSges(5,3) + t227 * t213 - t235 * t217;
t192 = t242 * t199 + t239 * t200;
t253 = -t239 * t199 + t242 * t200;
t214 = t243 * t220 + t261;
t229 = (t240 * mrSges(4,1) + t243 * mrSges(4,2)) * qJD(1);
t233 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t259;
t234 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t260;
t204 = -qJDD(3) * pkin(3) - t245 * pkin(6) - t261 + (qJD(1) * t230 - t220) * t243;
t248 = -m(5) * t204 + t211 * mrSges(5,1) - t212 * mrSges(5,2) + t227 * t216 - t228 * t217;
t251 = t240 * (m(4) * t215 - qJDD(3) * mrSges(4,2) + t231 * mrSges(4,3) - qJD(3) * t234 - t229 * t259 + t253) + t243 * (m(4) * t214 + qJDD(3) * mrSges(4,1) - t232 * mrSges(4,3) + qJD(3) * t233 - t229 * t260 + t248);
t207 = Ifges(5,4) * t228 + Ifges(5,2) * t227 + Ifges(5,6) * t235;
t208 = Ifges(5,1) * t228 + Ifges(5,4) * t227 + Ifges(5,5) * t235;
t247 = mrSges(5,1) * t201 - mrSges(5,2) * t202 + Ifges(5,5) * t212 + Ifges(5,6) * t211 + Ifges(5,3) * t226 + t228 * t207 - t227 * t208;
t225 = (Ifges(4,5) * qJD(3)) + (t243 * Ifges(4,1) - t240 * Ifges(4,4)) * qJD(1);
t224 = (Ifges(4,6) * qJD(3)) + (t243 * Ifges(4,4) - t240 * Ifges(4,2)) * qJD(1);
t222 = -qJDD(1) * pkin(1) + t249;
t221 = t246 * pkin(1) + t263;
t206 = Ifges(5,5) * t228 + Ifges(5,6) * t227 + Ifges(5,3) * t235;
t194 = mrSges(5,2) * t204 - mrSges(5,3) * t201 + Ifges(5,1) * t212 + Ifges(5,4) * t211 + Ifges(5,5) * t226 + t227 * t206 - t235 * t207;
t193 = -mrSges(5,1) * t204 + mrSges(5,3) * t202 + Ifges(5,4) * t212 + Ifges(5,2) * t211 + Ifges(5,6) * t226 - t228 * t206 + t235 * t208;
t190 = m(3) * t222 + qJDD(1) * mrSges(3,2) - (t246 * mrSges(3,3)) + t251;
t1 = [mrSges(2,1) * t254 - mrSges(2,2) * t252 + mrSges(3,2) * t222 - mrSges(3,3) * t221 + t243 * (mrSges(4,2) * t219 - mrSges(4,3) * t214 + Ifges(4,1) * t232 + Ifges(4,4) * t231 + Ifges(4,5) * qJDD(3) - pkin(6) * t192 - qJD(3) * t224 - t239 * t193 + t242 * t194) - t240 * (-mrSges(4,1) * t219 + mrSges(4,3) * t215 + Ifges(4,4) * t232 + Ifges(4,2) * t231 + Ifges(4,6) * qJDD(3) - pkin(3) * t192 + qJD(3) * t225 - t247) - pkin(5) * t251 - pkin(1) * t190 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1) + (-m(3) * t221 + m(4) * t219 - t231 * mrSges(4,1) + t246 * mrSges(3,2) + t232 * mrSges(4,2) + t192 + qJDD(1) * mrSges(3,3) + (t233 * t240 + t234 * t243) * qJD(1)) * qJ(2); t190; Ifges(4,5) * t232 + Ifges(4,6) * t231 + Ifges(4,3) * qJDD(3) + mrSges(4,1) * t214 - mrSges(4,2) * t215 + t239 * t194 + t242 * t193 + pkin(3) * t248 + pkin(6) * t253 + (t243 * t224 + t240 * t225) * qJD(1); t247;];
tauJ = t1;
