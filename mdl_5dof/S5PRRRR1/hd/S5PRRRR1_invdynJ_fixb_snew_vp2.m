% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5PRRRR1
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
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
% Datum: 2019-12-05 17:03
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5PRRRR1_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(2,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR1_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR1_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR1_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR1_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5PRRRR1_invdynJ_fixb_snew_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR1_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR1_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR1_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:02:59
% EndTime: 2019-12-05 17:03:00
% DurationCPUTime: 0.54s
% Computational Cost: add. (2353->181), mult. (4502->239), div. (0->0), fcn. (3209->8), ass. (0->69)
t248 = sin(qJ(4));
t249 = sin(qJ(3));
t252 = cos(qJ(4));
t253 = cos(qJ(3));
t232 = (t249 * t248 - t253 * t252) * qJD(2);
t246 = -g(3) + qJDD(1);
t250 = sin(qJ(2));
t254 = cos(qJ(2));
t239 = -t254 * g(1) + t250 * t246;
t227 = t249 * g(2) + t253 * t239;
t255 = qJD(2) ^ 2;
t221 = (-t253 ^ 2 * t255 - qJD(3) ^ 2) * pkin(2) + t227;
t226 = t253 * g(2) - t249 * t239;
t259 = (t249 * t253 * t255 + qJDD(3)) * pkin(2) + t226;
t205 = t252 * t221 + t248 * t259;
t262 = qJD(2) * qJD(3);
t261 = t249 * t262;
t237 = t253 * qJDD(2) - t261;
t238 = t250 * g(1) + t254 * t246;
t220 = (-t237 + t261) * pkin(2) - t238;
t247 = sin(qJ(5));
t251 = cos(qJ(5));
t196 = -t247 * t205 + t251 * t220;
t236 = t249 * qJDD(2) + t253 * t262;
t212 = -t232 * qJD(4) + t252 * t236 + t248 * t237;
t233 = (t253 * t248 + t249 * t252) * qJD(2);
t245 = qJD(3) + qJD(4);
t222 = -t247 * t233 + t251 * t245;
t244 = qJDD(3) + qJDD(4);
t199 = t222 * qJD(5) + t251 * t212 + t247 * t244;
t223 = t251 * t233 + t247 * t245;
t206 = -t222 * mrSges(6,1) + t223 * mrSges(6,2);
t211 = -t233 * qJD(4) - t248 * t236 + t252 * t237;
t210 = qJDD(5) - t211;
t228 = qJD(5) + t232;
t213 = -t228 * mrSges(6,2) + t222 * mrSges(6,3);
t194 = m(6) * t196 + t210 * mrSges(6,1) - t199 * mrSges(6,3) - t223 * t206 + t228 * t213;
t197 = t251 * t205 + t247 * t220;
t198 = -t223 * qJD(5) - t247 * t212 + t251 * t244;
t214 = t228 * mrSges(6,1) - t223 * mrSges(6,3);
t195 = m(6) * t197 - t210 * mrSges(6,2) + t198 * mrSges(6,3) + t222 * t206 - t228 * t214;
t218 = t232 * mrSges(5,1) + t233 * mrSges(5,2);
t225 = t245 * mrSges(5,1) - t233 * mrSges(5,3);
t186 = m(5) * t205 - t244 * mrSges(5,2) + t211 * mrSges(5,3) - t247 * t194 + t251 * t195 - t232 * t218 - t245 * t225;
t204 = t248 * t221 - t252 * t259;
t224 = -t245 * mrSges(5,2) - t232 * mrSges(5,3);
t193 = t244 * mrSges(5,1) + t198 * mrSges(6,1) - t199 * mrSges(6,2) - t212 * mrSges(5,3) + t222 * t213 - t223 * t214 - t233 * t218 + t245 * t224 + (-m(5) - m(6)) * t204;
t265 = t248 * t186 + t252 * t193;
t264 = qJD(2) * t249;
t263 = qJD(2) * t253;
t200 = Ifges(6,5) * t223 + Ifges(6,6) * t222 + Ifges(6,3) * t228;
t202 = Ifges(6,1) * t223 + Ifges(6,4) * t222 + Ifges(6,5) * t228;
t190 = -mrSges(6,1) * t204 + mrSges(6,3) * t197 + Ifges(6,4) * t199 + Ifges(6,2) * t198 + Ifges(6,6) * t210 - t223 * t200 + t228 * t202;
t201 = Ifges(6,4) * t223 + Ifges(6,2) * t222 + Ifges(6,6) * t228;
t191 = mrSges(6,2) * t204 - mrSges(6,3) * t196 + Ifges(6,1) * t199 + Ifges(6,4) * t198 + Ifges(6,5) * t210 + t222 * t200 - t228 * t201;
t216 = Ifges(5,4) * t233 - Ifges(5,2) * t232 + Ifges(5,6) * t245;
t217 = Ifges(5,1) * t233 - Ifges(5,4) * t232 + Ifges(5,5) * t245;
t258 = -mrSges(5,1) * t204 - mrSges(5,2) * t205 + Ifges(5,5) * t212 + Ifges(5,6) * t211 + Ifges(5,3) * t244 + t251 * t190 + t247 * t191 + t233 * t216 + t232 * t217;
t257 = mrSges(6,1) * t196 - mrSges(6,2) * t197 + Ifges(6,5) * t199 + Ifges(6,6) * t198 + Ifges(6,3) * t210 + t223 * t201 - t222 * t202;
t256 = m(5) * t220 - t211 * mrSges(5,1) + t212 * mrSges(5,2) + t251 * t194 + t247 * t195 + t232 * t224 + t233 * t225;
t241 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t263;
t240 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t264;
t235 = (-t253 * mrSges(4,1) + t249 * mrSges(4,2)) * qJD(2);
t231 = Ifges(4,5) * qJD(3) + (t249 * Ifges(4,1) + t253 * Ifges(4,4)) * qJD(2);
t230 = Ifges(4,6) * qJD(3) + (t249 * Ifges(4,4) + t253 * Ifges(4,2)) * qJD(2);
t215 = Ifges(5,5) * t233 - Ifges(5,6) * t232 + Ifges(5,3) * t245;
t187 = -mrSges(5,1) * t220 + mrSges(5,3) * t205 + Ifges(5,4) * t212 + Ifges(5,2) * t211 + Ifges(5,6) * t244 - t233 * t215 + t245 * t217 - t257;
t184 = mrSges(5,2) * t220 + mrSges(5,3) * t204 + Ifges(5,1) * t212 + Ifges(5,4) * t211 + Ifges(5,5) * t244 - t247 * t190 + t251 * t191 - t232 * t215 - t245 * t216;
t1 = [m(2) * t246 + t250 * (m(3) * t239 - qJDD(2) * mrSges(3,2) - t255 * mrSges(3,1) + t253 * (m(4) * t227 - qJDD(3) * mrSges(4,2) + t237 * mrSges(4,3) - qJD(3) * t240 + t252 * t186 - t248 * t193 + t235 * t263) - t249 * (m(4) * t226 + qJDD(3) * mrSges(4,1) - t236 * mrSges(4,3) + qJD(3) * t241 - t235 * t264 + t265)) + t254 * (qJDD(2) * mrSges(3,1) + t237 * mrSges(4,1) - t255 * mrSges(3,2) - t236 * mrSges(4,2) + (m(3) + m(4)) * t238 + (-t240 * t249 + t241 * t253) * qJD(2) - t256); Ifges(3,3) * qJDD(2) + mrSges(3,1) * t238 - mrSges(3,2) * t239 + t249 * (-mrSges(4,2) * t238 - mrSges(4,3) * t226 + Ifges(4,1) * t236 + Ifges(4,4) * t237 + Ifges(4,5) * qJDD(3) - qJD(3) * t230 + t252 * t184 - t248 * t187) + t253 * (mrSges(4,1) * t238 + mrSges(4,3) * t227 + Ifges(4,4) * t236 + Ifges(4,2) * t237 + Ifges(4,6) * qJDD(3) - pkin(2) * t256 + qJD(3) * t231 + t248 * t184 + t252 * t187); (t249 * t230 - t253 * t231) * qJD(2) + pkin(2) * t265 + t258 + Ifges(4,3) * qJDD(3) + Ifges(4,5) * t236 + Ifges(4,6) * t237 + mrSges(4,1) * t226 - mrSges(4,2) * t227; t258; t257;];
tauJ = t1;
