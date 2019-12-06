% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5PRRRR2
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,d5]';
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
% Datum: 2019-12-05 17:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5PRRRR2_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR2_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR2_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR2_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR2_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5PRRRR2_invdynJ_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR2_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR2_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR2_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:04:48
% EndTime: 2019-12-05 17:04:48
% DurationCPUTime: 0.20s
% Computational Cost: add. (1636->94), mult. (1855->126), div. (0->0), fcn. (984->8), ass. (0->49)
t229 = qJD(2) + qJD(3);
t225 = qJD(4) + t229;
t231 = sin(qJ(5));
t247 = t225 * t231;
t235 = cos(qJ(5));
t246 = t225 * t235;
t234 = sin(qJ(2));
t238 = cos(qJ(2));
t243 = t234 * g(1) - t238 * g(2);
t219 = qJDD(2) * pkin(2) + t243;
t241 = -t238 * g(1) - t234 * g(2);
t220 = -qJD(2) ^ 2 * pkin(2) + t241;
t233 = sin(qJ(3));
t237 = cos(qJ(3));
t206 = t237 * t219 - t233 * t220;
t228 = qJDD(2) + qJDD(3);
t203 = t228 * pkin(3) + t206;
t207 = t233 * t219 + t237 * t220;
t227 = t229 ^ 2;
t204 = -t227 * pkin(3) + t207;
t232 = sin(qJ(4));
t236 = cos(qJ(4));
t200 = t232 * t203 + t236 * t204;
t223 = t225 ^ 2;
t224 = qJDD(4) + t228;
t196 = t224 * pkin(6) + t200;
t230 = -g(3) + qJDD(1);
t194 = -t231 * t196 + t235 * t230;
t211 = (-mrSges(6,1) * t235 + mrSges(6,2) * t231) * t225;
t244 = qJD(5) * t225;
t212 = t231 * t224 + t235 * t244;
t218 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t246;
t192 = m(6) * t194 + qJDD(5) * mrSges(6,1) - t212 * mrSges(6,3) + qJD(5) * t218 - t211 * t247;
t195 = t235 * t196 + t231 * t230;
t213 = t235 * t224 - t231 * t244;
t217 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t247;
t193 = m(6) * t195 - qJDD(5) * mrSges(6,2) + t213 * mrSges(6,3) - qJD(5) * t217 + t211 * t246;
t242 = -t231 * t192 + t235 * t193;
t185 = m(5) * t200 - t223 * mrSges(5,1) - t224 * mrSges(5,2) + t242;
t199 = t236 * t203 - t232 * t204;
t197 = -t223 * pkin(6) - t199;
t190 = m(5) * t199 - m(6) * t197 + t224 * mrSges(5,1) + t213 * mrSges(6,1) - t223 * mrSges(5,2) - t212 * mrSges(6,2) + (-t217 * t231 + t218 * t235) * t225;
t245 = t232 * t185 + t236 * t190;
t208 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t231 + Ifges(6,6) * t235) * t225;
t209 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t231 + Ifges(6,2) * t235) * t225;
t210 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t231 + Ifges(6,4) * t235) * t225;
t240 = -mrSges(5,2) * t200 + pkin(6) * t242 + t231 * (mrSges(6,2) * t197 - mrSges(6,3) * t194 + Ifges(6,1) * t212 + Ifges(6,4) * t213 + Ifges(6,5) * qJDD(5) - qJD(5) * t209 + t208 * t246) + t235 * (-mrSges(6,1) * t197 + mrSges(6,3) * t195 + Ifges(6,4) * t212 + Ifges(6,2) * t213 + Ifges(6,6) * qJDD(5) + qJD(5) * t210 - t208 * t247) + mrSges(5,1) * t199 + Ifges(5,3) * t224;
t239 = mrSges(4,1) * t206 - mrSges(4,2) * t207 + Ifges(4,3) * t228 + pkin(3) * t245 + t240;
t1 = [t235 * t192 + t231 * t193 + (m(2) + m(3) + m(4) + m(5)) * t230; t239 + Ifges(3,3) * qJDD(2) + pkin(2) * (t233 * (m(4) * t207 - t227 * mrSges(4,1) - t228 * mrSges(4,2) + t236 * t185 - t232 * t190) + t237 * (m(4) * t206 + t228 * mrSges(4,1) - t227 * mrSges(4,2) + t245)) + mrSges(3,1) * t243 - mrSges(3,2) * t241; t239; t240; mrSges(6,1) * t194 - mrSges(6,2) * t195 + Ifges(6,5) * t212 + Ifges(6,6) * t213 + Ifges(6,3) * qJDD(5) + (t209 * t231 - t210 * t235) * t225;];
tauJ = t1;
