% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5PRPPR5
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
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta4]';
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
% Datum: 2019-12-31 17:38
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5PRPPR5_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR5_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR5_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPPR5_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR5_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR5_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR5_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPPR5_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPPR5_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:38:07
% EndTime: 2019-12-31 17:38:07
% DurationCPUTime: 0.24s
% Computational Cost: add. (1113->101), mult. (1729->130), div. (0->0), fcn. (768->8), ass. (0->47)
t248 = -pkin(2) - pkin(3);
t238 = qJD(2) ^ 2;
t231 = sin(pkin(7));
t233 = cos(pkin(7));
t220 = -t233 * g(1) - t231 * g(2);
t228 = -g(3) + qJDD(1);
t235 = sin(qJ(2));
t237 = cos(qJ(2));
t208 = t237 * t220 + t235 * t228;
t243 = qJDD(2) * qJ(3) + (2 * qJD(3) * qJD(2)) + t208;
t202 = t248 * t238 + t243;
t207 = -t235 * t220 + t237 * t228;
t240 = -t238 * qJ(3) + qJDD(3) - t207;
t204 = t248 * qJDD(2) + t240;
t230 = sin(pkin(8));
t232 = cos(pkin(8));
t199 = t232 * t202 + t230 * t204;
t234 = sin(qJ(5));
t247 = qJD(2) * t234;
t236 = cos(qJ(5));
t246 = qJD(2) * t236;
t245 = qJD(2) * qJD(5);
t197 = -t238 * pkin(4) - qJDD(2) * pkin(6) + t199;
t219 = t231 * g(1) - t233 * g(2) + qJDD(4);
t194 = -t234 * t197 + t236 * t219;
t216 = (t236 * mrSges(6,1) - t234 * mrSges(6,2)) * qJD(2);
t217 = -t234 * qJDD(2) - t236 * t245;
t222 = -qJD(5) * mrSges(6,2) - mrSges(6,3) * t246;
t192 = m(6) * t194 + qJDD(5) * mrSges(6,1) - t217 * mrSges(6,3) + qJD(5) * t222 + t216 * t247;
t195 = t236 * t197 + t234 * t219;
t218 = -t236 * qJDD(2) + t234 * t245;
t221 = qJD(5) * mrSges(6,1) + mrSges(6,3) * t247;
t193 = m(6) * t195 - qJDD(5) * mrSges(6,2) + t218 * mrSges(6,3) - qJD(5) * t221 - t216 * t246;
t244 = -t234 * t192 + t236 * t193;
t189 = m(5) * t199 - t238 * mrSges(5,1) + qJDD(2) * mrSges(5,2) + t244;
t198 = -t230 * t202 + t232 * t204;
t196 = qJDD(2) * pkin(4) - t238 * pkin(6) - t198;
t239 = -m(6) * t196 + t218 * mrSges(6,1) - t217 * mrSges(6,2) + t221 * t247 - t222 * t246;
t190 = m(5) * t198 - qJDD(2) * mrSges(5,1) - t238 * mrSges(5,2) + t239;
t242 = t230 * t189 + t232 * t190;
t205 = -t238 * pkin(2) + t243;
t241 = m(4) * t205 + qJDD(2) * mrSges(4,3) + t232 * t189 - t230 * t190;
t206 = -qJDD(2) * pkin(2) + t240;
t187 = m(4) * t206 - qJDD(2) * mrSges(4,1) - t238 * mrSges(4,3) + t242;
t211 = (Ifges(6,5) * qJD(5)) + (-t234 * Ifges(6,1) - t236 * Ifges(6,4)) * qJD(2);
t210 = (Ifges(6,6) * qJD(5)) + (-t234 * Ifges(6,4) - t236 * Ifges(6,2)) * qJD(2);
t1 = [m(2) * t228 + t235 * (m(3) * t208 - qJDD(2) * mrSges(3,2) + (-mrSges(3,1) - mrSges(4,1)) * t238 + t241) + t237 * (m(3) * t207 + qJDD(2) * mrSges(3,1) - t238 * mrSges(3,2) - t187); -pkin(2) * t187 + qJ(3) * (-t238 * mrSges(4,1) + t241) + mrSges(3,1) * t207 - mrSges(3,2) * t208 - pkin(3) * t242 - mrSges(4,1) * t206 + mrSges(4,3) * t205 - mrSges(5,1) * t198 + mrSges(5,2) * t199 - t234 * (mrSges(6,2) * t196 - mrSges(6,3) * t194 + Ifges(6,1) * t217 + Ifges(6,4) * t218 + Ifges(6,5) * qJDD(5) - qJD(5) * t210) - t236 * (-mrSges(6,1) * t196 + mrSges(6,3) * t195 + Ifges(6,4) * t217 + Ifges(6,2) * t218 + Ifges(6,6) * qJDD(5) + qJD(5) * t211) - pkin(4) * t239 - pkin(6) * t244 + (Ifges(3,3) + Ifges(4,2) + Ifges(5,3)) * qJDD(2); t187; m(5) * t219 + t236 * t192 + t234 * t193; mrSges(6,1) * t194 - mrSges(6,2) * t195 + Ifges(6,5) * t217 + Ifges(6,6) * t218 + Ifges(6,3) * qJDD(5) + (-t234 * t210 + t236 * t211) * qJD(2);];
tauJ = t1;
