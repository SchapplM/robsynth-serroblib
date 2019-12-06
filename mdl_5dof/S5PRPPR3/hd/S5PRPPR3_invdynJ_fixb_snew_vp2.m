% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5PRPPR3
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
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:27
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5PRPPR3_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR3_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR3_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPPR3_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR3_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR3_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR3_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPPR3_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPPR3_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:26:25
% EndTime: 2019-12-05 15:26:25
% DurationCPUTime: 0.25s
% Computational Cost: add. (832->99), mult. (1289->127), div. (0->0), fcn. (690->8), ass. (0->47)
t246 = -pkin(3) - pkin(6);
t229 = sin(pkin(7));
t231 = cos(pkin(7));
t219 = -t231 * g(1) - t229 * g(2);
t225 = -g(3) + qJDD(1);
t233 = sin(qJ(2));
t235 = cos(qJ(2));
t206 = -t233 * t219 + t235 * t225;
t204 = qJDD(2) * pkin(2) + t206;
t207 = t235 * t219 + t233 * t225;
t236 = qJD(2) ^ 2;
t205 = -t236 * pkin(2) + t207;
t228 = sin(pkin(8));
t230 = cos(pkin(8));
t199 = t230 * t204 - t228 * t205;
t240 = -t236 * qJ(4) + qJDD(4) - t199;
t198 = -qJDD(2) * pkin(3) + t240;
t196 = t246 * qJDD(2) + t240;
t218 = -t229 * g(1) + t231 * g(2) + qJDD(3);
t232 = sin(qJ(5));
t234 = cos(qJ(5));
t192 = t234 * t196 - t232 * t218;
t215 = (t232 * mrSges(6,1) + t234 * mrSges(6,2)) * qJD(2);
t242 = qJD(2) * qJD(5);
t217 = t234 * qJDD(2) - t232 * t242;
t244 = qJD(2) * t232;
t220 = -qJD(5) * mrSges(6,2) - mrSges(6,3) * t244;
t243 = qJD(2) * t234;
t190 = m(6) * t192 + qJDD(5) * mrSges(6,1) - t217 * mrSges(6,3) + qJD(5) * t220 - t215 * t243;
t193 = t232 * t196 + t234 * t218;
t216 = -t232 * qJDD(2) - t234 * t242;
t221 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t243;
t191 = m(6) * t193 - qJDD(5) * mrSges(6,2) + t216 * mrSges(6,3) - qJD(5) * t221 - t215 * t244;
t241 = t234 * t190 + t232 * t191;
t238 = -m(5) * t198 + t236 * mrSges(5,3) - t241;
t186 = m(4) * t199 - t236 * mrSges(4,2) + (mrSges(4,1) - mrSges(5,2)) * qJDD(2) + t238;
t200 = t228 * t204 + t230 * t205;
t239 = qJDD(2) * qJ(4) + 0.2e1 * qJD(4) * qJD(2) + t200;
t195 = t246 * t236 + t239;
t197 = t236 * pkin(3) - t239;
t237 = -m(5) * t197 + m(6) * t195 - t216 * mrSges(6,1) + t236 * mrSges(5,2) + t217 * mrSges(6,2) + qJDD(2) * mrSges(5,3) + t220 * t244 + t221 * t243;
t189 = m(4) * t200 - t236 * mrSges(4,1) - qJDD(2) * mrSges(4,2) + t237;
t245 = t230 * t186 + t228 * t189;
t210 = Ifges(6,5) * qJD(5) + (t234 * Ifges(6,1) - t232 * Ifges(6,4)) * qJD(2);
t209 = Ifges(6,6) * qJD(5) + (t234 * Ifges(6,4) - t232 * Ifges(6,2)) * qJD(2);
t187 = qJDD(2) * mrSges(5,2) - t238;
t1 = [m(2) * t225 + t233 * (m(3) * t207 - t236 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t228 * t186 + t230 * t189) + t235 * (m(3) * t206 + qJDD(2) * mrSges(3,1) - t236 * mrSges(3,2) + t245); pkin(2) * t245 + mrSges(3,1) * t206 - mrSges(3,2) * t207 - pkin(3) * t187 + qJ(4) * t237 + t234 * (mrSges(6,2) * t195 - mrSges(6,3) * t192 + Ifges(6,1) * t217 + Ifges(6,4) * t216 + Ifges(6,5) * qJDD(5) - qJD(5) * t209) - t232 * (-mrSges(6,1) * t195 + mrSges(6,3) * t193 + Ifges(6,4) * t217 + Ifges(6,2) * t216 + Ifges(6,6) * qJDD(5) + qJD(5) * t210) - pkin(6) * t241 + mrSges(4,1) * t199 - mrSges(4,2) * t200 + mrSges(5,2) * t198 - mrSges(5,3) * t197 + (Ifges(3,3) + Ifges(4,3) + Ifges(5,1)) * qJDD(2); -t232 * t190 + t234 * t191 + (m(4) + m(5)) * t218; t187; mrSges(6,1) * t192 - mrSges(6,2) * t193 + Ifges(6,5) * t217 + Ifges(6,6) * t216 + Ifges(6,3) * qJDD(5) + (t234 * t209 + t232 * t210) * qJD(2);];
tauJ = t1;
