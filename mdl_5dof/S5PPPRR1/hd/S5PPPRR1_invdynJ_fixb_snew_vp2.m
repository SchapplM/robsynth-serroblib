% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5PPPRR1
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d4,d5,theta1,theta2,theta3]';
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
% Datum: 2019-12-05 14:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5PPPRR1_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR1_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPPRR1_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPPRR1_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPPRR1_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR1_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPPRR1_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPPRR1_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPPRR1_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 14:57:59
% EndTime: 2019-12-05 14:57:59
% DurationCPUTime: 0.22s
% Computational Cost: add. (874->81), mult. (1342->116), div. (0->0), fcn. (978->10), ass. (0->45)
t192 = sin(pkin(7));
t195 = cos(pkin(7));
t186 = -t195 * g(1) - t192 * g(2);
t189 = -g(3) + qJDD(1);
t191 = sin(pkin(8));
t194 = cos(pkin(8));
t176 = t194 * t186 + t191 * t189;
t185 = -t192 * g(1) + t195 * g(2) + qJDD(2);
t190 = sin(pkin(9));
t193 = cos(pkin(9));
t173 = -t190 * t176 + t193 * t185;
t174 = t193 * t176 + t190 * t185;
t197 = sin(qJ(4));
t199 = cos(qJ(4));
t170 = t197 * t173 + t199 * t174;
t200 = qJD(4) ^ 2;
t168 = -t200 * pkin(4) + qJDD(4) * pkin(6) + t170;
t202 = -t191 * t186 + t194 * t189;
t175 = qJDD(3) - t202;
t196 = sin(qJ(5));
t198 = cos(qJ(5));
t165 = -t196 * t168 + t198 * t175;
t182 = (-t198 * mrSges(6,1) + t196 * mrSges(6,2)) * qJD(4);
t205 = qJD(4) * qJD(5);
t183 = t196 * qJDD(4) + t198 * t205;
t206 = qJD(4) * t198;
t188 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t206;
t207 = qJD(4) * t196;
t163 = m(6) * t165 + qJDD(5) * mrSges(6,1) - t183 * mrSges(6,3) + qJD(5) * t188 - t182 * t207;
t166 = t198 * t168 + t196 * t175;
t184 = t198 * qJDD(4) - t196 * t205;
t187 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t207;
t164 = m(6) * t166 - qJDD(5) * mrSges(6,2) + t184 * mrSges(6,3) - qJD(5) * t187 + t182 * t206;
t209 = t198 * t163 + t196 * t164 + (m(4) + m(5)) * t175;
t204 = -t196 * t163 + t198 * t164;
t169 = t199 * t173 - t197 * t174;
t167 = -qJDD(4) * pkin(4) - t200 * pkin(6) - t169;
t201 = -m(6) * t167 + t184 * mrSges(6,1) - t183 * mrSges(6,2) - t187 * t207 + t188 * t206;
t179 = Ifges(6,5) * qJD(5) + (t196 * Ifges(6,1) + t198 * Ifges(6,4)) * qJD(4);
t178 = Ifges(6,6) * qJD(5) + (t196 * Ifges(6,4) + t198 * Ifges(6,2)) * qJD(4);
t161 = m(5) * t169 + qJDD(4) * mrSges(5,1) - t200 * mrSges(5,2) + t201;
t160 = m(5) * t170 - t200 * mrSges(5,1) - qJDD(4) * mrSges(5,2) + t204;
t159 = m(4) * t174 + t199 * t160 - t197 * t161;
t158 = m(4) * t173 + t197 * t160 + t199 * t161;
t1 = [m(2) * t189 + t191 * (m(3) * t176 - t190 * t158 + t193 * t159) + t194 * (m(3) * t202 - t209); m(3) * t185 + t193 * t158 + t190 * t159; t209; Ifges(5,3) * qJDD(4) + mrSges(5,1) * t169 - mrSges(5,2) * t170 + t196 * (mrSges(6,2) * t167 - mrSges(6,3) * t165 + Ifges(6,1) * t183 + Ifges(6,4) * t184 + Ifges(6,5) * qJDD(5) - qJD(5) * t178) + t198 * (-mrSges(6,1) * t167 + mrSges(6,3) * t166 + Ifges(6,4) * t183 + Ifges(6,2) * t184 + Ifges(6,6) * qJDD(5) + qJD(5) * t179) + pkin(4) * t201 + pkin(6) * t204; mrSges(6,1) * t165 - mrSges(6,2) * t166 + Ifges(6,5) * t183 + Ifges(6,6) * t184 + Ifges(6,3) * qJDD(5) + (t196 * t178 - t198 * t179) * qJD(4);];
tauJ = t1;
