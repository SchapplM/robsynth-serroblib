% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S4PRRR3
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
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
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
% Datum: 2019-12-31 16:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S4PRRR3_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR3_invdynJ_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR3_invdynJ_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRR3_invdynJ_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR3_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR3_invdynJ_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR3_invdynJ_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRR3_invdynJ_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRR3_invdynJ_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:31:37
% EndTime: 2019-12-31 16:31:37
% DurationCPUTime: 0.15s
% Computational Cost: add. (658->80), mult. (907->111), div. (0->0), fcn. (550->8), ass. (0->41)
t186 = qJD(2) + qJD(3);
t190 = sin(qJ(4));
t203 = t186 * t190;
t193 = cos(qJ(4));
t202 = t186 * t193;
t188 = sin(pkin(7));
t189 = cos(pkin(7));
t181 = t188 * g(1) - t189 * g(2);
t182 = -t189 * g(1) - t188 * g(2);
t192 = sin(qJ(2));
t195 = cos(qJ(2));
t198 = t195 * t181 - t192 * t182;
t166 = qJDD(2) * pkin(2) + t198;
t201 = t192 * t181 + t195 * t182;
t167 = -qJD(2) ^ 2 * pkin(2) + t201;
t191 = sin(qJ(3));
t194 = cos(qJ(3));
t163 = t191 * t166 + t194 * t167;
t200 = qJD(4) * t186;
t184 = t186 ^ 2;
t185 = qJDD(2) + qJDD(3);
t160 = -t184 * pkin(3) + t185 * pkin(6) + t163;
t187 = -g(3) + qJDD(1);
t157 = -t190 * t160 + t193 * t187;
t173 = (-mrSges(5,1) * t193 + mrSges(5,2) * t190) * t186;
t174 = t190 * t185 + t193 * t200;
t180 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t202;
t155 = m(5) * t157 + qJDD(4) * mrSges(5,1) - t174 * mrSges(5,3) + qJD(4) * t180 - t173 * t203;
t158 = t193 * t160 + t190 * t187;
t175 = t193 * t185 - t190 * t200;
t179 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t203;
t156 = m(5) * t158 - qJDD(4) * mrSges(5,2) + t175 * mrSges(5,3) - qJD(4) * t179 + t173 * t202;
t199 = -t190 * t155 + t193 * t156;
t162 = t194 * t166 - t191 * t167;
t159 = -t185 * pkin(3) - t184 * pkin(6) - t162;
t168 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t190 + Ifges(5,6) * t193) * t186;
t169 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t190 + Ifges(5,2) * t193) * t186;
t170 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t190 + Ifges(5,4) * t193) * t186;
t196 = -m(5) * t159 + t175 * mrSges(5,1) - t174 * mrSges(5,2) - t179 * t203 + t180 * t202;
t197 = -mrSges(4,2) * t163 + pkin(6) * t199 + t190 * (mrSges(5,2) * t159 - mrSges(5,3) * t157 + Ifges(5,1) * t174 + Ifges(5,4) * t175 + Ifges(5,5) * qJDD(4) - qJD(4) * t169 + t168 * t202) + t193 * (-mrSges(5,1) * t159 + mrSges(5,3) * t158 + Ifges(5,4) * t174 + Ifges(5,2) * t175 + Ifges(5,6) * qJDD(4) + qJD(4) * t170 - t168 * t203) + pkin(3) * t196 + mrSges(4,1) * t162 + Ifges(4,3) * t185;
t1 = [t193 * t155 + t190 * t156 + (m(2) + m(3) + m(4)) * t187; Ifges(3,3) * qJDD(2) + mrSges(3,1) * t198 - mrSges(3,2) * t201 + pkin(2) * (t191 * (m(4) * t163 - t184 * mrSges(4,1) - t185 * mrSges(4,2) + t199) + t194 * (m(4) * t162 + t185 * mrSges(4,1) - t184 * mrSges(4,2) + t196)) + t197; t197; mrSges(5,1) * t157 - mrSges(5,2) * t158 + Ifges(5,5) * t174 + Ifges(5,6) * t175 + Ifges(5,3) * qJDD(4) + (t169 * t190 - t170 * t193) * t186;];
tauJ = t1;
