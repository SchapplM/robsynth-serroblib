% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S4PRRR2
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4]';
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
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 13:27
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S4PRRR2_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(2,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR2_invdynB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR2_invdynB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRR2_invdynB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR2_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S4PRRR2_invdynB_fixb_snew_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR2_invdynB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRR2_invdynB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRR2_invdynB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:27:23
% EndTime: 2019-07-18 13:27:23
% DurationCPUTime: 0.20s
% Computational Cost: add. (1645->98), mult. (2106->113), div. (0->0), fcn. (1072->6), ass. (0->44)
t205 = -m(1) - m(2);
t204 = m(4) + m(5);
t194 = sin(qJ(2));
t197 = cos(qJ(2));
t180 = t197 * g(1) + t194 * g(3);
t178 = qJDD(2) * pkin(1) + t180;
t181 = t194 * g(1) - t197 * g(3);
t198 = qJD(2) ^ 2;
t179 = -t198 * pkin(1) + t181;
t193 = sin(qJ(3));
t196 = cos(qJ(3));
t173 = t196 * t178 - t193 * t179;
t190 = qJD(2) + qJD(3);
t188 = t190 ^ 2;
t189 = qJDD(2) + qJDD(3);
t171 = t189 * pkin(2) + t173;
t174 = t193 * t178 + t196 * t179;
t172 = -t188 * pkin(2) + t174;
t192 = sin(qJ(4));
t195 = cos(qJ(4));
t169 = t195 * t171 - t192 * t172;
t184 = qJD(4) + t190;
t182 = t184 ^ 2;
t183 = qJDD(4) + t189;
t167 = m(5) * t169 + t183 * mrSges(5,1) - t182 * mrSges(5,2);
t170 = t192 * t171 + t195 * t172;
t168 = m(5) * t170 - t182 * mrSges(5,1) - t183 * mrSges(5,2);
t202 = t195 * t167 + t192 * t168;
t161 = m(4) * t173 + t189 * mrSges(4,1) - t188 * mrSges(4,2) + t202;
t162 = m(4) * t174 - t188 * mrSges(4,1) - t189 * mrSges(4,2) - t192 * t167 + t195 * t168;
t203 = t196 * t161 + t193 * t162;
t155 = m(3) * t180 + qJDD(2) * mrSges(3,1) - t198 * mrSges(3,2) + t203;
t156 = m(3) * t181 - t198 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t193 * t161 + t196 * t162;
t201 = -t194 * t155 + t197 * t156;
t200 = qJ(1) * m(2) - mrSges(1,2) + mrSges(2,3);
t199 = -t197 * t155 - t194 * t156;
t191 = g(2) + qJDD(1);
t164 = mrSges(5,2) * t191 - mrSges(5,3) * t169 + Ifges(5,5) * t183 - t182 * Ifges(5,6);
t163 = -mrSges(5,1) * t191 + mrSges(5,3) * t170 + t182 * Ifges(5,5) + Ifges(5,6) * t183;
t158 = mrSges(4,2) * t191 - mrSges(4,3) * t173 + Ifges(4,5) * t189 - t188 * Ifges(4,6) - t192 * t163 + t195 * t164;
t157 = mrSges(4,3) * t174 + t188 * Ifges(4,5) + Ifges(4,6) * t189 + t195 * t163 + t192 * t164 + (-m(5) * pkin(2) - mrSges(4,1)) * t191;
t153 = mrSges(3,2) * t191 - mrSges(3,3) * t180 + Ifges(3,5) * qJDD(2) - t198 * Ifges(3,6) - t193 * t157 + t196 * t158;
t152 = mrSges(3,3) * t181 + t198 * Ifges(3,5) + Ifges(3,6) * qJDD(2) + t196 * t157 + t193 * t158 + (-pkin(1) * t204 - mrSges(3,1)) * t191;
t1 = [t205 * g(1) + t199; -m(1) * g(2) + (-m(2) - m(3) - t204) * t191; t205 * g(3) + t201; mrSges(2,1) * t191 + mrSges(1,3) * g(2) + t200 * g(3) - qJ(1) * t201 - t197 * t152 - t194 * t153; -pkin(1) * t203 - mrSges(5,1) * t169 + mrSges(5,2) * t170 - Ifges(4,3) * t189 - mrSges(4,1) * t173 + mrSges(4,2) * t174 - pkin(2) * t202 - Ifges(3,3) * qJDD(2) - Ifges(5,3) * t183 - mrSges(3,1) * t180 + mrSges(3,2) * t181 + (mrSges(1,1) - mrSges(2,2)) * g(3) + (-mrSges(2,1) - mrSges(1,3)) * g(1); -mrSges(1,1) * g(2) + mrSges(2,2) * t191 - t200 * g(1) + qJ(1) * t199 - t194 * t152 + t197 * t153;];
tauB  = t1;
