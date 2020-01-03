% Calculate vector of cutting torques with Newton-Euler for
% S4RRRP3
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
%   pkin=[a2,a3,a4,d1,d2,d3]';
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
% m [3x5]
%   vector of cutting torques (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:14
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S4RRRP3_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP3_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP3_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP3_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP3_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP3_invdynm_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP3_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRP3_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRP3_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:14:02
% EndTime: 2019-12-31 17:14:04
% DurationCPUTime: 0.95s
% Computational Cost: add. (9741->190), mult. (12249->239), div. (0->0), fcn. (5457->6), ass. (0->71)
t161 = sin(qJ(1));
t164 = cos(qJ(1));
t146 = t161 * g(1) - t164 * g(2);
t139 = qJDD(1) * pkin(1) + t146;
t147 = -t164 * g(1) - t161 * g(2);
t166 = qJD(1) ^ 2;
t140 = -t166 * pkin(1) + t147;
t160 = sin(qJ(2));
t163 = cos(qJ(2));
t110 = t160 * t139 + t163 * t140;
t153 = qJD(1) + qJD(2);
t151 = t153 ^ 2;
t152 = qJDD(1) + qJDD(2);
t107 = -t151 * pkin(2) + t152 * pkin(6) + t110;
t159 = sin(qJ(3));
t162 = cos(qJ(3));
t184 = t162 * g(3);
t103 = -t159 * t107 - t184;
t104 = -t159 * g(3) + t162 * t107;
t114 = Ifges(5,6) * qJD(3) + (Ifges(5,5) * t159 - Ifges(5,3) * t162) * t153;
t117 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t159 + Ifges(4,2) * t162) * t153;
t130 = (-mrSges(5,1) * t162 - mrSges(5,3) * t159) * t153;
t178 = qJD(3) * t153;
t132 = t159 * t152 + t162 * t178;
t133 = t162 * t152 - t159 * t178;
t129 = (-pkin(3) * t162 - qJ(4) * t159) * t153;
t165 = qJD(3) ^ 2;
t181 = t153 * t162;
t100 = -t165 * pkin(3) + qJDD(3) * qJ(4) + 0.2e1 * qJD(4) * qJD(3) + t129 * t181 + t104;
t102 = -qJDD(3) * pkin(3) + t184 - t165 * qJ(4) + qJDD(4) + (t129 * t153 + t107) * t159;
t172 = -mrSges(5,1) * t102 + mrSges(5,3) * t100 + Ifges(5,4) * t132 + Ifges(5,2) * qJDD(3) - Ifges(5,6) * t133;
t144 = mrSges(5,2) * t181 + qJD(3) * mrSges(5,3);
t174 = -m(5) * t102 + qJDD(3) * mrSges(5,1) + qJD(3) * t144;
t182 = t153 * t159;
t142 = -qJD(3) * mrSges(5,1) + mrSges(5,2) * t182;
t175 = m(5) * t100 + qJDD(3) * mrSges(5,3) + qJD(3) * t142 + t130 * t181;
t118 = Ifges(5,4) * qJD(3) + (Ifges(5,1) * t159 - Ifges(5,5) * t162) * t153;
t179 = t118 + Ifges(4,5) * qJD(3) + (Ifges(4,1) * t159 + Ifges(4,4) * t162) * t153;
t186 = -((t114 - t117) * t159 + t162 * t179) * t153 + mrSges(4,1) * t103 - mrSges(4,2) * t104 + Ifges(4,5) * t132 + Ifges(4,6) * t133 + Ifges(4,3) * qJDD(3) + pkin(3) * (-t132 * mrSges(5,2) - t130 * t182 + t174) + qJ(4) * (t133 * mrSges(5,2) + t175) + t172;
t183 = mrSges(4,3) + mrSges(5,2);
t131 = (-mrSges(4,1) * t162 + mrSges(4,2) * t159) * t153;
t141 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t182;
t92 = m(4) * t104 - qJDD(3) * mrSges(4,2) - qJD(3) * t141 + t131 * t181 + t133 * t183 + t175;
t143 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t181;
t93 = m(4) * t103 + qJDD(3) * mrSges(4,1) + qJD(3) * t143 + (-t130 - t131) * t182 - t183 * t132 + t174;
t177 = -t159 * t93 + t162 * t92;
t83 = m(3) * t110 - t151 * mrSges(3,1) - t152 * mrSges(3,2) + t177;
t109 = t163 * t139 - t160 * t140;
t106 = -t152 * pkin(2) - t151 * pkin(6) - t109;
t97 = -t133 * pkin(3) - t132 * qJ(4) + (-0.2e1 * qJD(4) * t159 + (pkin(3) * t159 - qJ(4) * t162) * qJD(3)) * t153 + t106;
t94 = m(5) * t97 - t133 * mrSges(5,1) - t132 * mrSges(5,3) - t142 * t182 - t144 * t181;
t168 = -m(4) * t106 + t133 * mrSges(4,1) - t132 * mrSges(4,2) - t141 * t182 + t143 * t181 - t94;
t87 = m(3) * t109 + t152 * mrSges(3,1) - t151 * mrSges(3,2) + t168;
t76 = t160 * t83 + t163 * t87;
t85 = t159 * t92 + t162 * t93;
t176 = -t160 * t87 + t163 * t83;
t173 = -mrSges(5,1) * t97 + mrSges(5,2) * t100;
t115 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t159 + Ifges(4,6) * t162) * t153;
t116 = Ifges(5,2) * qJD(3) + (Ifges(5,4) * t159 - Ifges(5,6) * t162) * t153;
t79 = -mrSges(4,1) * t106 + mrSges(4,3) * t104 - pkin(3) * t94 + (-t115 - t116) * t182 + (Ifges(4,2) + Ifges(5,3)) * t133 + (Ifges(4,4) - Ifges(5,5)) * t132 + (Ifges(4,6) - Ifges(5,6)) * qJDD(3) + t179 * qJD(3) + t173;
t170 = mrSges(5,2) * t102 - mrSges(5,3) * t97 + Ifges(5,1) * t132 + Ifges(5,4) * qJDD(3) - Ifges(5,5) * t133 + qJD(3) * t114 + t116 * t181;
t80 = mrSges(4,2) * t106 - mrSges(4,3) * t103 + Ifges(4,1) * t132 + Ifges(4,4) * t133 + Ifges(4,5) * qJDD(3) - qJ(4) * t94 - qJD(3) * t117 + t115 * t181 + t170;
t171 = mrSges(3,1) * t109 - mrSges(3,2) * t110 + Ifges(3,3) * t152 + pkin(2) * t168 + pkin(6) * t177 + t159 * t80 + t162 * t79;
t169 = mrSges(2,1) * t146 - mrSges(2,2) * t147 + Ifges(2,3) * qJDD(1) + pkin(1) * t76 + t171;
t74 = m(2) * t147 - t166 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t176;
t73 = m(2) * t146 + qJDD(1) * mrSges(2,1) - t166 * mrSges(2,2) + t76;
t72 = mrSges(3,1) * g(3) + mrSges(3,3) * t110 + t151 * Ifges(3,5) + Ifges(3,6) * t152 - pkin(2) * t85 - t186;
t71 = -mrSges(3,2) * g(3) - mrSges(3,3) * t109 + Ifges(3,5) * t152 - t151 * Ifges(3,6) - pkin(6) * t85 - t159 * t79 + t162 * t80;
t70 = -mrSges(2,2) * g(3) - mrSges(2,3) * t146 + Ifges(2,5) * qJDD(1) - t166 * Ifges(2,6) - pkin(5) * t76 - t160 * t72 + t163 * t71;
t69 = Ifges(2,6) * qJDD(1) + t166 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t147 + t160 * t71 + t163 * t72 - pkin(1) * (-m(3) * g(3) + t85) + pkin(5) * t176;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t164 * t70 - t161 * t69 - pkin(4) * (t161 * t74 + t164 * t73), t70, t71, t80, t170; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t161 * t70 + t164 * t69 + pkin(4) * (-t161 * t73 + t164 * t74), t69, t72, t79, (-t159 * t114 - t162 * t118) * t153 + t172; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t169, t169, t171, t186, Ifges(5,5) * t132 + Ifges(5,6) * qJDD(3) - Ifges(5,3) * t133 - qJD(3) * t118 + t116 * t182 - t173;];
m_new = t1;
