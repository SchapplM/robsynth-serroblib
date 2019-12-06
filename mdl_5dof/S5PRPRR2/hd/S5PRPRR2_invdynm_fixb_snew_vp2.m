% Calculate vector of cutting torques with Newton-Euler for
% S5PRPRR2
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
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
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
% m [3x6]
%   vector of cutting torques (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:45
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5PRPRR2_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR2_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR2_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR2_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR2_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR2_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR2_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR2_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR2_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:44:52
% EndTime: 2019-12-05 15:44:56
% DurationCPUTime: 2.35s
% Computational Cost: add. (36458->170), mult. (48237->219), div. (0->0), fcn. (29126->10), ass. (0->80)
t163 = sin(pkin(8));
t165 = cos(pkin(8));
t152 = -t165 * g(1) - t163 * g(2);
t161 = -g(3) + qJDD(1);
t168 = sin(qJ(2));
t171 = cos(qJ(2));
t135 = -t168 * t152 + t171 * t161;
t133 = qJDD(2) * pkin(2) + t135;
t136 = t171 * t152 + t168 * t161;
t172 = qJD(2) ^ 2;
t134 = -t172 * pkin(2) + t136;
t162 = sin(pkin(9));
t164 = cos(pkin(9));
t128 = t164 * t133 - t162 * t134;
t125 = qJDD(2) * pkin(3) + t128;
t129 = t162 * t133 + t164 * t134;
t126 = -t172 * pkin(3) + t129;
t167 = sin(qJ(4));
t170 = cos(qJ(4));
t122 = t167 * t125 + t170 * t126;
t159 = qJD(2) + qJD(4);
t157 = t159 ^ 2;
t158 = qJDD(2) + qJDD(4);
t119 = -t157 * pkin(4) + t158 * pkin(7) + t122;
t151 = t163 * g(1) - t165 * g(2);
t150 = qJDD(3) - t151;
t166 = sin(qJ(5));
t169 = cos(qJ(5));
t116 = -t166 * t119 + t169 * t150;
t117 = t169 * t119 + t166 * t150;
t138 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t166 + Ifges(6,2) * t169) * t159;
t139 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t166 + Ifges(6,4) * t169) * t159;
t186 = qJD(5) * t159;
t143 = t166 * t158 + t169 * t186;
t144 = t169 * t158 - t166 * t186;
t189 = mrSges(6,1) * t116 - mrSges(6,2) * t117 + Ifges(6,5) * t143 + Ifges(6,6) * t144 + Ifges(6,3) * qJDD(5) + (t138 * t166 - t139 * t169) * t159;
t121 = t170 * t125 - t167 * t126;
t118 = -t158 * pkin(4) - t157 * pkin(7) - t121;
t188 = t159 * t166;
t147 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t188;
t187 = t159 * t169;
t148 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t187;
t176 = -m(6) * t118 + t144 * mrSges(6,1) - t143 * mrSges(6,2) - t147 * t188 + t148 * t187;
t108 = m(5) * t121 + t158 * mrSges(5,1) - t157 * mrSges(5,2) + t176;
t142 = (-mrSges(6,1) * t169 + mrSges(6,2) * t166) * t159;
t112 = m(6) * t116 + qJDD(5) * mrSges(6,1) - t143 * mrSges(6,3) + qJD(5) * t148 - t142 * t188;
t113 = m(6) * t117 - qJDD(5) * mrSges(6,2) + t144 * mrSges(6,3) - qJD(5) * t147 + t142 * t187;
t181 = -t166 * t112 + t169 * t113;
t97 = m(5) * t122 - t157 * mrSges(5,1) - t158 * mrSges(5,2) + t181;
t94 = t170 * t108 + t167 * t97;
t91 = m(4) * t128 + qJDD(2) * mrSges(4,1) - t172 * mrSges(4,2) + t94;
t182 = -t167 * t108 + t170 * t97;
t92 = m(4) * t129 - t172 * mrSges(4,1) - qJDD(2) * mrSges(4,2) + t182;
t85 = t162 * t92 + t164 * t91;
t101 = t169 * t112 + t166 * t113;
t185 = m(5) * t150 + t101;
t184 = -t162 * t91 + t164 * t92;
t83 = m(3) * t135 + qJDD(2) * mrSges(3,1) - t172 * mrSges(3,2) + t85;
t84 = m(3) * t136 - t172 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t184;
t183 = -t168 * t83 + t171 * t84;
t179 = m(4) * t150 + t185;
t137 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t166 + Ifges(6,6) * t169) * t159;
t105 = -mrSges(6,1) * t118 + mrSges(6,3) * t117 + Ifges(6,4) * t143 + Ifges(6,2) * t144 + Ifges(6,6) * qJDD(5) + qJD(5) * t139 - t137 * t188;
t106 = mrSges(6,2) * t118 - mrSges(6,3) * t116 + Ifges(6,1) * t143 + Ifges(6,4) * t144 + Ifges(6,5) * qJDD(5) - qJD(5) * t138 + t137 * t187;
t86 = mrSges(5,2) * t150 - mrSges(5,3) * t121 + Ifges(5,5) * t158 - t157 * Ifges(5,6) - pkin(7) * t101 - t166 * t105 + t169 * t106;
t90 = -mrSges(5,1) * t150 + mrSges(5,3) * t122 + t157 * Ifges(5,5) + Ifges(5,6) * t158 - pkin(4) * t101 - t189;
t80 = -mrSges(4,1) * t150 + mrSges(4,3) * t129 + t172 * Ifges(4,5) + Ifges(4,6) * qJDD(2) - pkin(3) * t185 + pkin(6) * t182 + t167 * t86 + t170 * t90;
t81 = mrSges(4,2) * t150 - mrSges(4,3) * t128 + Ifges(4,5) * qJDD(2) - t172 * Ifges(4,6) - pkin(6) * t94 - t167 * t90 + t170 * t86;
t73 = mrSges(3,1) * t151 + mrSges(3,3) * t136 + t172 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t179 + qJ(3) * t184 + t162 * t81 + t164 * t80;
t75 = -mrSges(3,2) * t151 - mrSges(3,3) * t135 + Ifges(3,5) * qJDD(2) - t172 * Ifges(3,6) - qJ(3) * t85 - t162 * t80 + t164 * t81;
t178 = -mrSges(2,2) * t152 + pkin(5) * t183 + t168 * t75 + t171 * t73 + mrSges(2,1) * t151 + pkin(1) * (m(3) * t151 - t179);
t177 = mrSges(5,1) * t121 - mrSges(5,2) * t122 + Ifges(5,3) * t158 + pkin(4) * t176 + pkin(7) * t181 + t169 * t105 + t166 * t106;
t174 = mrSges(4,1) * t128 - mrSges(4,2) * t129 + Ifges(4,3) * qJDD(2) + pkin(3) * t94 + t177;
t173 = mrSges(3,1) * t135 - mrSges(3,2) * t136 + Ifges(3,3) * qJDD(2) + pkin(2) * t85 + t174;
t98 = (m(2) + m(3)) * t151 - t179;
t79 = t168 * t84 + t171 * t83;
t77 = m(2) * t152 + t183;
t76 = -mrSges(2,1) * t161 + mrSges(2,3) * t152 - pkin(1) * t79 - t173;
t71 = mrSges(2,2) * t161 - mrSges(2,3) * t151 - pkin(5) * t79 - t168 * t73 + t171 * t75;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t165 * t71 - t163 * t76 - qJ(1) * (t163 * t77 + t165 * t98), t71, t75, t81, t86, t106; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t163 * t71 + t165 * t76 + qJ(1) * (-t163 * t98 + t165 * t77), t76, t73, t80, t90, t105; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t178, t178, t173, t174, t177, t189;];
m_new = t1;
