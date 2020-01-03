% Calculate vector of cutting torques with Newton-Euler for
% S5PRRPR8
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
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-31 17:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5PRRPR8_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR8_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR8_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR8_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR8_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR8_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR8_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR8_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR8_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:42:23
% EndTime: 2019-12-31 17:42:27
% DurationCPUTime: 2.36s
% Computational Cost: add. (38513->171), mult. (48237->219), div. (0->0), fcn. (29126->10), ass. (0->80)
t163 = sin(pkin(8));
t187 = cos(pkin(8));
t152 = -t187 * g(1) - t163 * g(2);
t161 = -g(3) + qJDD(1);
t167 = sin(qJ(2));
t170 = cos(qJ(2));
t135 = -t167 * t152 + t170 * t161;
t133 = qJDD(2) * pkin(2) + t135;
t136 = t170 * t152 + t167 * t161;
t171 = qJD(2) ^ 2;
t134 = -t171 * pkin(2) + t136;
t166 = sin(qJ(3));
t169 = cos(qJ(3));
t128 = t169 * t133 - t166 * t134;
t159 = qJDD(2) + qJDD(3);
t125 = t159 * pkin(3) + t128;
t129 = t166 * t133 + t169 * t134;
t160 = qJD(2) + qJD(3);
t158 = t160 ^ 2;
t126 = -t158 * pkin(3) + t129;
t162 = sin(pkin(9));
t164 = cos(pkin(9));
t122 = t162 * t125 + t164 * t126;
t119 = -t158 * pkin(4) + t159 * pkin(7) + t122;
t151 = t163 * g(1) - t187 * g(2);
t150 = qJDD(4) - t151;
t165 = sin(qJ(5));
t168 = cos(qJ(5));
t116 = -t165 * t119 + t168 * t150;
t117 = t168 * t119 + t165 * t150;
t138 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t165 + Ifges(6,2) * t168) * t160;
t139 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t165 + Ifges(6,4) * t168) * t160;
t184 = qJD(5) * t160;
t143 = t165 * t159 + t168 * t184;
t144 = t168 * t159 - t165 * t184;
t189 = mrSges(6,1) * t116 - mrSges(6,2) * t117 + Ifges(6,5) * t143 + Ifges(6,6) * t144 + Ifges(6,3) * qJDD(5) + (t138 * t165 - t139 * t168) * t160;
t188 = m(3) + m(4);
t121 = t164 * t125 - t162 * t126;
t118 = -t159 * pkin(4) - t158 * pkin(7) - t121;
t186 = t160 * t165;
t147 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t186;
t185 = t160 * t168;
t148 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t185;
t175 = -m(6) * t118 + t144 * mrSges(6,1) - t143 * mrSges(6,2) - t147 * t186 + t148 * t185;
t108 = m(5) * t121 + t159 * mrSges(5,1) - t158 * mrSges(5,2) + t175;
t142 = (-mrSges(6,1) * t168 + mrSges(6,2) * t165) * t160;
t112 = m(6) * t116 + qJDD(5) * mrSges(6,1) - t143 * mrSges(6,3) + qJD(5) * t148 - t142 * t186;
t113 = m(6) * t117 - qJDD(5) * mrSges(6,2) + t144 * mrSges(6,3) - qJD(5) * t147 + t142 * t185;
t179 = -t165 * t112 + t168 * t113;
t97 = m(5) * t122 - t158 * mrSges(5,1) - t159 * mrSges(5,2) + t179;
t94 = t164 * t108 + t162 * t97;
t90 = m(4) * t128 + t159 * mrSges(4,1) - t158 * mrSges(4,2) + t94;
t180 = -t162 * t108 + t164 * t97;
t91 = m(4) * t129 - t158 * mrSges(4,1) - t159 * mrSges(4,2) + t180;
t85 = t166 * t91 + t169 * t90;
t101 = t168 * t112 + t165 * t113;
t183 = m(5) * t150 + t101;
t182 = -t166 * t90 + t169 * t91;
t83 = m(3) * t135 + qJDD(2) * mrSges(3,1) - t171 * mrSges(3,2) + t85;
t84 = m(3) * t136 - t171 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t182;
t181 = -t167 * t83 + t170 * t84;
t137 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t165 + Ifges(6,6) * t168) * t160;
t105 = -mrSges(6,1) * t118 + mrSges(6,3) * t117 + Ifges(6,4) * t143 + Ifges(6,2) * t144 + Ifges(6,6) * qJDD(5) + qJD(5) * t139 - t137 * t186;
t106 = mrSges(6,2) * t118 - mrSges(6,3) * t116 + Ifges(6,1) * t143 + Ifges(6,4) * t144 + Ifges(6,5) * qJDD(5) - qJD(5) * t138 + t137 * t185;
t86 = mrSges(5,2) * t150 - mrSges(5,3) * t121 + Ifges(5,5) * t159 - t158 * Ifges(5,6) - pkin(7) * t101 - t165 * t105 + t168 * t106;
t92 = -mrSges(5,1) * t150 + mrSges(5,3) * t122 + t158 * Ifges(5,5) + Ifges(5,6) * t159 - pkin(4) * t101 - t189;
t80 = mrSges(4,1) * t151 + mrSges(4,3) * t129 + t158 * Ifges(4,5) + Ifges(4,6) * t159 - pkin(3) * t183 + qJ(4) * t180 + t162 * t86 + t164 * t92;
t81 = -mrSges(4,2) * t151 - mrSges(4,3) * t128 + Ifges(4,5) * t159 - t158 * Ifges(4,6) - qJ(4) * t94 - t162 * t92 + t164 * t86;
t73 = Ifges(3,6) * qJDD(2) + t171 * Ifges(3,5) + mrSges(3,1) * t151 + mrSges(3,3) * t136 + t166 * t81 + t169 * t80 - pkin(2) * (-m(4) * t151 + t183) + pkin(6) * t182;
t75 = -mrSges(3,2) * t151 - mrSges(3,3) * t135 + Ifges(3,5) * qJDD(2) - t171 * Ifges(3,6) - pkin(6) * t85 - t166 * t80 + t169 * t81;
t177 = -mrSges(2,2) * t152 + pkin(5) * t181 + t167 * t75 + t170 * t73 + mrSges(2,1) * t151 + pkin(1) * (t188 * t151 - t183);
t176 = mrSges(5,1) * t121 - mrSges(5,2) * t122 + Ifges(5,3) * t159 + pkin(4) * t175 + pkin(7) * t179 + t168 * t105 + t165 * t106;
t173 = mrSges(4,1) * t128 - mrSges(4,2) * t129 + Ifges(4,3) * t159 + pkin(3) * t94 + t176;
t172 = mrSges(3,1) * t135 - mrSges(3,2) * t136 + Ifges(3,3) * qJDD(2) + pkin(2) * t85 + t173;
t98 = (m(2) + t188) * t151 - t183;
t79 = t167 * t84 + t170 * t83;
t77 = m(2) * t152 + t181;
t76 = -mrSges(2,1) * t161 + mrSges(2,3) * t152 - pkin(1) * t79 - t172;
t71 = mrSges(2,2) * t161 - mrSges(2,3) * t151 - pkin(5) * t79 - t167 * t73 + t170 * t75;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t187 * t71 - t163 * t76 - qJ(1) * (t163 * t77 + t187 * t98), t71, t75, t81, t86, t106; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t163 * t71 + t187 * t76 + qJ(1) * (-t163 * t98 + t187 * t77), t76, t73, t80, t92, t105; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t177, t177, t172, t173, t176, t189;];
m_new = t1;
