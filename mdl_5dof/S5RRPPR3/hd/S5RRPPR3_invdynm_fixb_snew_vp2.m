% Calculate vector of cutting torques with Newton-Euler for
% S5RRPPR3
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
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
% Datum: 2019-12-31 19:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRPPR3_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR3_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR3_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR3_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR3_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR3_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR3_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR3_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR3_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:26:31
% EndTime: 2019-12-31 19:26:32
% DurationCPUTime: 1.38s
% Computational Cost: add. (22374->180), mult. (27659->219), div. (0->0), fcn. (13400->8), ass. (0->81)
t190 = -pkin(3) - pkin(7);
t189 = mrSges(4,1) - mrSges(5,2);
t188 = -Ifges(5,4) + Ifges(4,5);
t187 = Ifges(5,5) - Ifges(4,6);
t164 = sin(qJ(2));
t167 = cos(qJ(2));
t165 = sin(qJ(1));
t168 = cos(qJ(1));
t144 = t165 * g(1) - g(2) * t168;
t139 = qJDD(1) * pkin(1) + t144;
t145 = -g(1) * t168 - g(2) * t165;
t169 = qJD(1) ^ 2;
t140 = -pkin(1) * t169 + t145;
t121 = t167 * t139 - t140 * t164;
t157 = qJD(1) + qJD(2);
t155 = t157 ^ 2;
t156 = qJDD(1) + qJDD(2);
t161 = sin(pkin(8));
t162 = cos(pkin(8));
t118 = pkin(2) * t156 + t121;
t122 = t164 * t139 + t167 * t140;
t119 = -pkin(2) * t155 + t122;
t113 = t162 * t118 - t161 * t119;
t179 = -t155 * qJ(4) + qJDD(4) - t113;
t111 = -pkin(3) * t156 + t179;
t108 = t190 * t156 + t179;
t160 = -g(3) + qJDD(3);
t163 = sin(qJ(5));
t166 = cos(qJ(5));
t104 = t108 * t166 - t160 * t163;
t133 = (mrSges(6,1) * t163 + mrSges(6,2) * t166) * t157;
t184 = qJD(5) * t157;
t135 = t156 * t166 - t163 * t184;
t186 = t157 * t163;
t141 = -qJD(5) * mrSges(6,2) - mrSges(6,3) * t186;
t185 = t157 * t166;
t101 = m(6) * t104 + qJDD(5) * mrSges(6,1) - mrSges(6,3) * t135 + qJD(5) * t141 - t133 * t185;
t105 = t108 * t163 + t160 * t166;
t134 = -t156 * t163 - t166 * t184;
t142 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t185;
t102 = m(6) * t105 - qJDD(5) * mrSges(6,2) + t134 * mrSges(6,3) - qJD(5) * t142 - t133 * t186;
t93 = t166 * t101 + t163 * t102;
t178 = -m(5) * t111 + t155 * mrSges(5,3) - t93;
t87 = m(4) * t113 - t155 * mrSges(4,2) + t189 * t156 + t178;
t114 = t161 * t118 + t162 * t119;
t180 = t156 * qJ(4) + 0.2e1 * qJD(4) * t157 + t114;
t109 = pkin(3) * t155 - t180;
t107 = t190 * t155 + t180;
t99 = -m(6) * t107 + t134 * mrSges(6,1) - t135 * mrSges(6,2) - t141 * t186 - t142 * t185;
t175 = -m(5) * t109 + t155 * mrSges(5,2) + t156 * mrSges(5,3) - t99;
t91 = m(4) * t114 - mrSges(4,1) * t155 - mrSges(4,2) * t156 + t175;
t85 = t161 * t91 + t162 * t87;
t82 = m(3) * t121 + mrSges(3,1) * t156 - mrSges(3,2) * t155 + t85;
t183 = -t161 * t87 + t162 * t91;
t83 = m(3) * t122 - mrSges(3,1) * t155 - mrSges(3,2) * t156 + t183;
t76 = t164 * t83 + t167 * t82;
t182 = -t164 * t82 + t167 * t83;
t94 = -t101 * t163 + t166 * t102;
t92 = m(5) * t160 + t94;
t181 = m(4) * t160 + t92;
t126 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t166 - Ifges(6,2) * t163) * t157;
t127 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t166 - Ifges(6,4) * t163) * t157;
t177 = mrSges(6,1) * t104 - mrSges(6,2) * t105 + Ifges(6,5) * t135 + Ifges(6,6) * t134 + Ifges(6,3) * qJDD(5) + t126 * t185 + t127 * t186;
t125 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t166 - Ifges(6,6) * t163) * t157;
t96 = -mrSges(6,1) * t107 + mrSges(6,3) * t105 + Ifges(6,4) * t135 + Ifges(6,2) * t134 + Ifges(6,6) * qJDD(5) + qJD(5) * t127 - t125 * t185;
t97 = mrSges(6,2) * t107 - mrSges(6,3) * t104 + Ifges(6,1) * t135 + Ifges(6,4) * t134 + Ifges(6,5) * qJDD(5) - qJD(5) * t126 - t125 * t186;
t176 = mrSges(5,2) * t111 - mrSges(5,3) * t109 + Ifges(5,1) * t156 - pkin(7) * t93 - t163 * t96 + t166 * t97;
t174 = -mrSges(5,1) * t109 - pkin(4) * t99 - pkin(7) * t94 - t163 * t97 - t166 * t96;
t173 = mrSges(5,1) * t111 + pkin(4) * t93 + t177;
t172 = -mrSges(4,2) * t114 + mrSges(4,1) * t113 + Ifges(4,3) * t156 + t176 + pkin(3) * (-mrSges(5,2) * t156 + t178) + qJ(4) * t175;
t171 = mrSges(3,1) * t121 - mrSges(3,2) * t122 + Ifges(3,3) * t156 + pkin(2) * t85 + t172;
t170 = mrSges(2,1) * t144 - mrSges(2,2) * t145 + Ifges(2,3) * qJDD(1) + pkin(1) * t76 + t171;
t78 = -mrSges(4,3) * t113 - qJ(4) * t92 + (mrSges(4,2) - mrSges(5,3)) * t160 + t188 * t156 + t187 * t155 + t173;
t77 = mrSges(4,3) * t114 - pkin(3) * t92 + t188 * t155 - t187 * t156 - t189 * t160 + t174;
t74 = m(2) * t145 - mrSges(2,1) * t169 - qJDD(1) * mrSges(2,2) + t182;
t73 = m(2) * t144 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t169 + t76;
t72 = -mrSges(3,2) * g(3) - mrSges(3,3) * t121 + Ifges(3,5) * t156 - Ifges(3,6) * t155 - qJ(3) * t85 - t161 * t77 + t162 * t78;
t71 = mrSges(3,1) * g(3) + mrSges(3,3) * t122 + t155 * Ifges(3,5) + Ifges(3,6) * t156 - pkin(2) * t181 + qJ(3) * t183 + t161 * t78 + t162 * t77;
t70 = -mrSges(2,2) * g(3) - mrSges(2,3) * t144 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t169 - pkin(6) * t76 - t164 * t71 + t167 * t72;
t69 = Ifges(2,6) * qJDD(1) + t169 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t145 + t164 * t72 + t167 * t71 - pkin(1) * (-m(3) * g(3) + t181) + pkin(6) * t182;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t168 * t70 - t165 * t69 - pkin(5) * (t165 * t74 + t168 * t73), t70, t72, t78, t176, t97; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t165 * t70 + t168 * t69 + pkin(5) * (-t165 * t73 + t168 * t74), t69, t71, t77, mrSges(5,3) * t160 + Ifges(5,4) * t156 - Ifges(5,5) * t155 - t173, t96; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t170, t170, t171, t172, -mrSges(5,2) * t160 + t155 * Ifges(5,4) + Ifges(5,5) * t156 - t174, t177;];
m_new = t1;
