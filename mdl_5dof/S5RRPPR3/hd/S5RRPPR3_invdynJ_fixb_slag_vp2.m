% Calculate vector of inverse dynamics joint torques for
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPPR3_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR3_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR3_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR3_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR3_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR3_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR3_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR3_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR3_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:26:28
% EndTime: 2019-12-31 19:26:31
% DurationCPUTime: 1.52s
% Computational Cost: add. (1480->214), mult. (2180->271), div. (0->0), fcn. (1050->12), ass. (0->116)
t94 = sin(qJ(5));
t193 = -t94 / 0.2e1;
t97 = cos(qJ(5));
t174 = t97 / 0.2e1;
t192 = mrSges(4,1) - mrSges(5,2);
t113 = mrSges(6,1) * t94 + mrSges(6,2) * t97;
t191 = t113 + mrSges(5,3);
t190 = -mrSges(6,3) - t192;
t189 = mrSges(4,2) - t191;
t114 = mrSges(6,1) * t97 - mrSges(6,2) * t94;
t159 = Ifges(6,4) * t97;
t160 = Ifges(6,4) * t94;
t142 = pkin(1) * qJD(1);
t95 = sin(qJ(2));
t134 = t95 * t142;
t98 = cos(qJ(2));
t133 = t98 * t142;
t90 = qJD(1) + qJD(2);
t55 = pkin(2) * t90 + t133;
t92 = sin(pkin(8));
t93 = cos(pkin(8));
t26 = t93 * t134 + t55 * t92;
t20 = qJ(4) * t90 + t26;
t188 = ((-Ifges(6,1) * t94 - t159) * t174 + (-Ifges(6,2) * t97 - t160) * t193) * t90 + t20 * t114 + qJD(5) * (-Ifges(6,5) * t94 - Ifges(6,6) * t97) / 0.2e1;
t61 = t92 * t134;
t25 = t55 * t93 - t61;
t121 = qJD(4) - t25;
t172 = -pkin(3) - pkin(7);
t18 = t172 * t90 + t121;
t10 = qJD(3) * t97 + t18 * t94;
t9 = -qJD(3) * t94 + t18 * t97;
t119 = t10 * t94 + t9 * t97;
t150 = t93 * t95;
t169 = pkin(1) * t98;
t171 = pkin(1) * t95;
t149 = t94 * mrSges(6,3);
t53 = -qJD(5) * mrSges(6,2) - t90 * t149;
t146 = t97 * mrSges(6,3);
t54 = qJD(5) * mrSges(6,1) - t90 * t146;
t186 = (m(4) * t25 - m(5) * t121 - m(6) * t119 - t94 * t53 - t97 * t54 + (m(5) * pkin(3) + t192) * t90) * pkin(1) * (t92 * t98 + t150) + (mrSges(3,1) * t171 + mrSges(3,2) * t169) * t90;
t138 = qJD(5) * t94;
t89 = qJDD(1) + qJDD(2);
t46 = -t90 * t138 + t89 * t97;
t137 = qJD(5) * t97;
t47 = -t90 * t137 - t89 * t94;
t184 = -mrSges(6,1) * t47 + mrSges(6,2) * t46 + t89 * mrSges(5,3);
t45 = t113 * t90;
t183 = mrSges(5,3) * t90 + t45;
t27 = qJDD(5) * mrSges(6,1) - mrSges(6,3) * t46;
t28 = -qJDD(5) * mrSges(6,2) + mrSges(6,3) * t47;
t181 = t97 * t27 + t94 * t28;
t40 = t93 * t133 - t61;
t140 = qJD(4) * t90;
t132 = qJD(2) * t171;
t50 = -qJD(1) * t132 + qJDD(1) * t169;
t35 = pkin(2) * t89 + t50;
t141 = qJD(2) * t98;
t51 = (qJD(1) * t141 + qJDD(1) * t95) * pkin(1);
t16 = t35 * t92 + t51 * t93;
t5 = qJ(4) * t89 + t140 + t16;
t167 = pkin(2) * t92;
t74 = qJ(4) + t167;
t180 = t5 * t74 + (-t40 + qJD(4)) * t20;
t91 = qJ(1) + qJ(2);
t85 = pkin(8) + t91;
t76 = sin(t85);
t77 = cos(t85);
t179 = -g(1) * t76 + g(2) * t77;
t86 = sin(t91);
t87 = cos(t91);
t178 = t86 * mrSges(3,1) + t87 * mrSges(3,2) + t189 * t77 + (-m(6) * t172 - t190) * t76;
t177 = -t87 * mrSges(3,1) + t86 * mrSges(3,2) + t189 * t76 + t190 * t77;
t15 = t35 * t93 - t92 * t51;
t117 = qJDD(4) - t15;
t4 = t172 * t89 + t117;
t1 = t9 * qJD(5) + qJDD(3) * t97 + t4 * t94;
t120 = t10 * t97 - t9 * t94;
t2 = -t10 * qJD(5) - qJDD(3) * t94 + t4 * t97;
t101 = t120 * qJD(5) + t1 * t94 + t2 * t97;
t176 = t101 * m(6) + t53 * t137 - t54 * t138 + t181;
t173 = m(4) + m(5);
t96 = sin(qJ(1));
t170 = pkin(1) * t96;
t168 = pkin(2) * t86;
t79 = pkin(2) * t87;
t166 = pkin(2) * t93;
t99 = cos(qJ(1));
t88 = t99 * pkin(1);
t157 = t89 * mrSges(4,1);
t156 = t89 * mrSges(4,2);
t155 = t89 * mrSges(5,2);
t152 = t90 * mrSges(4,2);
t135 = t77 * pkin(3) + t76 * qJ(4) + t79;
t78 = -pkin(3) - t166;
t64 = t77 * qJ(4);
t129 = t64 - t168;
t80 = pkin(2) + t169;
t42 = -t92 * t171 + t80 * t93;
t123 = t88 + t135;
t37 = -pkin(3) - t42;
t122 = -t168 - t170;
t41 = pkin(1) * t93 * t141 - t92 * t132;
t33 = qJD(4) + t41;
t43 = pkin(1) * t150 + t80 * t92;
t36 = qJ(4) + t43;
t118 = t20 * t33 + t36 * t5;
t112 = Ifges(6,1) * t97 - t160;
t111 = -Ifges(6,2) * t94 + t159;
t108 = -pkin(3) * t76 + t129;
t103 = -t20 * t90 + t179;
t31 = Ifges(6,6) * qJD(5) + t111 * t90;
t32 = Ifges(6,5) * qJD(5) + t112 * t90;
t8 = -pkin(3) * t89 + t117;
t100 = -t2 * t146 - t1 * t149 - t31 * t137 / 0.2e1 - t32 * t138 / 0.2e1 + t47 * t111 / 0.2e1 + t46 * t112 / 0.2e1 + (Ifges(6,4) * t46 + Ifges(6,2) * t47) * t193 + t50 * mrSges(3,1) - t51 * mrSges(3,2) + t15 * mrSges(4,1) - t16 * mrSges(4,2) + t8 * mrSges(5,2) + (Ifges(6,1) * t46 + Ifges(6,4) * t47) * t174 + t191 * t5 + (-t10 * t137 + t9 * t138) * mrSges(6,3) + (0.2e1 * Ifges(6,5) * t174 - Ifges(6,6) * t94) * qJDD(5) + (Ifges(5,1) + Ifges(3,3) + Ifges(4,3)) * t89 + t188 * qJD(5);
t70 = t77 * pkin(7);
t3 = [m(4) * (t15 * t42 + t16 * t43 + t26 * t41) + m(3) * (t50 * t98 + t51 * t95) * pkin(1) - t41 * t152 - t43 * t156 + m(5) * (t37 * t8 + t118) + m(6) * t118 + t100 + Ifges(2,3) * qJDD(1) + t37 * t155 + t42 * t157 + (mrSges(3,1) * t169 - mrSges(3,2) * t171) * t89 + t184 * t36 + t183 * t33 + t176 * (-pkin(7) + t37) + (-m(3) * t88 - m(6) * (t70 + t123) - m(5) * t123 - m(4) * (t79 + t88) - mrSges(2,1) * t99 + t96 * mrSges(2,2) + t177) * g(2) + (m(3) * t170 - m(5) * (t108 - t170) - m(6) * (t122 + t64) - m(4) * t122 + t96 * mrSges(2,1) + mrSges(2,2) * t99 + t178) * g(1) - t186 * qJD(2); -t156 * t167 + m(4) * (t15 * t93 + t16 * t92) * pkin(2) + qJD(4) * t45 + t100 + mrSges(5,3) * t140 + t78 * t155 + t157 * t166 + t184 * t74 + t180 * m(6) + (t78 * t8 + t180) * m(5) + t176 * (-pkin(7) + t78) + (-m(4) * t26 + t152 - t183) * t40 + (-m(4) * t79 - m(5) * t135 - m(6) * (t70 + t135) + t177) * g(2) + (m(4) * t168 - m(5) * t108 - m(6) * t129 + t178) * g(1) + t186 * qJD(1); m(6) * (-t119 * qJD(5) + t1 * t97 - t2 * t94) - t53 * t138 + t97 * t28 - t54 * t137 - t94 * t27 + t173 * qJDD(3) + (-m(6) - t173) * g(3); t155 - t183 * t90 + (t97 * t53 - t94 * t54) * qJD(5) + (t101 + t103) * m(6) + (t103 + t8) * m(5) + t181; t2 * mrSges(6,1) - t1 * mrSges(6,2) + Ifges(6,5) * t46 + Ifges(6,6) * t47 + Ifges(6,3) * qJDD(5) + g(3) * t113 + t10 * t54 - t9 * t53 + (t94 * t32 / 0.2e1 + t31 * t174 + t120 * mrSges(6,3) - t188) * t90 + t179 * t114;];
tau = t3;
