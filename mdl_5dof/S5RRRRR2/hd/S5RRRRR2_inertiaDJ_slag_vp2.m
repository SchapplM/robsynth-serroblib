% Calculate time derivative of joint inertia matrix for
% S5RRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
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
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:54
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRR2_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(2,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR2_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR2_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5RRRRR2_inertiaDJ_slag_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR2_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR2_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR2_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:52:57
% EndTime: 2019-12-05 18:53:01
% DurationCPUTime: 1.24s
% Computational Cost: add. (1577->199), mult. (4632->318), div. (0->0), fcn. (3908->8), ass. (0->109)
t174 = 2 * pkin(2);
t92 = cos(qJ(3));
t173 = 0.2e1 * t92;
t86 = sin(qJ(5));
t90 = cos(qJ(5));
t172 = -mrSges(6,1) * t90 + mrSges(6,2) * t86 - mrSges(5,1);
t167 = qJD(3) + qJD(4);
t89 = sin(qJ(2));
t171 = t167 * t89;
t127 = qJD(5) * t90;
t87 = sin(qJ(4));
t88 = sin(qJ(3));
t138 = t87 * t88;
t91 = cos(qJ(4));
t63 = -t91 * t92 + t138;
t45 = t167 * t63;
t139 = t86 * t45;
t64 = t87 * t92 + t88 * t91;
t104 = t64 * t127 - t139;
t128 = qJD(5) * t86;
t146 = t45 * t90;
t103 = t64 * t128 + t146;
t93 = cos(qJ(2));
t134 = qJD(2) * t93;
t20 = (-t63 * t134 - t64 * t171) * pkin(1);
t162 = pkin(1) * t89;
t54 = t63 * t162;
t78 = -pkin(1) * t93 - t92 * pkin(2);
t34 = t54 * t86 + t78 * t90;
t133 = qJD(3) * t88;
t70 = pkin(2) * t133 + qJD(2) * t162;
t8 = qJD(5) * t34 + t20 * t90 + t70 * t86;
t35 = -t54 * t90 + t78 * t86;
t9 = -qJD(5) * t35 - t20 * t86 + t70 * t90;
t170 = t8 * t90 - t9 * t86;
t135 = t88 ^ 2 + t92 ^ 2;
t169 = t135 * mrSges(4,3) - mrSges(3,2);
t136 = t86 ^ 2 + t90 ^ 2;
t168 = t136 * mrSges(6,3) - mrSges(5,2);
t126 = qJD(5) * t92;
t46 = t167 * t64;
t12 = mrSges(6,1) * t46 + t103 * mrSges(6,3);
t13 = -mrSges(6,2) * t46 - t104 * mrSges(6,3);
t22 = mrSges(5,1) * t46 - mrSges(5,2) * t45;
t166 = -t12 * t90 - t13 * t86 - t22;
t129 = qJD(4) * t91;
t121 = t90 * t129;
t161 = pkin(2) * t87;
t125 = qJD(4) * t161;
t132 = qJD(3) * t92;
t110 = mrSges(6,1) * t86 + mrSges(6,2) * t90;
t37 = t110 * t64;
t144 = t64 * t86;
t40 = -mrSges(6,2) * t63 - mrSges(6,3) * t144;
t165 = pkin(2) * t40 * t121 + t90 * t13 * t161 + Ifges(4,5) * t132 - Ifges(4,6) * t133 + t37 * t125 + (-pkin(2) * t129 * t63 + t125 * t64 - t161 * t46) * mrSges(5,3);
t164 = 2 * m(6);
t160 = t8 * t40;
t143 = t64 * t90;
t41 = mrSges(6,1) * t63 - mrSges(6,3) * t143;
t158 = t9 * t41;
t156 = Ifges(4,4) * t88;
t154 = Ifges(6,4) * t86;
t153 = Ifges(6,4) * t90;
t122 = t92 * t134;
t123 = t88 * t134;
t21 = (t87 * t122 + (t92 * t171 + t123) * t91) * pkin(1) - t167 * t138 * t162;
t150 = t21 * t37;
t53 = t64 * t162;
t149 = t21 * t53;
t148 = t34 * t12;
t147 = t35 * t13;
t10 = t104 * mrSges(6,1) - t103 * mrSges(6,2);
t145 = t53 * t10;
t65 = t110 * qJD(5);
t142 = t65 * t91;
t49 = mrSges(5,1) * t63 + mrSges(5,2) * t64;
t141 = t70 * t49;
t140 = t78 * t22;
t137 = t93 * (mrSges(4,1) * t88 + mrSges(4,2) * t92) * qJD(3);
t131 = qJD(4) * t53;
t130 = qJD(4) * t86;
t119 = t87 * t127;
t117 = (-mrSges(4,1) * t92 + mrSges(4,2) * t88 - mrSges(3,1)) * t89;
t116 = t172 * t87;
t112 = mrSges(5,3) * t45 - t10;
t67 = Ifges(6,5) * t127 - Ifges(6,6) * t128;
t109 = Ifges(6,1) * t90 - t154;
t108 = -Ifges(6,2) * t86 + t153;
t107 = t34 * t90 + t35 * t86;
t106 = -t86 * t40 - t90 * t41;
t105 = -t106 + t49;
t68 = t108 * qJD(5);
t69 = t109 * qJD(5);
t75 = Ifges(6,2) * t90 + t154;
t77 = Ifges(6,1) * t86 + t153;
t102 = t77 * t127 - t75 * t128 + t90 * t68 + t86 * t69;
t101 = (-t20 * t63 + t21 * t64 - t45 * t53 + t46 * t54) * mrSges(5,3);
t27 = Ifges(6,6) * t63 + t108 * t64;
t28 = Ifges(6,5) * t63 + t109 * t64;
t5 = -t103 * Ifges(6,4) - t104 * Ifges(6,2) + Ifges(6,6) * t46;
t6 = -t103 * Ifges(6,1) - t104 * Ifges(6,4) + Ifges(6,5) * t46;
t100 = -t68 * t144 / 0.2e1 - t77 * t146 / 0.2e1 + t69 * t143 / 0.2e1 + t46 * (Ifges(6,5) * t86 + Ifges(6,6) * t90) / 0.2e1 + t63 * t67 / 0.2e1 + t28 * t127 / 0.2e1 + t86 * t6 / 0.2e1 + t90 * t5 / 0.2e1 - Ifges(5,5) * t45 - Ifges(5,6) * t46 - t104 * t75 / 0.2e1 - (t64 * t77 + t27) * t128 / 0.2e1;
t99 = -t103 * Ifges(6,5) - t104 * Ifges(6,6) + Ifges(6,3) * t46;
t98 = t6 * t143 - t28 * t146 + t63 * t99 + (Ifges(4,1) * t92 - t156) * t133 + t27 * t139 + (Ifges(4,4) * t173 + (Ifges(4,1) - Ifges(4,2)) * t88) * t132 + (-0.2e1 * Ifges(5,1) * t45 - t86 * t5 + (-t27 * t90 - t28 * t86) * qJD(5)) * t64 + ((Ifges(6,5) * t90 - Ifges(6,6) * t86) * t64 + (Ifges(6,3) + (2 * Ifges(5,2))) * t63) * t46 + 0.2e1 * (t45 * t63 - t64 * t46) * Ifges(5,4);
t76 = Ifges(4,2) * t92 + t156;
t97 = -t76 * t133 + t98;
t96 = -t20 * mrSges(5,2) + t53 * t65 + t100 + t172 * t21 + (-t107 * qJD(5) + t170) * mrSges(6,3);
t94 = pkin(2) ^ 2;
t1 = [0.2e1 * m(5) * (-t20 * t54 + t78 * t70 + t149) + (t34 * t9 + t35 * t8 + t149) * t164 + 0.2e1 * t141 + 0.2e1 * t140 + 0.2e1 * t145 + 0.2e1 * t148 + 0.2e1 * t147 + 0.2e1 * t150 + 0.2e1 * t160 + 0.2e1 * t158 + t97 + 0.2e1 * t101 + (0.2e1 * (t117 + (m(4) * (-0.1e1 + t135) * t162 + t169) * t93) * qJD(2) - 0.2e1 * t137) * pkin(1); ((m(5) * t78 + m(6) * t107 + t105) * t133 + (m(6) * (-t35 * t127 + t34 * t128 - t8 * t86 - t9 * t90) - m(5) * t70 - t40 * t127 + t41 * t128 + t166) * t92) * pkin(2) + (-t137 + (t169 * t93 + t117) * qJD(2)) * pkin(1) + t141 + t140 + t145 + t148 + t147 + t150 + t160 + t158 + t97 + t101; 0.2e1 * (t166 * t92 + t105 * t133 + (-t40 * t90 + t86 * t41) * t126) * pkin(2) + t98 + (-t76 + (-m(6) * t136 - m(5)) * t94 * t173) * t133; ((-t41 * t130 + m(5) * (-qJD(4) * t54 - t21) + m(6) * (qJD(4) * t35 * t90 - t34 * t130 - t21) + t112) * t91 + (-t86 * t12 + t106 * qJD(5) + m(5) * (t20 + t131) + m(6) * (-t34 * t127 - t35 * t128 + t131 + t170)) * t87) * pkin(2) + t96 + ((t89 * t133 - t122) * mrSges(4,2) + (-t89 * t132 - t123) * mrSges(4,1)) * pkin(1) + t165; (-t41 * t119 + t112 * t91 + (-t41 * t129 + (-qJD(5) * t40 - t12) * t87) * t86) * pkin(2) + t100 + t165; -0.2e1 * pkin(2) * t142 + (t116 * t174 + ((-0.1e1 + t136) * t94 * t87 * t164 + t168 * t174) * t91) * qJD(4) + t102; t96; t100; (-t142 + (t168 * t91 + t116) * qJD(4)) * pkin(2) + t102; t102; mrSges(6,1) * t9 - mrSges(6,2) * t8 + t99; ((t90 * t126 - t86 * t133) * mrSges(6,2) + (t86 * t126 + t90 * t133) * mrSges(6,1)) * pkin(2) + t99; ((t87 * t128 - t121) * mrSges(6,2) + (-t86 * t129 - t119) * mrSges(6,1)) * pkin(2) + t67; t67; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
