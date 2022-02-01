% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% m [6x1]
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
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:00
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPPPR2_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR2_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR2_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR2_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPPR2_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 08:59:06
% EndTime: 2022-01-23 08:59:08
% DurationCPUTime: 1.04s
% Computational Cost: add. (3435->199), mult. (8207->314), div. (0->0), fcn. (8565->8), ass. (0->122)
t108 = sin(pkin(9));
t111 = cos(pkin(9));
t113 = cos(pkin(7));
t134 = t113 * t111;
t110 = sin(pkin(7));
t112 = cos(pkin(8));
t137 = t112 * t110;
t87 = t108 * t137 + t134;
t162 = t87 * mrSges(6,2);
t114 = sin(qJ(5));
t109 = sin(pkin(8));
t115 = cos(qJ(5));
t139 = t109 * t115;
t135 = t113 * t108;
t88 = t111 * t137 - t135;
t67 = t110 * t139 - t114 * t88;
t165 = t67 * mrSges(6,3);
t46 = -t162 + t165;
t180 = -t46 / 0.2e1;
t179 = t67 / 0.2e1;
t178 = t108 ^ 2;
t177 = t109 ^ 2;
t106 = t110 ^ 2;
t176 = t112 ^ 2;
t175 = m(4) / 0.4e1;
t174 = m(5) / 0.2e1;
t173 = m(6) / 0.2e1;
t172 = -mrSges(6,1) / 0.2e1;
t171 = mrSges(6,2) / 0.2e1;
t170 = t87 / 0.2e1;
t91 = t111 * t139 - t114 * t112;
t169 = -t91 / 0.2e1;
t168 = t110 / 0.2e1;
t167 = -t114 / 0.2e1;
t166 = -t115 / 0.2e1;
t140 = t109 * t114;
t68 = t110 * t140 + t88 * t115;
t164 = t68 * mrSges(6,3);
t163 = t87 * mrSges(6,1);
t142 = t109 * t110;
t161 = -mrSges(5,1) * t142 - t67 * mrSges(6,1) + t68 * mrSges(6,2) + t88 * mrSges(5,3);
t160 = t113 * mrSges(4,1) + t87 * mrSges(5,1) + t88 * mrSges(5,2) + mrSges(4,3) * t137;
t159 = Ifges(6,5) * t67 - Ifges(6,6) * t68;
t136 = t112 * t113;
t95 = -pkin(2) * t113 - t110 * qJ(3) - pkin(1);
t79 = qJ(2) * t136 + t109 * t95;
t71 = -t113 * qJ(4) + t79;
t84 = (pkin(3) * t109 - qJ(4) * t112 + qJ(2)) * t110;
t40 = t108 * t84 + t111 * t71;
t141 = t109 * t113;
t78 = -qJ(2) * t141 + t112 * t95;
t72 = t113 * pkin(3) - t78;
t32 = t87 * pkin(4) - t88 * pkin(6) + t72;
t34 = pkin(6) * t142 + t40;
t20 = -t114 * t34 + t115 * t32;
t21 = t114 * t32 + t115 * t34;
t39 = -t108 * t71 + t111 * t84;
t33 = -pkin(4) * t142 - t39;
t35 = t68 * mrSges(6,1) + t67 * mrSges(6,2);
t47 = t163 - t164;
t1 = t20 * t46 - t21 * t47 + t33 * t35 + t159 * t170 + (-t20 * mrSges(6,3) + Ifges(6,4) * t67 + Ifges(6,5) * t170) * t67 + (-t21 * mrSges(6,3) - Ifges(6,6) * t87 / 0.2e1 - Ifges(6,4) * t68 + (Ifges(6,1) - Ifges(6,2)) * t67) * t68;
t158 = t1 * qJD(1);
t138 = t111 * t110;
t86 = t112 * t135 - t138;
t157 = t108 * t86;
t156 = t108 * t87;
t155 = t111 * t86;
t154 = t111 * t88;
t153 = t114 * mrSges(6,1);
t152 = t114 * t47;
t151 = t115 * mrSges(6,2);
t150 = t115 * t46;
t104 = t106 * qJ(2);
t107 = t113 ^ 2;
t89 = t110 * t108 + t112 * t134;
t69 = t113 * t139 - t89 * t114;
t70 = t113 * t140 + t89 * t115;
t73 = -mrSges(5,2) * t142 - t87 * mrSges(5,3);
t93 = t113 * mrSges(4,2) - mrSges(4,3) * t142;
t2 = t70 * t46 + t69 * t47 + t89 * t73 + t161 * t86 + (mrSges(4,1) * t109 + mrSges(4,2) * t112) * t106 + (t107 + t106) * mrSges(3,3) + (t160 * t109 + t112 * t93) * t113 + m(6) * (t20 * t69 + t21 * t70 + t33 * t86) + m(5) * (t72 * t141 - t39 * t86 + t40 * t89) + m(4) * (t104 + (-t109 * t78 + t112 * t79) * t113) + m(3) * (t107 * qJ(2) + t104);
t149 = t2 * qJD(1);
t143 = t108 * t109;
t123 = t111 * t140 + t115 * t112;
t82 = t123 * t110;
t83 = t91 * t110;
t3 = m(6) * (t20 * t82 - t21 * t83) - t83 * t46 + t82 * t47 + (t160 * t112 + (-t161 * t108 - t111 * t73 - t93) * t109 - m(6) * t33 * t143 + m(5) * (-t109 * t111 * t40 + t112 * t72 + t39 * t143) + m(4) * (-t109 * t79 - t112 * t78)) * t110;
t148 = t3 * qJD(1);
t120 = (t123 * t179 + t68 * t169) * mrSges(6,3) + t123 * t180 + t47 * t169 + t35 * t143 / 0.2e1;
t125 = t70 * t171 + t69 * t172;
t5 = t120 + t125;
t147 = t5 * qJD(1);
t6 = t161 * t88 + (-t150 - t73 + t152) * t87 + m(6) * (t33 * t88 + (t114 * t20 - t115 * t21) * t87) + m(5) * (-t39 * t88 - t40 * t87);
t146 = t6 * qJD(1);
t118 = (t46 * t167 + t47 * t166 + (t114 * t179 + t68 * t166) * mrSges(6,3)) * t108 - t111 * t35 / 0.2e1;
t124 = -t83 * t171 + t82 * t172;
t9 = t118 + t124;
t145 = t9 * qJD(1);
t10 = (t162 / 0.2e1 + t180 + t165 / 0.2e1) * t115 + (t163 / 0.2e1 + t164 / 0.2e1 + t47 / 0.2e1) * t114;
t144 = t10 * qJD(1);
t130 = m(4) * t168;
t116 = t130 + (t108 * t89 - t155) * t174 + (-t155 + (-t114 * t69 + t115 * t70) * t108) * t173;
t126 = -t123 * t82 - t91 * t83;
t127 = -t111 ^ 2 - t178;
t12 = -m(6) * t126 / 0.2e1 + 0.2e1 * ((m(5) / 0.4e1 + t175) * t176 + (m(6) * t178 / 0.4e1 - m(5) * t127 / 0.4e1 + t175) * t177) * t110 + t116;
t133 = t12 * qJD(1);
t129 = t109 * t174;
t117 = (t88 * t143 + (-t114 * t123 - t115 * t91) * t87) * t173 + (t108 * t88 - t111 * t87) * t129;
t122 = (t114 * t70 + t115 * t69) * t173 + t113 * t129;
t17 = t117 - t122;
t132 = t17 * qJD(1);
t119 = (-t154 + (-t114 ^ 2 - t115 ^ 2) * t156) * t173 + (-t154 - t156) * t174;
t128 = m(5) * t168;
t121 = (-t114 * t83 + t115 * t82) * t173 + t112 * t128;
t19 = t119 - t121;
t131 = t19 * qJD(1);
t18 = t119 + t121;
t16 = t117 + t122;
t13 = (-t178 * t177 * t110 + t126) * t173 + (t127 * t177 - t176) * t128 + (-t176 - t177) * t130 + t116;
t11 = t150 / 0.2e1 + t164 * t167 - t152 / 0.2e1 + t165 * t166 + (t153 / 0.2e1 + t151 / 0.2e1) * t87;
t8 = t118 - t124;
t4 = t120 - t125;
t7 = [t2 * qJD(2) + t3 * qJD(3) + t6 * qJD(4) + t1 * qJD(5), t149 + 0.2e1 * ((-t123 * t69 + t91 * t70) * t173 + (t157 * t173 + (t111 * t89 - t136 + t157) * t174) * t109) * qJD(2) + t13 * qJD(3) + t16 * qJD(4) + t4 * qJD(5), t148 + t13 * qJD(2) + t18 * qJD(4) + t8 * qJD(5) + m(6) * (t109 * t138 - t114 * t82 - t115 * t83) * qJD(3) * t108, t16 * qJD(2) + t18 * qJD(3) + t11 * qJD(5) + t146, t158 + t4 * qJD(2) + t8 * qJD(3) + t11 * qJD(4) + (-t21 * mrSges(6,1) - t20 * mrSges(6,2) + t159) * qJD(5); -t12 * qJD(3) + t17 * qJD(4) + t5 * qJD(5) - t149, 0, -t133, t132, t147 + (-t91 * mrSges(6,1) + mrSges(6,2) * t123) * qJD(5); t12 * qJD(2) + t19 * qJD(4) + t9 * qJD(5) - t148, t133, 0, t131, t145 + (-mrSges(6,1) * t115 + mrSges(6,2) * t114) * qJD(5) * t108; -t17 * qJD(2) - t19 * qJD(3) - t10 * qJD(5) - t146, -t132, -t131, 0, -t144 + (-t151 - t153) * qJD(5); -t5 * qJD(2) - t9 * qJD(3) + t10 * qJD(4) - t158, -t147, -t145, t144, 0;];
Cq = t7;
