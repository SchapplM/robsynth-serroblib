% Calculate vector of inverse dynamics joint torques for
% S5RRRPR2
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-05 18:41
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRPR2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR2_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR2_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR2_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR2_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR2_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR2_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR2_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:40:44
% EndTime: 2019-12-05 18:40:49
% DurationCPUTime: 1.96s
% Computational Cost: add. (3151->278), mult. (5136->374), div. (0->0), fcn. (2732->16), ass. (0->139)
t118 = sin(qJ(5));
t187 = t118 / 0.2e1;
t116 = sin(pkin(9));
t117 = cos(pkin(9));
t123 = cos(qJ(3));
t119 = sin(qJ(3));
t164 = t116 * t119;
t173 = pkin(2) * qJD(3);
t124 = cos(qJ(2));
t120 = sin(qJ(2));
t161 = t120 * t123;
t140 = -t119 * t124 - t161;
t172 = qJD(1) * pkin(1);
t70 = t140 * t172;
t162 = t119 * t120;
t139 = t123 * t124 - t162;
t71 = t139 * t172;
t178 = t116 * t70 + t117 * t71 - (t117 * t123 - t164) * t173;
t122 = cos(qJ(5));
t145 = mrSges(6,1) * t118 + mrSges(6,2) * t122;
t114 = qJD(1) + qJD(2);
t109 = qJD(3) + t114;
t154 = t120 * t172;
t153 = t124 * t172;
t84 = pkin(2) * t114 + t153;
t50 = t119 * t84 + t123 * t154;
t171 = t116 * t50;
t49 = -t119 * t154 + t123 * t84;
t42 = pkin(3) * t109 + t49;
t23 = t117 * t42 - t171;
t18 = -pkin(4) * t109 - t23;
t191 = t18 * t145 + qJD(5) * (Ifges(6,5) * t122 - Ifges(6,6) * t118) / 0.2e1;
t163 = t117 * t119;
t179 = t116 * t71 - t117 * t70 - (t116 * t123 + t163) * t173;
t115 = qJ(1) + qJ(2);
t110 = sin(t115);
t111 = cos(t115);
t188 = m(5) + m(6);
t155 = m(4) + t188;
t146 = t155 * pkin(2) + mrSges(3,1);
t190 = -t110 * mrSges(3,2) + t146 * t111;
t45 = t117 * t50;
t24 = t116 * t42 + t45;
t19 = pkin(8) * t109 + t24;
t13 = qJD(4) * t122 - t118 * t19;
t14 = qJD(4) * t118 + t122 * t19;
t141 = -t118 * t13 + t122 * t14;
t112 = qJ(3) + t115;
t101 = pkin(9) + t112;
t94 = cos(t101);
t91 = t94 * pkin(8);
t93 = sin(t101);
t189 = t93 * mrSges(5,1) - m(6) * (-pkin(4) * t93 + t91);
t185 = mrSges(6,3) * t13;
t184 = pkin(1) * t124;
t183 = pkin(3) * t116;
t182 = pkin(3) * t117;
t113 = qJDD(1) + qJDD(2);
t108 = qJDD(3) + t113;
t77 = -qJD(2) * t154 + qJDD(1) * t184;
t62 = pkin(2) * t113 + t77;
t160 = qJD(2) * t124;
t78 = (qJD(1) * t160 + qJDD(1) * t120) * pkin(1);
t22 = -t50 * qJD(3) - t119 * t78 + t123 * t62;
t16 = pkin(3) * t108 + t22;
t21 = t49 * qJD(3) + t119 * t62 + t123 * t78;
t9 = t116 * t16 + t117 * t21;
t6 = pkin(8) * t108 + t9;
t2 = t13 * qJD(5) + qJDD(4) * t118 + t122 * t6;
t181 = t122 * t2;
t105 = pkin(2) + t184;
t72 = -pkin(1) * t162 + t123 * t105;
t67 = pkin(3) + t72;
t73 = pkin(1) * t161 + t105 * t119;
t35 = t116 * t67 + t117 * t73;
t104 = pkin(2) * t123 + pkin(3);
t69 = pkin(2) * t163 + t116 * t104;
t177 = mrSges(6,2) * t118;
t176 = Ifges(6,4) * t118;
t175 = Ifges(6,4) * t122;
t174 = Ifges(6,2) * t122;
t169 = t122 * mrSges(6,1);
t167 = qJD(5) * t14;
t166 = t109 * t118;
t165 = t109 * t122;
t159 = qJD(3) * t119;
t158 = qJD(3) * t123;
t157 = qJD(5) * t118;
t156 = qJD(5) * t122;
t61 = t108 * t118 + t109 * t156;
t48 = qJDD(5) * mrSges(6,1) - mrSges(6,3) * t61;
t80 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t165;
t151 = -qJD(5) * t80 - t48;
t149 = t188 * pkin(3) + mrSges(4,1);
t3 = qJDD(4) * t122 - t118 * t6 - t167;
t147 = (-t3 - t167) * mrSges(6,3);
t85 = -t169 + t177;
t144 = t174 + t176;
t8 = -t116 * t21 + t117 * t16;
t34 = -t116 * t73 + t117 * t67;
t142 = t118 * t14 + t122 * t13;
t138 = m(6) * pkin(4) + mrSges(5,1) + t169;
t137 = mrSges(2,1) + (m(3) + t155) * pkin(1);
t136 = (g(2) * t94 + g(3) * t93) * mrSges(6,1);
t68 = -pkin(2) * t164 + t104 * t117;
t134 = t118 * (Ifges(6,1) * t122 - t176);
t102 = sin(t112);
t103 = cos(t112);
t132 = -t102 * mrSges(4,2) - t93 * mrSges(5,2) + t149 * t103 - t94 * t177;
t131 = t103 * mrSges(4,2) + t149 * t102 - t93 * t177 + (mrSges(5,2) - mrSges(6,3)) * t94;
t130 = m(6) * (-t142 * qJD(5) - t118 * t3 + t181);
t5 = -pkin(4) * t108 - t8;
t51 = Ifges(6,6) * qJD(5) + t144 * t109;
t87 = Ifges(6,4) * t165;
t52 = Ifges(6,1) * t166 + Ifges(6,5) * qJD(5) + t87;
t60 = t108 * t122 - t109 * t157;
t129 = t22 * mrSges(4,1) + t8 * mrSges(5,1) - t21 * mrSges(4,2) + mrSges(6,3) * t181 + t5 * t85 + (Ifges(6,1) * t61 + Ifges(6,4) * t60) * t187 + t122 * (Ifges(6,4) * t61 + Ifges(6,2) * t60) / 0.2e1 + t60 * t144 / 0.2e1 + t61 * (Ifges(6,1) * t118 + t175) / 0.2e1 - t51 * t157 / 0.2e1 + (t52 + t109 * (-Ifges(6,2) * t118 + t175)) * t156 / 0.2e1 + (Ifges(5,3) + Ifges(4,3)) * t108 + (0.2e1 * Ifges(6,5) * t187 + Ifges(6,6) * t122) * qJDD(5) + (t134 * t109 / 0.2e1 + t191) * qJD(5);
t128 = -m(6) * (-t94 * pkin(4) - t93 * pkin(8)) + t93 * mrSges(6,3) + t94 * mrSges(5,1) + t132;
t127 = t77 * mrSges(3,1) - t9 * mrSges(5,2) + Ifges(3,3) * t113 + t129;
t126 = mrSges(3,2) * t111 + t146 * t110 + t131;
t125 = cos(qJ(1));
t121 = sin(qJ(1));
t100 = -pkin(4) - t182;
t99 = pkin(8) + t183;
t79 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t166;
t64 = pkin(8) + t69;
t63 = -pkin(4) - t68;
t59 = t85 * t109;
t47 = -qJDD(5) * mrSges(6,2) + mrSges(6,3) * t60;
t38 = -t105 * t159 + (t140 * qJD(2) - t120 * t158) * pkin(1);
t37 = t105 * t158 + (t139 * qJD(2) - t120 * t159) * pkin(1);
t33 = -mrSges(6,1) * t60 + mrSges(6,2) * t61;
t30 = pkin(8) + t35;
t29 = -pkin(4) - t34;
t26 = t117 * t49 - t171;
t25 = t116 * t49 + t45;
t12 = t116 * t38 + t117 * t37;
t11 = t116 * t37 - t117 * t38;
t1 = [t127 + m(4) * (t21 * t73 + t22 * t72 + t37 * t50 + t38 * t49) + m(5) * (-t11 * t23 + t12 * t24 + t34 * t8 + t35 * t9) + m(6) * (t11 * t18 + t141 * t12 + t29 * t5) - t78 * mrSges(3,2) + t11 * t59 + (-t121 * mrSges(2,2) + (m(6) * pkin(8) + mrSges(6,3)) * t93 + t138 * t94 + t137 * t125 + t132 + t190) * g(2) + t29 * t33 + (mrSges(4,1) * t38 - mrSges(5,1) * t11 - mrSges(4,2) * t37 - mrSges(5,2) * t12) * t109 + (mrSges(4,1) * t72 + mrSges(5,1) * t34 - mrSges(4,2) * t73 - mrSges(5,2) * t35) * t108 + (-m(6) * t91 + mrSges(2,2) * t125 + t137 * t121 + t138 * t93 + t126) * g(3) + Ifges(2,3) * qJDD(1) + t30 * t130 + (t12 * t80 + t30 * t47 + (-t30 * t79 - t185) * qJD(5)) * t122 + (-t12 * t79 + t151 * t30 + t147) * t118 + (m(3) * (t120 * t78 + t124 * t77) + (-t113 * t120 - t114 * t160) * mrSges(3,2) + (-qJD(2) * t114 * t120 + t113 * t124) * mrSges(3,1)) * pkin(1); t127 + t64 * t130 + (t68 * mrSges(5,1) - t69 * mrSges(5,2) + (mrSges(4,1) * t123 - mrSges(4,2) * t119) * pkin(2)) * t108 + (t64 * t47 - t178 * t80 + (-t64 * t79 - t185) * qJD(5) + t136) * t122 + (t114 * t153 - t78) * mrSges(3,2) - t179 * t59 + (t128 + t190) * g(2) + t63 * t33 + t114 * mrSges(3,1) * t154 + (t178 * mrSges(5,2) + (-pkin(2) * t158 + t71) * mrSges(4,2) + t179 * mrSges(5,1) + (-pkin(2) * t159 - t70) * mrSges(4,1)) * t109 + (t126 + t189) * g(3) + (t151 * t64 + t178 * t79 + t147) * t118 + (-t178 * t141 - t179 * t18 + t5 * t63) * m(6) + (-t178 * t24 + t179 * t23 + t68 * t8 + t69 * t9) * m(5) + ((t119 * t21 + t123 * t22 + (-t119 * t49 + t123 * t50) * qJD(3)) * pkin(2) - t49 * t70 - t50 * t71) * m(4); t129 + (mrSges(4,1) * t50 + mrSges(5,1) * t25 + mrSges(4,2) * t49 + mrSges(5,2) * t26) * t109 + t100 * t33 - t25 * t59 + t108 * mrSges(5,1) * t182 + t128 * g(2) + t99 * t130 + (-t108 * t183 - t9) * mrSges(5,2) + (t131 + t189) * g(3) + (-t26 * t80 + t99 * t47 + (-t79 * t99 - t185) * qJD(5) + t136) * t122 + (t151 * t99 + t26 * t79 + t147) * t118 + (t100 * t5 - t141 * t26 - t18 * t25) * m(6) + (t23 * t25 - t24 * t26 + (t116 * t9 + t117 * t8) * pkin(3)) * m(5); t118 * t47 + t122 * t48 + (-t118 * t79 + t122 * t80) * qJD(5) + (t141 * qJD(5) + t118 * t2 + t122 * t3 - g(1)) * m(6) + (qJDD(4) - g(1)) * m(5); t3 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,5) * t61 + Ifges(6,6) * t60 + Ifges(6,3) * qJDD(5) + g(1) * t85 - t13 * t80 + t14 * t79 + (t51 * t187 + (-t134 / 0.2e1 + t174 * t187) * t109 + t142 * mrSges(6,3) - (t52 + t87) * t122 / 0.2e1 - t191) * t109 + (-g(2) * t93 + g(3) * t94) * t145;];
tau = t1;
