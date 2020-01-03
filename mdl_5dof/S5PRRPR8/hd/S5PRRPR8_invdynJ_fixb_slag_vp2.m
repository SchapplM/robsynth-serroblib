% Calculate vector of inverse dynamics joint torques for
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRRPR8_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR8_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR8_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR8_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR8_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR8_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR8_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR8_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR8_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:42:19
% EndTime: 2019-12-31 17:42:24
% DurationCPUTime: 2.20s
% Computational Cost: add. (2098->254), mult. (3771->354), div. (0->0), fcn. (2545->14), ass. (0->130)
t114 = sin(qJ(5));
t186 = mrSges(6,2) * t114;
t109 = qJ(2) + qJ(3);
t103 = pkin(9) + t109;
t94 = sin(t103);
t95 = cos(t103);
t214 = t94 * t186 + t95 * (m(6) * pkin(7) + mrSges(6,3));
t197 = t114 / 0.2e1;
t104 = sin(t109);
t105 = cos(t109);
t213 = mrSges(4,1) * t104 + mrSges(5,1) * t94 + mrSges(4,2) * t105 + mrSges(5,2) * t95;
t117 = cos(qJ(5));
t170 = t117 * mrSges(6,1);
t212 = t170 - t186;
t108 = qJD(2) + qJD(3);
t110 = sin(pkin(9));
t115 = sin(qJ(3));
t118 = cos(qJ(3));
t116 = sin(qJ(2));
t159 = qJD(1) * t116;
t119 = cos(qJ(2));
t86 = qJD(2) * pkin(2) + qJD(1) * t119;
t49 = -t115 * t159 + t118 * t86;
t42 = pkin(3) * t108 + t49;
t112 = cos(pkin(9));
t50 = t115 * t86 + t118 * t159;
t43 = t112 * t50;
t24 = t110 * t42 + t43;
t19 = pkin(7) * t108 + t24;
t13 = qJD(4) * t117 - t114 * t19;
t107 = qJDD(2) + qJDD(3);
t156 = qJD(1) * qJD(2);
t147 = t116 * t156;
t76 = t119 * qJDD(1) - t147;
t68 = qJDD(2) * pkin(2) + t76;
t146 = t119 * t156;
t77 = qJDD(1) * t116 + t146;
t22 = -t50 * qJD(3) - t115 * t77 + t118 * t68;
t16 = pkin(3) * t107 + t22;
t21 = t49 * qJD(3) + t115 * t68 + t118 * t77;
t9 = t110 * t16 + t112 * t21;
t6 = pkin(7) * t107 + t9;
t2 = qJD(5) * t13 + qJDD(4) * t114 + t117 * t6;
t14 = qJD(4) * t114 + t117 * t19;
t3 = -qJD(5) * t14 + qJDD(4) * t117 - t114 * t6;
t211 = -t114 * t3 + t117 * t2;
t111 = sin(pkin(8));
t113 = cos(pkin(8));
t210 = g(1) * t113 + g(2) * t111;
t136 = mrSges(6,1) * t114 + mrSges(6,2) * t117;
t175 = t110 * t50;
t23 = t112 * t42 - t175;
t18 = -pkin(4) * t108 - t23;
t208 = t18 * t136 + qJD(5) * (Ifges(6,5) * t117 - Ifges(6,6) * t114) / 0.2e1;
t167 = t108 * t114;
t74 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t167;
t166 = t108 * t117;
t75 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t166;
t131 = -t114 * t74 + t117 * t75;
t158 = qJD(5) * t114;
t63 = t107 * t117 - t108 * t158;
t45 = -qJDD(5) * mrSges(6,2) + mrSges(6,3) * t63;
t157 = qJD(5) * t117;
t64 = t107 * t114 + t108 * t157;
t46 = qJDD(5) * mrSges(6,1) - mrSges(6,3) * t64;
t206 = -t114 * t46 + t117 * t45;
t31 = -mrSges(6,1) * t63 + mrSges(6,2) * t64;
t8 = -t110 * t21 + t112 * t16;
t5 = -pkin(4) * t107 - t8;
t204 = m(6) * t5 + t31;
t202 = -t105 * mrSges(4,1) + mrSges(4,2) * t104 + (-mrSges(5,1) - t212) * t95 + (mrSges(5,2) - mrSges(6,3)) * t94;
t61 = t212 * t108;
t201 = m(5) * t23 - m(6) * t18 + t108 * mrSges(5,1) + t61;
t133 = t114 * t14 + t117 * t13;
t122 = -qJD(5) * t133 + t211;
t200 = m(6) * t122 - t74 * t157 - t75 * t158 + t206;
t132 = -t114 * t13 + t117 * t14;
t199 = m(5) * t24 + m(6) * t132 - t108 * mrSges(5,2) + t131;
t198 = pkin(4) * t94;
t195 = pkin(2) * t116;
t194 = pkin(2) * t118;
t193 = pkin(3) * t104;
t98 = pkin(3) * t105;
t192 = pkin(3) * t110;
t191 = pkin(3) * t112;
t106 = t119 * pkin(2);
t162 = t112 * t115;
t99 = pkin(3) + t194;
t60 = pkin(2) * t162 + t110 * t99;
t185 = Ifges(6,4) * t114;
t184 = Ifges(6,4) * t117;
t183 = Ifges(6,2) * t117;
t182 = pkin(2) * qJD(3);
t181 = t107 * mrSges(5,1);
t180 = t107 * mrSges(5,2);
t179 = t108 * mrSges(4,1);
t177 = t108 * mrSges(4,2);
t165 = t110 * t115;
t164 = t111 * t114;
t163 = t111 * t117;
t161 = t113 * t114;
t160 = t113 * t117;
t155 = t95 * pkin(4) + t94 * pkin(7) + t98;
t153 = t94 * t170;
t152 = t108 * t182;
t143 = t214 * t111;
t142 = t214 * t113;
t140 = -g(1) * t111 + g(2) * t113;
t135 = t183 + t185;
t73 = t115 * t119 + t116 * t118;
t72 = -t115 * t116 + t118 * t119;
t59 = -pkin(2) * t165 + t112 * t99;
t128 = t114 * (Ifges(6,1) * t117 - t185);
t79 = -t193 - t195;
t124 = m(6) * (t79 - t198) - t153;
t123 = m(6) * (-t193 - t198) - t153;
t51 = Ifges(6,6) * qJD(5) + t108 * t135;
t87 = Ifges(6,4) * t166;
t52 = Ifges(6,1) * t167 + Ifges(6,5) * qJD(5) + t87;
t121 = -t21 * mrSges(4,2) - t9 * mrSges(5,2) + t22 * mrSges(4,1) + (Ifges(6,1) * t64 + Ifges(6,4) * t63) * t197 + t117 * (Ifges(6,4) * t64 + Ifges(6,2) * t63) / 0.2e1 + t63 * t135 / 0.2e1 + t64 * (Ifges(6,1) * t114 + t184) / 0.2e1 - t5 * t212 - t51 * t158 / 0.2e1 + t8 * mrSges(5,1) + (t52 + t108 * (-Ifges(6,2) * t114 + t184)) * t157 / 0.2e1 + (Ifges(5,3) + Ifges(4,3)) * t107 + (-t13 * t157 - t14 * t158 + t211) * mrSges(6,3) + (0.2e1 * Ifges(6,5) * t197 + Ifges(6,6) * t117) * qJDD(5) + (t128 * t108 / 0.2e1 + t208) * qJD(5);
t120 = qJD(2) ^ 2;
t67 = t72 * qJD(1);
t66 = t73 * qJD(1);
t36 = t108 * t73;
t35 = t108 * t72;
t34 = t110 * t72 + t112 * t73;
t33 = t110 * t73 - t112 * t72;
t12 = -t110 * t36 + t112 * t35;
t11 = t110 * t35 + t112 * t36;
t1 = [m(2) * qJDD(1) - t11 * t61 + t33 * t31 + t131 * t12 + (-qJDD(2) * t116 - t119 * t120) * mrSges(3,2) + (qJDD(2) * t119 - t116 * t120) * mrSges(3,1) + (mrSges(4,1) * t72 - mrSges(5,1) * t33 - mrSges(4,2) * t73) * t107 + (-mrSges(4,1) * t36 - mrSges(5,1) * t11 - mrSges(4,2) * t35 - mrSges(5,2) * t12) * t108 + (-t180 + (-t114 * t75 - t117 * t74) * qJD(5) + t206) * t34 + (-m(2) - m(3) - m(4) - m(5) - m(6)) * g(3) + m(6) * (t11 * t18 + t12 * t132 + t122 * t34 + t33 * t5) + m(3) * (t116 * t77 + t119 * t76) + m(4) * (t21 * t73 + t22 * t72 + t35 * t50 - t36 * t49) + m(5) * (-t11 * t23 + t12 * t24 - t33 * t8 + t34 * t9); m(5) * (t59 * t8 + t60 * t9) - g(1) * (t113 * t124 + t142) - g(2) * (t111 * t124 + t143) + t66 * t179 - t60 * t180 + t121 + Ifges(3,3) * qJDD(2) + t67 * t177 + t59 * t181 + t204 * (-pkin(4) - t59) + (-t115 * pkin(2) * t107 - t118 * t152) * mrSges(4,2) + (-t77 + t146) * mrSges(3,2) + (t107 * t194 - t115 * t152) * mrSges(4,1) + (t76 + t147) * mrSges(3,1) + ((t115 * t21 + t118 * t22 + (-t115 * t49 + t118 * t50) * qJD(3)) * pkin(2) + t49 * t66 - t50 * t67) * m(4) + t201 * (t110 * t67 + t112 * t66) - t199 * (-t110 * t66 + t112 * t67) + t200 * (pkin(7) + t60) + (-m(6) * (t106 + t155) - mrSges(3,1) * t119 + mrSges(3,2) * t116 - m(5) * (t98 + t106) - m(4) * t106 + t202) * g(3) + t210 * (m(4) * t195 - m(5) * t79 + mrSges(3,1) * t116 + mrSges(3,2) * t119 + t213) + (-t201 * (t110 * t118 + t162) + t199 * (t112 * t118 - t165)) * t182; m(5) * (t110 * t9 + t112 * t8) * pkin(3) - t180 * t192 - g(1) * (t113 * t123 + t142) - g(2) * (t111 * t123 + t143) + t121 + t49 * t177 + t50 * t179 + t181 * t191 + t204 * (-pkin(4) - t191) + t201 * (t110 * t49 + t43) - t199 * (t112 * t49 - t175) + (-m(5) * t98 - m(6) * t155 + t202) * g(3) + t200 * (pkin(7) + t192) + (m(5) * t193 + t213) * t210; t114 * t45 + t117 * t46 + t131 * qJD(5) + (qJD(5) * t132 + t114 * t2 + t117 * t3 + t140) * m(6) + (qJDD(4) + t140) * m(5); Ifges(6,5) * t64 + Ifges(6,6) * t63 + Ifges(6,3) * qJDD(5) - t2 * mrSges(6,2) + t3 * mrSges(6,1) - t13 * t75 + t14 * t74 - g(1) * ((-t161 * t95 + t163) * mrSges(6,1) + (-t160 * t95 - t164) * mrSges(6,2)) - g(2) * ((-t164 * t95 - t160) * mrSges(6,1) + (-t163 * t95 + t161) * mrSges(6,2)) + g(3) * t136 * t94 + (t51 * t197 + (-t128 / 0.2e1 + t183 * t197) * t108 + t133 * mrSges(6,3) - (t52 + t87) * t117 / 0.2e1 - t208) * t108;];
tau = t1;
