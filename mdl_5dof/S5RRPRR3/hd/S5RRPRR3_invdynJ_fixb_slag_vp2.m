% Calculate vector of inverse dynamics joint torques for
% S5RRPRR3
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 10:34
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRR3_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR3_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR3_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR3_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR3_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR3_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR3_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR3_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR3_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:33:59
% EndTime: 2022-01-20 10:34:04
% DurationCPUTime: 1.88s
% Computational Cost: add. (2960->271), mult. (4825->358), div. (0->0), fcn. (2646->16), ass. (0->135)
t118 = sin(qJ(5));
t184 = t118 / 0.2e1;
t119 = sin(qJ(4));
t123 = cos(qJ(4));
t116 = sin(pkin(9));
t124 = cos(qJ(2));
t117 = cos(pkin(9));
t120 = sin(qJ(2));
t155 = t117 * t120;
t135 = pkin(1) * (-t116 * t124 - t155);
t62 = qJD(1) * t135;
t156 = t116 * t120;
t134 = pkin(1) * (t117 * t124 - t156);
t64 = qJD(1) * t134;
t179 = pkin(2) * t117;
t100 = pkin(3) + t179;
t180 = pkin(2) * t116;
t67 = t100 * t123 - t119 * t180;
t174 = -t67 * qJD(4) + t119 * t62 + t123 * t64;
t122 = cos(qJ(5));
t170 = mrSges(6,1) * t122;
t195 = -mrSges(5,1) - t170;
t114 = qJD(1) + qJD(2);
t165 = pkin(1) * qJD(1);
t147 = t120 * t165;
t80 = pkin(2) * t114 + t124 * t165;
t42 = -t116 * t147 + t117 * t80;
t37 = pkin(3) * t114 + t42;
t43 = t116 * t80 + t117 * t147;
t17 = -t119 * t43 + t123 * t37;
t115 = qJ(1) + qJ(2);
t109 = pkin(9) + t115;
t102 = qJ(4) + t109;
t93 = sin(t102);
t94 = cos(t102);
t194 = g(1) * t94 + g(2) * t93;
t142 = mrSges(6,1) * t118 + mrSges(6,2) * t122;
t108 = qJD(4) + t114;
t15 = -pkin(4) * t108 - t17;
t193 = t15 * t142 + qJD(5) * (Ifges(6,5) * t122 - Ifges(6,6) * t118) / 0.2e1;
t69 = t100 * t119 + t123 * t180;
t173 = qJD(4) * t69 - t119 * t64 + t123 * t62;
t169 = mrSges(6,2) * t118;
t192 = t93 * mrSges(5,2) + t94 * t169;
t18 = t119 * t37 + t123 * t43;
t16 = pkin(8) * t108 + t18;
t13 = qJD(3) * t122 - t118 * t16;
t14 = qJD(3) * t118 + t122 * t16;
t138 = -t118 * t13 + t122 * t14;
t189 = m(6) * pkin(4);
t89 = t94 * pkin(8);
t191 = (t189 - t195) * t93 + t94 * mrSges(5,2) - m(6) * t89;
t113 = qJDD(1) + qJDD(2);
t181 = pkin(1) * t124;
t73 = -qJD(2) * t147 + qJDD(1) * t181;
t58 = pkin(2) * t113 + t73;
t153 = qJD(2) * t124;
t74 = (qJD(1) * t153 + qJDD(1) * t120) * pkin(1);
t32 = -t116 * t74 + t117 * t58;
t22 = pkin(3) * t113 + t32;
t33 = t116 * t58 + t117 * t74;
t9 = -qJD(4) * t18 - t119 * t33 + t123 * t22;
t190 = m(3) * pkin(1);
t188 = m(5) + m(6);
t182 = mrSges(6,3) * t13;
t111 = cos(t115);
t101 = pkin(2) * t111;
t107 = qJDD(4) + t113;
t8 = qJD(4) * t17 + t119 * t22 + t123 * t33;
t5 = pkin(8) * t107 + t8;
t2 = qJD(5) * t13 + qJDD(3) * t118 + t122 * t5;
t178 = t122 * t2;
t159 = qJD(5) * t14;
t3 = qJDD(3) * t122 - t118 * t5 - t159;
t177 = t3 * t118;
t103 = pkin(2) + t181;
t66 = -pkin(1) * t156 + t103 * t117;
t59 = pkin(3) + t66;
t68 = pkin(1) * t155 + t103 * t116;
t29 = t119 * t59 + t123 * t68;
t172 = -mrSges(6,3) * t94 - t169 * t93;
t171 = pkin(4) * t94 + pkin(8) * t93;
t168 = Ifges(6,4) * t118;
t167 = Ifges(6,4) * t122;
t166 = Ifges(6,2) * t122;
t99 = cos(t109);
t92 = pkin(3) * t99;
t160 = t92 + t101;
t158 = t108 * t118;
t157 = t108 * t122;
t125 = cos(qJ(1));
t112 = t125 * pkin(1);
t154 = t101 + t112;
t152 = qJD(5) * t118;
t151 = qJD(5) * t122;
t150 = m(4) + t188;
t149 = -t93 * mrSges(6,3) + t195 * t94;
t57 = t107 * t118 + t108 * t151;
t41 = qJDD(5) * mrSges(6,1) - mrSges(6,3) * t57;
t76 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t157;
t145 = -qJD(5) * t76 - t41;
t144 = t160 + t171;
t143 = (-t3 - t159) * mrSges(6,3);
t81 = t169 - t170;
t141 = t166 + t168;
t139 = t118 * t14 + t122 * t13;
t28 = -t119 * t68 + t123 * t59;
t133 = t118 * (Ifges(6,1) * t122 - t168);
t110 = sin(t115);
t98 = sin(t109);
t131 = -mrSges(3,1) * t111 - mrSges(4,1) * t99 + t110 * mrSges(3,2) + t98 * mrSges(4,2) + t149;
t130 = -qJD(5) * t139 - t177;
t45 = Ifges(6,6) * qJD(5) + t108 * t141;
t84 = Ifges(6,4) * t157;
t46 = Ifges(6,1) * t158 + Ifges(6,5) * qJD(5) + t84;
t56 = t107 * t122 - t108 * t152;
t6 = -pkin(4) * t107 - t9;
t129 = t9 * mrSges(5,1) + mrSges(6,3) * t178 + t6 * t81 + (Ifges(6,1) * t57 + Ifges(6,4) * t56) * t184 + t122 * (Ifges(6,4) * t57 + Ifges(6,2) * t56) / 0.2e1 + t56 * t141 / 0.2e1 + t57 * (Ifges(6,1) * t118 + t167) / 0.2e1 - t45 * t152 / 0.2e1 + Ifges(5,3) * t107 + (t46 + t108 * (-Ifges(6,2) * t118 + t167)) * t151 / 0.2e1 + (0.2e1 * Ifges(6,5) * t184 + Ifges(6,6) * t122) * qJDD(5) + (t133 * t108 / 0.2e1 + t193) * qJD(5);
t128 = m(6) * (t130 + t178);
t127 = t73 * mrSges(3,1) + t32 * mrSges(4,1) - t74 * mrSges(3,2) + t129 + (Ifges(3,3) + Ifges(4,3)) * t113;
t126 = mrSges(3,2) * t111 + t99 * mrSges(4,2) + (pkin(2) * t150 + mrSges(3,1)) * t110 + (pkin(3) * t188 + mrSges(4,1)) * t98 + t172;
t121 = sin(qJ(1));
t75 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t158;
t65 = qJD(2) * t134;
t63 = qJD(2) * t135;
t61 = pkin(8) + t69;
t60 = -pkin(4) - t67;
t55 = t81 * t108;
t40 = -qJDD(5) * mrSges(6,2) + mrSges(6,3) * t56;
t30 = -mrSges(6,1) * t56 + mrSges(6,2) * t57;
t24 = pkin(8) + t29;
t23 = -pkin(4) - t28;
t11 = qJD(4) * t29 + t119 * t65 - t123 * t63;
t10 = qJD(4) * t28 + t119 * t63 + t123 * t65;
t1 = [t127 + (mrSges(2,2) * t125 + (mrSges(2,1) + (m(3) + t150) * pkin(1)) * t121 + t126 + t191) * g(1) + (-t10 * t108 - t107 * t29 - t8) * mrSges(5,2) + (-t113 * t68 - t114 * t65 - t33) * mrSges(4,2) + (t107 * t28 - t108 * t11) * mrSges(5,1) + (t113 * t66 + t114 * t63) * mrSges(4,1) + t24 * t128 + m(6) * (t138 * t10 + t11 * t15 + t23 * t6) + (t121 * mrSges(2,2) - m(6) * (t112 + t144) - m(5) * (t92 + t154) - m(4) * t154 + (-mrSges(2,1) - t190) * t125 + t131 + t192) * g(2) + (t120 * t74 + t124 * t73) * t190 + t11 * t55 + t23 * t30 + (-t10 * t75 + t145 * t24 + t143) * t118 + ((-t113 * t120 - t114 * t153) * mrSges(3,2) + (-qJD(2) * t114 * t120 + t113 * t124) * mrSges(3,1)) * pkin(1) + (t10 * t76 + t24 * t40 + (-t24 * t75 - t182) * qJD(5)) * t122 + m(4) * (t32 * t66 + t33 * t68 + t42 * t63 + t43 * t65) + m(5) * (t10 * t18 - t11 * t17 + t28 * t9 + t29 * t8) + Ifges(2,3) * qJDD(1); t127 + (t93 * mrSges(5,1) + t126) * g(1) + t173 * t55 + (-t113 * t180 - t33) * mrSges(4,2) + t113 * mrSges(4,1) * t179 + (g(2) * mrSges(6,2) * t94 + t145 * t61 + t174 * t75 + t143) * t118 + (g(1) * mrSges(6,1) * t93 + t61 * t40 - t174 * t76 + (-t61 * t75 - t182) * qJD(5)) * t122 + (-t69 * t107 + t108 * t174 + t194 - t8) * mrSges(5,2) + t60 * t30 + (-t62 * mrSges(4,1) + t64 * mrSges(4,2) + (mrSges(3,1) * t120 + mrSges(3,2) * t124) * t165) * t114 + (t107 * t67 - t108 * t173) * mrSges(5,1) + t131 * g(2) + t61 * t128 + ((pkin(4) * t93 - t89) * g(1) + t6 * t60 - t144 * g(2) + t173 * t15 - t174 * t138) * m(6) + (-t160 * g(2) - t17 * t173 - t174 * t18 + t67 * t9 + t69 * t8) * m(5) + ((t116 * t33 + t117 * t32) * pkin(2) - t42 * t62 - t43 * t64 - t101 * g(2)) * m(4); m(6) * (qJD(5) * t138 + t118 * t2 + t122 * t3) + t76 * t151 + t118 * t40 - t75 * t152 + t122 * t41 + (m(4) + m(5)) * qJDD(3) - t150 * g(3); t129 + t130 * mrSges(6,3) + (t172 + t191) * g(1) + (t108 * mrSges(5,1) - t55) * t18 + (mrSges(5,2) * t108 + t118 * t75 - t122 * t76) * t17 - t6 * t189 + (-m(6) * t171 + t149 + t192) * g(2) + (m(6) * (-t13 * t151 - t14 * t152 - t177 + t178) + t122 * t40 - t118 * t41 - t76 * t152 - t75 * t151) * pkin(8) - pkin(4) * t30 - t8 * mrSges(5,2) - m(6) * (t138 * t17 + t15 * t18); t3 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,5) * t57 + Ifges(6,6) * t56 + Ifges(6,3) * qJDD(5) + g(3) * t81 - t13 * t76 + t14 * t75 + (t45 * t184 + (-t133 / 0.2e1 + t166 * t184) * t108 + t139 * mrSges(6,3) - (t46 + t84) * t122 / 0.2e1 - t193) * t108 + t194 * t142;];
tau = t1;
