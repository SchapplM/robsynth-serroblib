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
% Datum: 2019-12-05 18:31
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 18:30:18
% EndTime: 2019-12-05 18:30:22
% DurationCPUTime: 1.80s
% Computational Cost: add. (2960->263), mult. (4825->351), div. (0->0), fcn. (2646->16), ass. (0->131)
t112 = sin(qJ(5));
t177 = t112 / 0.2e1;
t113 = sin(qJ(4));
t117 = cos(qJ(4));
t110 = sin(pkin(9));
t118 = cos(qJ(2));
t111 = cos(pkin(9));
t114 = sin(qJ(2));
t150 = t111 * t114;
t130 = pkin(1) * (-t110 * t118 - t150);
t62 = qJD(1) * t130;
t151 = t110 * t114;
t129 = pkin(1) * (t111 * t118 - t151);
t64 = qJD(1) * t129;
t173 = pkin(2) * t110;
t172 = pkin(2) * t111;
t96 = pkin(3) + t172;
t67 = -t113 * t173 + t117 * t96;
t168 = -t67 * qJD(4) + t113 * t62 + t117 * t64;
t108 = qJD(1) + qJD(2);
t159 = qJD(1) * pkin(1);
t145 = t114 * t159;
t80 = pkin(2) * t108 + t118 * t159;
t42 = -t110 * t145 + t111 * t80;
t37 = pkin(3) * t108 + t42;
t43 = t110 * t80 + t111 * t145;
t17 = -t113 * t43 + t117 * t37;
t116 = cos(qJ(5));
t139 = mrSges(6,1) * t112 + mrSges(6,2) * t116;
t103 = qJD(4) + t108;
t15 = -pkin(4) * t103 - t17;
t184 = t15 * t139 + qJD(5) * (Ifges(6,5) * t116 - Ifges(6,6) * t112) / 0.2e1;
t69 = t113 * t96 + t117 * t173;
t167 = t69 * qJD(4) - t113 * t64 + t117 * t62;
t18 = t113 * t37 + t117 * t43;
t16 = pkin(8) * t103 + t18;
t13 = qJD(3) * t116 - t112 * t16;
t14 = qJD(3) * t112 + t116 * t16;
t135 = -t112 * t13 + t116 * t14;
t164 = mrSges(6,1) * t116;
t182 = m(6) * pkin(4);
t133 = mrSges(5,1) + t164 + t182;
t109 = qJ(1) + qJ(2);
t104 = pkin(9) + t109;
t97 = qJ(4) + t104;
t90 = cos(t97);
t87 = t90 * pkin(8);
t89 = sin(t97);
t183 = -m(6) * t87 + t90 * mrSges(5,2) + t133 * t89;
t107 = qJDD(1) + qJDD(2);
t174 = pkin(1) * t118;
t73 = -qJD(2) * t145 + qJDD(1) * t174;
t58 = pkin(2) * t107 + t73;
t149 = qJD(2) * t118;
t74 = (qJD(1) * t149 + qJDD(1) * t114) * pkin(1);
t32 = -t110 * t74 + t111 * t58;
t22 = pkin(3) * t107 + t32;
t33 = t110 * t58 + t111 * t74;
t9 = -t18 * qJD(4) - t113 * t33 + t117 * t22;
t181 = m(5) + m(6);
t179 = g(2) * t89;
t178 = g(3) * t90;
t175 = mrSges(6,3) * t13;
t102 = qJDD(4) + t107;
t8 = t17 * qJD(4) + t113 * t22 + t117 * t33;
t5 = pkin(8) * t102 + t8;
t2 = t13 * qJD(5) + qJDD(3) * t112 + t116 * t5;
t171 = t116 * t2;
t154 = qJD(5) * t14;
t3 = qJDD(3) * t116 - t112 * t5 - t154;
t170 = t3 * t112;
t98 = pkin(2) + t174;
t66 = -pkin(1) * t151 + t111 * t98;
t59 = pkin(3) + t66;
t68 = pkin(1) * t150 + t110 * t98;
t29 = t113 * t59 + t117 * t68;
t163 = mrSges(6,2) * t112;
t166 = -t90 * mrSges(6,3) - t89 * t163;
t165 = -t89 * mrSges(5,2) - t90 * t163;
t162 = Ifges(6,4) * t112;
t161 = Ifges(6,4) * t116;
t160 = Ifges(6,2) * t116;
t153 = t103 * t112;
t152 = t103 * t116;
t148 = qJD(5) * t112;
t147 = qJD(5) * t116;
t146 = m(4) + t181;
t57 = t102 * t112 + t103 * t147;
t41 = qJDD(5) * mrSges(6,1) - mrSges(6,3) * t57;
t76 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t152;
t143 = -qJD(5) * t76 - t41;
t142 = t181 * pkin(3) + mrSges(4,1);
t141 = (-t3 - t154) * mrSges(6,3);
t140 = t146 * pkin(2) + mrSges(3,1);
t81 = t163 - t164;
t138 = t160 + t162;
t136 = t112 * t14 + t116 * t13;
t28 = -t113 * t68 + t117 * t59;
t132 = mrSges(2,1) + (m(3) + t146) * pkin(1);
t128 = t112 * (Ifges(6,1) * t116 - t162);
t126 = t133 * t90;
t125 = -t136 * qJD(5) - t170;
t45 = Ifges(6,6) * qJD(5) + t138 * t103;
t84 = Ifges(6,4) * t152;
t46 = Ifges(6,1) * t153 + Ifges(6,5) * qJD(5) + t84;
t56 = t102 * t116 - t103 * t148;
t6 = -pkin(4) * t102 - t9;
t124 = t9 * mrSges(5,1) + mrSges(6,3) * t171 + t6 * t81 + (Ifges(6,1) * t57 + Ifges(6,4) * t56) * t177 + t116 * (Ifges(6,4) * t57 + Ifges(6,2) * t56) / 0.2e1 + t56 * t138 / 0.2e1 + t57 * (Ifges(6,1) * t112 + t161) / 0.2e1 - t45 * t148 / 0.2e1 + Ifges(5,3) * t102 + (t46 + t103 * (-Ifges(6,2) * t112 + t161)) * t147 / 0.2e1 + (0.2e1 * Ifges(6,5) * t177 + Ifges(6,6) * t116) * qJDD(5) + (t128 * t103 / 0.2e1 + t184) * qJD(5);
t123 = m(6) * (t125 + t171);
t105 = sin(t109);
t106 = cos(t109);
t94 = sin(t104);
t95 = cos(t104);
t122 = -t105 * mrSges(3,2) - t94 * mrSges(4,2) + t140 * t106 + t142 * t95 + t165;
t121 = t73 * mrSges(3,1) + t32 * mrSges(4,1) - t74 * mrSges(3,2) + t124 + (Ifges(3,3) + Ifges(4,3)) * t107;
t120 = mrSges(3,2) * t106 + t95 * mrSges(4,2) + t140 * t105 + t142 * t94 + t166;
t119 = cos(qJ(1));
t115 = sin(qJ(1));
t75 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t153;
t65 = qJD(2) * t129;
t63 = qJD(2) * t130;
t61 = pkin(8) + t69;
t60 = -pkin(4) - t67;
t55 = t81 * t103;
t40 = -qJDD(5) * mrSges(6,2) + mrSges(6,3) * t56;
t30 = -mrSges(6,1) * t56 + mrSges(6,2) * t57;
t24 = pkin(8) + t29;
t23 = -pkin(4) - t28;
t11 = t29 * qJD(4) + t113 * t65 - t117 * t63;
t10 = t28 * qJD(4) + t113 * t63 + t117 * t65;
t1 = [t121 + (-t115 * mrSges(2,2) + (m(6) * pkin(8) + mrSges(6,3)) * t89 + t126 + t132 * t119 + t122) * g(2) + (-t10 * t75 + t143 * t24 + t141) * t112 + (t10 * t76 + t24 * t40 + (-t24 * t75 - t175) * qJD(5)) * t116 + t11 * t55 + t23 * t30 + (-t10 * t103 - t102 * t29 - t8) * mrSges(5,2) + (-t107 * t68 - t108 * t65 - t33) * mrSges(4,2) + (t102 * t28 - t103 * t11) * mrSges(5,1) + (t107 * t66 + t108 * t63) * mrSges(4,1) + m(4) * (t32 * t66 + t33 * t68 + t42 * t63 + t43 * t65) + m(5) * (t10 * t18 - t11 * t17 + t28 * t9 + t29 * t8) + (mrSges(2,2) * t119 + t132 * t115 + t120 + t183) * g(3) + t24 * t123 + m(6) * (t135 * t10 + t11 * t15 + t23 * t6) + Ifges(2,3) * qJDD(1) + (m(3) * (t114 * t74 + t118 * t73) + (-t107 * t114 - t108 * t149) * mrSges(3,2) + (-qJD(2) * t108 * t114 + t107 * t118) * mrSges(3,1)) * pkin(1); t121 + (t143 * t61 + t168 * t75 + t141) * t112 + (t61 * t40 - t168 * t76 + (-t61 * t75 - t175) * qJD(5) + (g(2) * t90 + g(3) * t89) * mrSges(6,1)) * t116 + (-t102 * t69 + t168 * t103 + t178 - t8) * mrSges(5,2) + (t90 * mrSges(5,1) + t89 * mrSges(6,3) + t122) * g(2) + (-t62 * mrSges(4,1) + t64 * mrSges(4,2) + (mrSges(3,1) * t114 + mrSges(3,2) * t118) * t159) * t108 + (t102 * t67 - t167 * t103) * mrSges(5,1) + (-t107 * t173 - t33) * mrSges(4,2) + t61 * t123 + t107 * mrSges(4,1) * t172 + t167 * t55 + t60 * t30 + (t89 * mrSges(5,1) + t120) * g(3) + ((t90 * pkin(4) + t89 * pkin(8)) * g(2) + (pkin(4) * t89 - t87) * g(3) + t6 * t60 + t167 * t15 - t168 * t135) * m(6) + (-t167 * t17 - t168 * t18 + t67 * t9 + t69 * t8) * m(5) + (-t42 * t62 - t43 * t64 + (t110 * t33 + t111 * t32) * pkin(2)) * m(4); m(6) * (t135 * qJD(5) + t112 * t2 + t116 * t3) + t76 * t147 + t112 * t40 - t75 * t148 + t116 * t41 + (m(4) + m(5)) * qJDD(3) - t146 * g(1); t124 + (t125 + t179) * mrSges(6,3) + (t126 + t165) * g(2) + (t166 + t183) * g(3) + (t103 * mrSges(5,1) - t55) * t18 + (t103 * mrSges(5,2) + t112 * t75 - t116 * t76) * t17 - t6 * t182 + (-t112 * t41 + t116 * t40 + (-t112 * t76 - t116 * t75) * qJD(5) + (-t13 * t147 - t14 * t148 - t170 + t171 + t179) * m(6)) * pkin(8) - m(6) * (t135 * t17 + t15 * t18) - pkin(4) * t30 - t8 * mrSges(5,2); t3 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,5) * t57 + Ifges(6,6) * t56 + Ifges(6,3) * qJDD(5) + g(1) * t81 - t13 * t76 + t14 * t75 + (t45 * t177 + (-t128 / 0.2e1 + t160 * t177) * t103 + t136 * mrSges(6,3) - (t46 + t84) * t116 / 0.2e1 - t184) * t103 + (t178 - t179) * t139;];
tau = t1;
