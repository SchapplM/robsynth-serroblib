% Calculate vector of inverse dynamics joint torques for
% S5RPPRR1
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
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
% Datum: 2019-12-05 17:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPRR1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR1_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR1_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR1_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR1_invdynJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR1_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR1_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR1_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:37:59
% EndTime: 2019-12-05 17:38:08
% DurationCPUTime: 3.47s
% Computational Cost: add. (1401->276), mult. (2598->359), div. (0->0), fcn. (1341->8), ass. (0->129)
t181 = m(4) + m(5);
t184 = -m(6) - t181;
t185 = m(3) - t184;
t102 = sin(qJ(4));
t183 = -t102 / 0.2e1;
t105 = cos(qJ(4));
t163 = t105 / 0.2e1;
t122 = t102 * mrSges(5,1) + t105 * mrSges(5,2);
t98 = qJ(4) + qJ(5);
t85 = sin(t98);
t86 = cos(t98);
t124 = t85 * mrSges(6,1) + t86 * mrSges(6,2);
t182 = t124 + t122;
t180 = mrSges(3,2) - mrSges(4,3);
t179 = -mrSges(4,2) - mrSges(3,3);
t94 = qJD(1) * qJD(2);
t178 = qJDD(1) * qJ(2) + t94;
t101 = sin(qJ(5));
t104 = cos(qJ(5));
t140 = qJD(1) * t102;
t83 = qJD(1) * qJ(2) + qJD(3);
t69 = -qJD(1) * pkin(6) + t83;
t40 = -pkin(7) * t140 + t102 * t69;
t145 = t104 * t40;
t139 = qJD(1) * t105;
t41 = -pkin(7) * t139 + t105 * t69;
t37 = qJD(4) * pkin(4) + t41;
t18 = t101 * t37 + t145;
t177 = qJD(5) * t18;
t136 = qJD(4) * t102;
t66 = qJDD(3) + t178;
t55 = -qJDD(1) * pkin(6) + t66;
t30 = t105 * t55 - t136 * t69;
t135 = qJD(4) * t105;
t31 = t102 * t55 + t69 * t135;
t118 = t102 * t31 + t105 * t30;
t103 = sin(qJ(1));
t106 = cos(qJ(1));
t176 = g(1) * t106 + g(2) * t103;
t149 = Ifges(5,4) * t105;
t150 = Ifges(5,4) * t102;
t174 = (-Ifges(5,2) * t105 - t150) * t183 + (-Ifges(5,1) * t102 - t149) * t163;
t141 = t104 * t105;
t115 = t101 * t102 - t141;
t148 = t101 * t40;
t17 = t104 * t37 - t148;
t131 = qJD(1) * qJD(4);
t58 = qJDD(1) * t105 - t102 * t131;
t19 = qJDD(4) * pkin(4) - pkin(7) * t58 + t30;
t59 = -qJDD(1) * t102 - t105 * t131;
t22 = pkin(7) * t59 + t31;
t2 = qJD(5) * t17 + t101 * t19 + t104 * t22;
t134 = qJD(5) * t101;
t92 = qJD(4) + qJD(5);
t28 = -t101 * t136 - t102 * t134 + t141 * t92;
t53 = -t101 * t105 - t104 * t102;
t111 = t53 * qJD(5);
t29 = qJD(4) * t53 + t111;
t3 = -t101 * t22 + t104 * t19 - t177;
t173 = t115 * t3 - t17 * t29 - t18 * t28 + t2 * t53;
t123 = mrSges(5,1) * t105 - mrSges(5,2) * t102;
t100 = pkin(1) + qJ(3);
t70 = qJD(1) * t100 - qJD(2);
t172 = t70 * t123 + qJD(4) * (-Ifges(5,5) * t102 - Ifges(5,6) * t105) / 0.2e1;
t171 = m(5) * t118 + t105 * (qJDD(4) * mrSges(5,1) - mrSges(5,3) * t58) + t102 * (-qJDD(4) * mrSges(5,2) + mrSges(5,3) * t59);
t170 = mrSges(2,1) - t180 + t182;
t169 = m(5) * pkin(6) - m(6) * (-pkin(7) - pkin(6)) + mrSges(2,2) + mrSges(5,3) + mrSges(6,3) + t179;
t50 = -t101 * t140 + t104 * t139;
t167 = t50 / 0.2e1;
t99 = qJ(2) - pkin(6);
t162 = pkin(7) - t99;
t161 = mrSges(6,1) * t86;
t160 = mrSges(6,2) * t85;
t49 = -t101 * t139 - t104 * t140;
t159 = mrSges(6,3) * t49;
t158 = pkin(4) * t102;
t153 = t50 * mrSges(6,3);
t152 = t50 * Ifges(6,4);
t151 = t106 * pkin(1) + t103 * qJ(2);
t142 = qJDD(1) * pkin(1);
t138 = qJD(2) * t102;
t137 = qJD(2) * t105;
t132 = qJD(1) * qJD(3);
t27 = -mrSges(6,1) * t49 + mrSges(6,2) * t50;
t78 = t100 + t158;
t51 = qJD(1) * t78 - qJD(2);
t128 = -m(6) * t51 - t27;
t62 = t162 * t105;
t127 = (t102 ^ 2 + t105 ^ 2) * t69;
t126 = -t59 * mrSges(5,1) + t58 * mrSges(5,2);
t15 = qJD(1) * t111 + t101 * t59 + t104 * t58;
t16 = qJD(1) * qJD(5) * t115 - t101 * t58 + t104 * t59;
t125 = -t16 * mrSges(6,1) + t15 * mrSges(6,2);
t121 = t105 * Ifges(5,1) - t150;
t120 = -t102 * Ifges(5,2) + t149;
t61 = t162 * t102;
t33 = -t101 * t62 - t104 * t61;
t32 = t101 * t61 - t104 * t62;
t63 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t140;
t64 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t139;
t117 = t102 * t64 - t105 * t63;
t56 = qJDD(1) * t100 - qJDD(2) + t132;
t116 = qJD(3) * t70 + t100 * t56;
t110 = -t123 + t160;
t24 = t49 * Ifges(6,2) + t92 * Ifges(6,6) + t152;
t44 = Ifges(6,4) * t49;
t25 = t50 * Ifges(6,1) + t92 * Ifges(6,5) + t44;
t91 = qJDD(4) + qJDD(5);
t109 = t3 * mrSges(6,1) - t2 * mrSges(6,2) + t17 * t159 + t24 * t167 - t51 * (mrSges(6,1) * t50 + mrSges(6,2) * t49) - t50 * (Ifges(6,1) * t49 - t152) / 0.2e1 + Ifges(6,6) * t16 + Ifges(6,5) * t15 - t92 * (Ifges(6,5) * t49 - Ifges(6,6) * t50) / 0.2e1 + Ifges(6,3) * t91 - (-Ifges(6,2) * t50 + t25 + t44) * t49 / 0.2e1;
t108 = qJD(1) ^ 2;
t84 = qJDD(2) - t142;
t72 = pkin(4) * t135 + qJD(3);
t68 = t106 * t161;
t67 = t103 * t161;
t57 = t122 * qJD(1);
t48 = Ifges(5,5) * qJD(4) + qJD(1) * t121;
t47 = Ifges(5,6) * qJD(4) + qJD(1) * t120;
t39 = -qJD(4) * t62 + t138;
t38 = t136 * t162 + t137;
t36 = mrSges(6,1) * t92 - t153;
t35 = -mrSges(6,2) * t92 + t159;
t34 = -pkin(4) * t59 + t56;
t21 = t104 * t41 - t148;
t20 = -t101 * t41 - t145;
t10 = -mrSges(6,2) * t91 + mrSges(6,3) * t16;
t9 = mrSges(6,1) * t91 - mrSges(6,3) * t15;
t5 = -qJD(5) * t33 - t101 * t39 + t104 * t38;
t4 = qJD(5) * t32 + t101 * t38 + t104 * t39;
t1 = [m(3) * (-pkin(1) * t84 + (t178 + t94) * qJ(2)) + (Ifges(6,1) * t29 - Ifges(6,4) * t28) * t167 + (mrSges(4,3) * t100 + Ifges(3,1) + Ifges(4,1) + Ifges(2,3)) * qJDD(1) - (mrSges(6,2) * t34 + Ifges(6,1) * t15 + Ifges(6,4) * t16 + Ifges(6,5) * t91) * t115 + t56 * t122 + (Ifges(5,1) * t58 + Ifges(5,4) * t59) * t163 + (0.2e1 * Ifges(5,5) * t163 - Ifges(5,6) * t102) * qJDD(4) + (t66 + t178) * mrSges(4,2) - t118 * mrSges(5,3) + ((m(3) * pkin(1) + m(6) * t78 + t181 * t100 + t170) * t103 + (-t185 * qJ(2) + t169) * t106) * g(1) + m(5) * (qJD(2) * t127 + t116) + (-m(3) * t151 + t184 * (t106 * qJ(3) + t151) + (-m(6) * t158 - t170) * t106 + t169 * t103) * g(2) + (t56 + t132) * mrSges(4,3) + t174 * t131 + (t135 * t63 - t136 * t64 + t171) * t99 + t172 * qJD(4) + t173 * mrSges(6,3) + (t84 - t142) * mrSges(3,2) + t92 * (Ifges(6,5) * t29 - Ifges(6,6) * t28) / 0.2e1 + t72 * t27 + t51 * (mrSges(6,1) * t28 + mrSges(6,2) * t29) + qJD(3) * t57 + t49 * (Ifges(6,4) * t29 - Ifges(6,2) * t28) / 0.2e1 + t29 * t25 / 0.2e1 + t32 * t9 + t33 * t10 + t4 * t35 + t5 * t36 - t28 * t24 / 0.2e1 + m(6) * (t17 * t5 + t18 * t4 + t2 * t33 + t3 * t32 + t34 * t78 + t51 * t72) + (-mrSges(6,1) * t34 + Ifges(6,4) * t15 + Ifges(6,2) * t16 + Ifges(6,6) * t91) * t53 + 0.2e1 * t178 * mrSges(3,3) + (Ifges(5,4) * t58 + Ifges(5,2) * t59) * t183 + t64 * t137 + t63 * t138 + m(4) * (qJ(2) * t66 + qJD(2) * t83 + t116) + t59 * t120 / 0.2e1 + t58 * t121 / 0.2e1 + t78 * t125 + t100 * t126 - t47 * t135 / 0.2e1 - t48 * t136 / 0.2e1; t49 * t35 - t50 * t36 + t180 * qJDD(1) + (-m(3) * qJ(2) + t179) * t108 + (-m(4) * t83 - t102 * t63 - t105 * t64) * qJD(1) + m(3) * t84 - m(4) * t56 - t125 - t126 + (-g(1) * t103 + g(2) * t106) * t185 + (-t17 * t50 + t18 * t49 - t34) * m(6) + (-qJD(1) * t127 - t56) * m(5); qJDD(1) * mrSges(4,2) - t108 * mrSges(4,3) - t53 * t10 + t28 * t35 + t29 * t36 - t115 * t9 - t117 * qJD(4) + m(4) * t66 - m(6) * t173 + (-t181 * t70 + t128 - t57) * qJD(1) + t176 * t184 + t171; (t47 * t163 + t102 * t48 / 0.2e1 - t174 * qJD(1) + t128 * t105 * pkin(4) - t172) * qJD(1) + ((g(3) * t102 - t176 * t105 - t17 * t134) * m(6) + (m(6) * t2 - qJD(5) * t36 + t10) * t101 + (t9 + t35 * qJD(5) + (t3 + t177) * m(6)) * t104) * pkin(4) + (t106 * t110 - t68) * g(1) + (t103 * t110 - t67) * g(2) + t182 * g(3) - m(6) * (t17 * t20 + t18 * t21) + t117 * t69 + t109 + t18 * t153 + Ifges(5,6) * t59 + Ifges(5,5) * t58 + t30 * mrSges(5,1) - t31 * mrSges(5,2) - t21 * t35 - t20 * t36 + Ifges(5,3) * qJDD(4); (t36 + t153) * t18 - g(1) * (-t106 * t160 + t68) + t109 - g(2) * (-t103 * t160 + t67) + g(3) * t124 - t17 * t35;];
tau = t1;
