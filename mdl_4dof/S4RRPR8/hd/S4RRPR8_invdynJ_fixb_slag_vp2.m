% Calculate vector of inverse dynamics joint torques for
% S4RRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RRPR8_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR8_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR8_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR8_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR8_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR8_invdynJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR8_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR8_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR8_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:07:48
% EndTime: 2019-12-31 17:08:00
% DurationCPUTime: 6.66s
% Computational Cost: add. (1256->321), mult. (2799->431), div. (0->0), fcn. (1531->6), ass. (0->153)
t106 = sin(qJ(2));
t109 = cos(qJ(2));
t174 = t109 * mrSges(4,3);
t94 = t106 * qJ(3);
t144 = pkin(1) + t94;
t97 = t109 * pkin(2);
t122 = -t144 - t97;
t50 = t122 * qJD(1);
t163 = t106 * qJD(1);
t91 = pkin(5) * t163;
t203 = qJD(3) + t91;
t68 = -qJD(2) * pkin(2) + t203;
t104 = qJD(2) * qJ(3);
t162 = t109 * qJD(1);
t92 = pkin(5) * t162;
t72 = t92 + t104;
t221 = -pkin(5) * (-t72 * t106 + t68 * t109) * m(4) - t50 * (t106 * mrSges(4,1) - t174);
t219 = Ifges(4,4) + Ifges(3,5);
t218 = Ifges(3,6) - Ifges(4,6);
t105 = sin(qJ(4));
t108 = cos(qJ(4));
t48 = -t105 * t162 + t108 * t163;
t196 = t48 / 0.2e1;
t220 = -pkin(6) * t163 + t203;
t208 = m(4) + m(5);
t177 = Ifges(4,5) * t109;
t132 = t106 * Ifges(4,1) - t177;
t90 = Ifges(3,4) * t162;
t217 = Ifges(3,1) * t163 + qJD(1) * t132 + qJD(2) * t219 + t90;
t153 = mrSges(4,2) * t163;
t216 = mrSges(3,3) * t163 + t153 + (-mrSges(3,1) - mrSges(4,1)) * qJD(2);
t215 = -t218 * t106 + t109 * t219;
t133 = t109 * mrSges(4,1) + t106 * mrSges(4,3);
t135 = mrSges(3,1) * t109 - mrSges(3,2) * t106;
t214 = -t133 - t135;
t161 = qJD(1) * qJD(2);
t65 = qJDD(1) * t106 + t109 * t161;
t64 = -t109 * qJDD(1) + t106 * t161;
t53 = t64 * pkin(5);
t54 = t65 * pkin(5);
t213 = t106 * t54 - t109 * t53;
t28 = qJDD(2) * qJ(3) + qJD(2) * qJD(3) - t53;
t151 = qJDD(3) + t54;
t32 = -qJDD(2) * pkin(2) + t151;
t212 = t106 * t32 + t109 * t28;
t101 = -qJD(2) + qJD(4);
t211 = -mrSges(2,1) + t214;
t107 = sin(qJ(1));
t110 = cos(qJ(1));
t202 = g(1) * t110 + g(2) * t107;
t210 = m(5) * pkin(6) + mrSges(2,2) - mrSges(3,3) + mrSges(5,3);
t194 = pkin(2) + pkin(3);
t29 = -qJD(2) * t194 + t220;
t61 = -pkin(6) * t162 + t92;
t49 = t104 + t61;
t10 = -t105 * t49 + t108 * t29;
t16 = -pkin(6) * t65 - qJDD(2) * t194 + t151;
t17 = pkin(6) * t64 + t28;
t1 = qJD(4) * t10 + t105 * t16 + t108 * t17;
t100 = -qJDD(2) + qJDD(4);
t41 = Ifges(5,4) * t48;
t47 = -t105 * t163 - t108 * t162;
t12 = Ifges(5,2) * t47 + Ifges(5,6) * t101 + t41;
t40 = Ifges(5,4) * t47;
t13 = Ifges(5,1) * t48 + Ifges(5,5) * t101 + t40;
t11 = t105 * t29 + t108 * t49;
t2 = -qJD(4) * t11 - t105 * t17 + t108 * t16;
t201 = t109 * t194 + t144;
t27 = t201 * qJD(1);
t123 = t105 * t106 + t108 * t109;
t115 = t123 * qJD(4);
t8 = -qJD(1) * t115 + t105 * t64 + t108 * t65;
t124 = t105 * t109 - t106 * t108;
t116 = t124 * qJD(4);
t9 = qJD(1) * t116 - t105 * t65 + t108 * t64;
t209 = t2 * mrSges(5,1) - t1 * mrSges(5,2) + Ifges(5,5) * t8 + Ifges(5,6) * t9 + Ifges(5,3) * t100 - t101 * (Ifges(5,5) * t47 - Ifges(5,6) * t48) / 0.2e1 - t27 * (mrSges(5,1) * t48 + mrSges(5,2) * t47) - (-Ifges(5,2) * t48 + t13 + t40) * t47 / 0.2e1 + (-Ifges(5,1) * t47 + t12 + t41) * t196;
t66 = -t105 * qJ(3) - t108 * t194;
t207 = qJD(4) * t66 - t105 * t61 + t108 * t220;
t67 = t108 * qJ(3) - t105 * t194;
t206 = -qJD(4) * t67 - t105 * t220 - t108 * t61;
t150 = mrSges(5,1) * t123 - mrSges(5,2) * t124;
t193 = pkin(5) - pkin(6);
t192 = t106 / 0.2e1;
t189 = pkin(5) * t106;
t188 = pkin(5) * t109;
t185 = t11 * t48;
t184 = t47 * mrSges(5,3);
t182 = t97 + t94;
t181 = t110 * pkin(1) + t107 * pkin(5);
t180 = Ifges(3,4) * t106;
t179 = Ifges(3,4) * t109;
t178 = Ifges(4,5) * t106;
t170 = qJ(3) * t109;
t169 = qJDD(1) * pkin(1);
t168 = t106 * t110;
t167 = t109 * t110;
t166 = qJD(2) * t106;
t164 = qJD(3) * t106;
t157 = t109 * pkin(3) + t182;
t154 = t194 * t106;
t77 = t193 * t109;
t152 = mrSges(4,2) * t162;
t143 = -t161 / 0.2e1;
t141 = pkin(2) * t167 + qJ(3) * t168 + t181;
t38 = -qJDD(2) * mrSges(4,1) + t65 * mrSges(4,2);
t34 = t124 * t107;
t35 = t123 * t107;
t140 = -t34 * mrSges(5,1) - t35 * mrSges(5,2);
t36 = t105 * t167 - t108 * t168;
t37 = t123 * t110;
t139 = -t36 * mrSges(5,1) - t37 * mrSges(5,2);
t134 = mrSges(3,1) * t106 + mrSges(3,2) * t109;
t131 = t109 * Ifges(3,2) + t180;
t128 = pkin(2) * t106 - t170;
t24 = -mrSges(5,2) * t101 + t184;
t25 = mrSges(5,1) * t101 - mrSges(5,3) * t48;
t126 = -t105 * t25 + t108 * t24;
t76 = t193 * t106;
t22 = -t105 * t77 + t108 * t76;
t23 = t105 * t76 + t108 * t77;
t121 = pkin(1) * t134;
t120 = -t154 + t170;
t118 = t106 * (Ifges(3,1) * t109 - t180);
t117 = t109 * (Ifges(4,3) * t106 + t177);
t114 = t174 + (-m(4) * pkin(2) - mrSges(4,1)) * t106;
t113 = qJ(3) * t65 + qJD(3) * t163 + t169;
t89 = Ifges(4,5) * t163;
t83 = qJ(3) * t167;
t82 = t107 * t170;
t74 = qJD(2) * mrSges(4,3) + t152;
t73 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t162;
t69 = -pkin(1) - t182;
t63 = qJD(2) * t77;
t62 = t193 * t166;
t59 = t128 * qJD(1);
t58 = t133 * qJD(1);
t51 = pkin(1) + t157;
t44 = Ifges(3,6) * qJD(2) + qJD(1) * t131;
t43 = Ifges(4,6) * qJD(2) - Ifges(4,3) * t162 + t89;
t42 = qJD(2) * t128 - t164;
t39 = -mrSges(4,2) * t64 + qJDD(2) * mrSges(4,3);
t33 = t120 * qJD(1);
t26 = qJD(2) * t120 + t164;
t21 = qJD(2) * t123 - t115;
t20 = -qJD(2) * t124 + t116;
t15 = -mrSges(5,1) * t47 + mrSges(5,2) * t48;
t14 = pkin(2) * t64 - t113;
t7 = -t194 * t64 + t113;
t6 = -mrSges(5,2) * t100 + mrSges(5,3) * t9;
t5 = mrSges(5,1) * t100 - mrSges(5,3) * t8;
t4 = -qJD(4) * t23 + t105 * t62 + t108 * t63;
t3 = qJD(4) * t22 + t105 * t63 - t108 * t62;
t18 = [(-t202 + t212) * mrSges(4,2) + (-t10 * t21 + t11 * t20) * mrSges(5,3) + (t35 * mrSges(5,1) - t34 * mrSges(5,2) + (m(3) * pkin(1) - m(4) * t122 + t201 * m(5) - t211) * t107 + ((-m(3) - t208) * pkin(5) + t210) * t110) * g(1) + (t219 * t106 + t218 * t109) * qJDD(2) / 0.2e1 + ((Ifges(3,1) + Ifges(4,1)) * t65 + t219 * qJDD(2)) * t192 + (t189 * t65 + t213) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + t213 * pkin(5)) + (-m(4) * t141 - m(5) * (pkin(3) * t167 + t141) - t37 * mrSges(5,1) + t36 * mrSges(5,2) - m(3) * t181 + t211 * t110 + t210 * t107) * g(2) + m(4) * (pkin(5) * t212 + t14 * t69 + t42 * t50) + (Ifges(5,1) * t21 + Ifges(5,4) * t20) * t196 + t39 * t188 + (-qJDD(2) * mrSges(3,1) + t38) * t189 + (-t72 * mrSges(4,2) - t44 / 0.2e1 + t43 / 0.2e1 + (-t74 - t73) * pkin(5)) * t166 + ((t216 * pkin(5) + t217 / 0.2e1 + t68 * mrSges(4,2)) * t109 + t215 * qJD(2) / 0.2e1 - t221) * qJD(2) + (-mrSges(3,3) * t188 - pkin(1) * mrSges(3,1) + t69 * mrSges(4,1) - t131 / 0.2e1 + t178 / 0.2e1 + (Ifges(4,5) - Ifges(3,4)) * t192 + (-Ifges(3,2) / 0.2e1 - Ifges(4,3)) * t109) * t64 - t109 * (Ifges(4,5) * t65 + Ifges(4,6) * qJDD(2)) / 0.2e1 + t101 * (Ifges(5,5) * t21 + Ifges(5,6) * t20) / 0.2e1 - t42 * t58 + t51 * (-mrSges(5,1) * t9 + mrSges(5,2) * t8) + t47 * (Ifges(5,4) * t21 + Ifges(5,2) * t20) / 0.2e1 + t27 * (-mrSges(5,1) * t20 + mrSges(5,2) * t21) + t21 * t13 / 0.2e1 + t22 * t5 + t23 * t6 + t3 * t24 + t4 * t25 + t26 * t15 + t20 * t12 / 0.2e1 + (t106 * (Ifges(4,1) * t109 + t178) + t109 * (-Ifges(3,2) * t106 + t179) + t118) * t161 / 0.2e1 + (t106 * Ifges(3,1) + t132 + t179) * t65 / 0.2e1 + (-mrSges(5,3) * t1 - Ifges(5,4) * t8 - Ifges(5,2) * t9 - Ifges(5,6) * t100) * t123 + (mrSges(5,3) * t2 - Ifges(5,1) * t8 - Ifges(5,4) * t9 - Ifges(5,5) * t100) * t124 + t135 * t169 - t69 * mrSges(4,3) * t65 + m(5) * (t1 * t23 + t10 * t4 + t11 * t3 + t2 * t22 + t26 * t27 + t51 * t7) + t117 * t143 - t14 * t133 + (-pkin(1) * t65 - qJDD(2) * t188) * mrSges(3,2) + t7 * t150 - t121 * t161 + t109 * (Ifges(3,4) * t65 + Ifges(3,6) * qJDD(2)) / 0.2e1 + Ifges(2,3) * qJDD(1); -t209 - t216 * t92 - (-Ifges(3,2) * t163 + t217 + t90) * t162 / 0.2e1 - t218 * t64 + t219 * t65 + (-m(4) * t182 - m(5) * t157 - t150 + t214) * g(3) + t215 * t143 + (Ifges(4,2) + Ifges(3,3)) * qJDD(2) + t206 * t25 + t207 * t24 + (-t27 * t33 + t1 * t67 + t2 * t66 - g(1) * (-t110 * t154 + t83) - g(2) * (-t107 * t154 + t82) + t207 * t11 + t206 * t10) * m(5) + ((t117 / 0.2e1 - t118 / 0.2e1 + t121) * qJD(1) + t221) * qJD(1) + (-t114 * t107 + t140) * g(2) + t66 * t5 + t67 * t6 + t59 * t58 + t53 * mrSges(3,2) - t54 * mrSges(3,1) - pkin(2) * t38 + qJ(3) * t39 + t28 * mrSges(4,3) - t32 * mrSges(4,1) - t33 * t15 + t72 * t153 + t73 * t91 + (-t114 * t110 + t139) * g(1) + t203 * t74 + t202 * t134 + (-pkin(2) * t32 - g(1) * t83 - g(2) * t82 + qJ(3) * t28 + qJD(3) * t72 - t50 * t59) * m(4) - (Ifges(4,1) * t162 + t43 + t89) * t163 / 0.2e1 - t68 * t152 + t44 * t163 / 0.2e1 - t10 * t184 - mrSges(5,3) * t185; t105 * t6 + t108 * t5 + t126 * qJD(4) + t208 * t109 * g(3) + (-t126 - t74) * qJD(2) + ((-t15 - t58) * qJD(1) - t208 * t202) * t106 + t38 + (t1 * t105 + t108 * t2 - t163 * t27 + t101 * (-t10 * t105 + t108 * t11)) * m(5) + (-qJD(2) * t72 + t163 * t50 + t32) * m(4); (t10 * t47 + t185) * mrSges(5,3) + g(3) * t150 - g(1) * t139 - g(2) * t140 - t10 * t24 + t11 * t25 + t209;];
tau = t18;
