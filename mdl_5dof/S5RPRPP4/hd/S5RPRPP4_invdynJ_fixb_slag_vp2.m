% Calculate vector of inverse dynamics joint torques for
% S5RPRPP4
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
%   pkin=[a2,a3,a4,a5,d1,d3,theta4]';
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
% Datum: 2019-12-31 18:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRPP4_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP4_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP4_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPP4_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP4_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP4_invdynJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP4_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPP4_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPP4_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:14:18
% EndTime: 2019-12-31 18:14:29
% DurationCPUTime: 6.45s
% Computational Cost: add. (1717->334), mult. (3309->423), div. (0->0), fcn. (1856->8), ass. (0->143)
t220 = mrSges(5,1) + mrSges(6,1);
t219 = mrSges(5,2) - mrSges(6,3);
t218 = m(5) + m(6);
t213 = Ifges(5,1) + Ifges(6,1);
t211 = Ifges(5,5) + Ifges(6,4);
t109 = cos(qJ(3));
t216 = t109 / 0.2e1;
t212 = Ifges(5,4) - Ifges(6,5);
t108 = sin(qJ(1));
t110 = cos(qJ(1));
t206 = -g(1) * t108 + g(2) * t110;
t107 = sin(qJ(3));
t126 = t107 * mrSges(4,1) + t109 * mrSges(4,2);
t103 = qJ(3) + pkin(7);
t96 = sin(t103);
t97 = cos(t103);
t215 = t219 * t97 + t220 * t96 + t126;
t188 = -m(3) - m(4);
t210 = Ifges(6,6) - Ifges(5,6);
t105 = sin(pkin(7));
t154 = cos(pkin(7));
t67 = t105 * t109 + t107 * t154;
t58 = t67 * qJD(1);
t177 = Ifges(6,5) * t58;
t55 = Ifges(5,4) * t58;
t135 = t154 * t109;
t148 = t107 * qJD(1);
t61 = qJD(1) * t135 - t105 * t148;
t209 = t211 * qJD(3) + t213 * t61 + t177 - t55;
t145 = qJD(1) * qJD(3);
t72 = qJDD(1) * t109 - t107 * t145;
t73 = -qJDD(1) * t107 - t109 * t145;
t35 = t105 * t72 - t154 * t73;
t25 = -qJDD(3) * mrSges(5,2) - mrSges(5,3) * t35;
t28 = -mrSges(6,2) * t35 + qJDD(3) * mrSges(6,3);
t208 = t25 + t28;
t36 = t105 * t73 + t154 * t72;
t26 = qJDD(3) * mrSges(5,1) - mrSges(5,3) * t36;
t27 = -qJDD(3) * mrSges(6,1) + t36 * mrSges(6,2);
t207 = t27 - t26;
t179 = mrSges(5,3) * t61;
t180 = mrSges(6,2) * t61;
t167 = -qJD(3) * t220 + t179 + t180;
t170 = t58 * mrSges(5,3);
t181 = mrSges(6,2) * t58;
t47 = qJD(3) * mrSges(6,3) - t181;
t166 = qJD(3) * mrSges(5,2) + t170 - t47;
t146 = qJD(1) * qJD(2);
t82 = qJDD(1) * qJ(2) + t146;
t152 = qJD(3) * t107;
t111 = -pkin(1) - pkin(6);
t80 = qJDD(1) * t111 + qJDD(2);
t81 = qJD(1) * t111 + qJD(2);
t39 = t109 * t80 - t152 * t81;
t151 = qJD(3) * t109;
t40 = t107 * t80 + t81 * t151;
t122 = t40 * t107 + t39 * t109;
t127 = mrSges(4,1) * t109 - mrSges(4,2) * t107;
t164 = Ifges(4,4) * t109;
t204 = qJ(2) * t127 + (-Ifges(4,1) * t107 - t164) * t216;
t147 = t109 * qJD(1);
t161 = t107 * (qJD(3) * mrSges(4,1) - mrSges(4,3) * t147);
t77 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t148;
t203 = (t109 * t77 - t161) * qJD(3);
t52 = (-qJ(4) * qJD(1) + t81) * t107;
t163 = t105 * t52;
t74 = t109 * t81;
t53 = -qJ(4) * t147 + t74;
t50 = qJD(3) * pkin(3) + t53;
t12 = t154 * t50 - t163;
t48 = t154 * t52;
t13 = t105 * t50 + t48;
t144 = qJD(1) * qJD(4);
t11 = qJDD(3) * pkin(3) - qJ(4) * t72 - t109 * t144 + t39;
t15 = qJ(4) * t73 - t107 * t144 + t40;
t4 = -t105 * t15 + t11 * t154;
t5 = t105 * t11 + t154 * t15;
t59 = -qJD(3) * t135 + t105 * t152;
t60 = t67 * qJD(3);
t66 = t105 * t107 - t135;
t202 = t12 * t60 + t13 * t59 + t4 * t66 - t5 * t67;
t1 = qJDD(3) * qJ(5) + qJD(3) * qJD(5) + t5;
t3 = -qJDD(3) * pkin(4) + qJDD(5) - t4;
t8 = -qJD(3) * pkin(4) + qJD(5) - t12;
t9 = qJD(3) * qJ(5) + t13;
t201 = -t1 * t67 - t3 * t66 + t9 * t59 - t8 * t60;
t200 = mrSges(3,2) - mrSges(2,1) - mrSges(6,2) - mrSges(4,3) - mrSges(5,3);
t199 = m(4) * t122 + t109 * (qJDD(3) * mrSges(4,1) - mrSges(4,3) * t72) + t107 * (-qJDD(3) * mrSges(4,2) + mrSges(4,3) * t73);
t129 = pkin(4) * t96 - qJ(5) * t97;
t198 = -m(6) * t129 + mrSges(2,2) - mrSges(3,3) - t215;
t196 = qJD(1) ^ 2;
t194 = -t58 / 0.2e1;
t193 = t58 / 0.2e1;
t190 = t61 / 0.2e1;
t189 = -t66 / 0.2e1;
t178 = Ifges(5,4) * t61;
t176 = pkin(3) * t105;
t101 = t107 * pkin(3);
t165 = Ifges(4,4) * t107;
t155 = t110 * pkin(1) + t108 * qJ(2);
t91 = qJ(2) + t101;
t153 = qJ(4) - t111;
t83 = pkin(3) * t151 + qJD(2);
t149 = qJDD(1) * mrSges(3,2);
t142 = pkin(3) * t147;
t75 = pkin(3) * t148 + qJD(1) * qJ(2) + qJD(4);
t140 = t154 * pkin(3);
t139 = t35 * mrSges(5,1) + t36 * mrSges(5,2);
t138 = t35 * mrSges(6,1) - t36 * mrSges(6,3);
t100 = t110 * qJ(2);
t137 = -pkin(1) * t108 + t100;
t136 = (t146 + t82) * qJ(2);
t134 = -t145 / 0.2e1;
t133 = t153 * t109;
t128 = -g(1) * t110 - g(2) * t108;
t125 = t109 * Ifges(4,1) - t165;
t124 = -t107 * Ifges(4,2) + t164;
t123 = -Ifges(4,5) * t107 - Ifges(4,6) * t109;
t118 = t107 * (-Ifges(4,2) * t109 - t165);
t114 = -qJD(4) * t109 + t152 * t153;
t42 = -pkin(3) * t73 + qJDD(4) + t82;
t106 = -qJ(4) - pkin(6);
t95 = -qJDD(1) * pkin(1) + qJDD(2);
t90 = -t140 - pkin(4);
t88 = qJ(5) + t176;
t76 = t153 * t107;
t69 = t126 * qJD(1);
t63 = Ifges(4,5) * qJD(3) + qJD(1) * t125;
t62 = Ifges(4,6) * qJD(3) + qJD(1) * t124;
t54 = Ifges(6,5) * t61;
t51 = -qJD(3) * t133 - qJD(4) * t107;
t31 = pkin(4) * t67 + qJ(5) * t66 + t91;
t30 = mrSges(5,1) * t58 + mrSges(5,2) * t61;
t29 = mrSges(6,1) * t58 - mrSges(6,3) * t61;
t22 = -t58 * Ifges(5,2) + Ifges(5,6) * qJD(3) + t178;
t21 = Ifges(6,6) * qJD(3) + t58 * Ifges(6,3) + t54;
t20 = pkin(4) * t61 + qJ(5) * t58 + t142;
t19 = t154 * t53 - t163;
t18 = t105 * t53 + t48;
t14 = pkin(4) * t58 - qJ(5) * t61 + t75;
t6 = -pkin(4) * t59 + qJ(5) * t60 + qJD(5) * t66 + t83;
t2 = pkin(4) * t35 - qJ(5) * t36 - qJD(5) * t61 + t42;
t7 = [(Ifges(4,1) * t72 + Ifges(4,4) * t73 + Ifges(4,5) * qJDD(3)) * t216 + t72 * t125 / 0.2e1 + qJD(3) ^ 2 * t123 / 0.2e1 + t73 * t124 / 0.2e1 + t31 * t138 + t91 * t139 + ((Ifges(5,2) + Ifges(6,3)) * t67 + t212 * (t66 / 0.2e1 - t189)) * t35 - t107 * qJDD(3) * Ifges(4,6) + (Ifges(2,3) + Ifges(3,1)) * qJDD(1) + m(3) * (-pkin(1) * t95 + t136) + t67 * (Ifges(6,5) * t36 + Ifges(6,6) * qJDD(3)) / 0.2e1 + t202 * mrSges(5,3) + (-t75 * mrSges(5,2) + t14 * mrSges(6,3) - Ifges(5,4) * t194 - Ifges(6,5) * t193) * t60 + (t126 + 0.2e1 * mrSges(3,3)) * t82 - t67 * (Ifges(5,4) * t36 + Ifges(5,6) * qJDD(3)) / 0.2e1 + m(4) * t136 + (-m(3) * t137 - m(4) * t100 - t218 * (t110 * t101 + t108 * t106 + t137) + t198 * t110 + (-m(4) * t111 - t200) * t108) * g(1) + (t188 * t155 - t218 * (t108 * t101 - t106 * t110 + t155) + (-m(4) * pkin(6) + t200) * t110 + t198 * t108) * g(2) - t209 * t60 / 0.2e1 + (-t210 * t59 - t211 * t60) * qJD(3) / 0.2e1 + (t211 * qJDD(3) + t213 * t36) * t189 + (t212 * t59 - t213 * t60) * t190 + (-t212 * t67 - t213 * t66) * t36 / 0.2e1 + t83 * t30 + t95 * mrSges(3,2) + qJD(2) * t69 + qJ(2) * (-mrSges(4,1) * t73 + mrSges(4,2) * t72) + t2 * (mrSges(6,1) * t67 + mrSges(6,3) * t66) + t42 * (mrSges(5,1) * t67 - mrSges(5,2) * t66) + (Ifges(4,5) * t109 + t210 * t67 - t211 * t66) * qJDD(3) / 0.2e1 + t6 * t29 + m(5) * (t42 * t91 + t75 * t83) + m(6) * (t14 * t6 + t2 * t31) + t199 * t111 + (-t75 * mrSges(5,1) + t22 / 0.2e1 - t21 / 0.2e1 - t14 * mrSges(6,1) - Ifges(6,3) * t193 + Ifges(5,2) * t194) * t59 + t201 * mrSges(6,2) - pkin(1) * t149 - t62 * t151 / 0.2e1 - t63 * t152 / 0.2e1 + (m(5) * t13 + m(6) * t9 - t166) * (t105 * t114 + t154 * t51) + (-m(5) * t12 + m(6) * t8 + t167) * (t105 * t51 - t114 * t154) + (-m(5) * t4 + m(6) * t3 + t207) * (-t105 * t76 + t133 * t154) + (m(5) * t5 + m(6) * t1 + t208) * (-t105 * t133 - t154 * t76) + t118 * t134 + t111 * t203 + t204 * t145 - t122 * mrSges(4,3) - t107 * (Ifges(4,4) * t72 + Ifges(4,2) * t73) / 0.2e1; t149 + t208 * t67 + t207 * t66 + t167 * t60 + t166 * t59 + t203 + (qJ(2) * t188 - mrSges(3,3)) * t196 + (-m(5) * t75 - m(6) * t14 - t29 - t30 - t69) * qJD(1) + m(3) * t95 - m(5) * t202 - m(6) * t201 + t206 * (-t188 + t218) + t199; ((t105 * t5 + t154 * t4) * pkin(3) + t12 * t18 - t13 * t19 - t75 * t142) * m(5) + (Ifges(5,3) + Ifges(4,3) + Ifges(6,2)) * qJDD(3) + (-Ifges(5,2) * t61 + t209 - t55) * t193 + t210 * t35 + t211 * t36 - (t210 * t61 - t211 * t58) * qJD(3) / 0.2e1 - (-t213 * t58 - t178 + t21 + t54) * t61 / 0.2e1 + t206 * (t127 + t218 * pkin(3) * t109 + (m(6) * pkin(4) + t220) * t97 + (m(6) * qJ(5) - t219) * t96) + t88 * t28 + t90 * t27 + Ifges(4,5) * t72 + Ifges(4,6) * t73 - t75 * (mrSges(5,1) * t61 - mrSges(5,2) * t58) - t14 * (mrSges(6,1) * t61 + mrSges(6,3) * t58) - t12 * t170 + qJD(5) * t47 + t39 * mrSges(4,1) - t40 * mrSges(4,2) - t20 * t29 + t1 * mrSges(6,3) - t3 * mrSges(6,1) + t4 * mrSges(5,1) - t5 * mrSges(5,2) + t22 * t190 + (Ifges(6,3) * t61 - t177) * t194 + t25 * t176 + t13 * t179 + t9 * t180 + t8 * t181 + t26 * t140 + t81 * t161 + t63 * t148 / 0.2e1 + t62 * t147 / 0.2e1 - t30 * t142 + (-m(6) * (-t101 - t129) + m(5) * t101 + t215) * g(3) + t166 * t19 - t167 * t18 - t77 * t74 + (t1 * t88 - t14 * t20 - t18 * t8 + t3 * t90 + (-t19 + qJD(5)) * t9) * m(6) + t123 * t134 + (t118 / 0.2e1 - t204) * t196; -t166 * t58 - t167 * t61 + t138 + t139 + (t58 * t9 - t61 * t8 + t128 + t2) * m(6) + (t12 * t61 + t13 * t58 + t128 + t42) * m(5); -qJD(3) * t47 + t61 * t29 + (-g(3) * t96 - t9 * qJD(3) + t14 * t61 - t206 * t97 + t3) * m(6) + t27;];
tau = t7;
