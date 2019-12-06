% Calculate vector of inverse dynamics joint torques for
% S5RPPRR2
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
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
% Datum: 2019-12-05 17:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPRR2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR2_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR2_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR2_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR2_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR2_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR2_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR2_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:39:34
% EndTime: 2019-12-05 17:39:46
% DurationCPUTime: 4.80s
% Computational Cost: add. (3070->323), mult. (6277->424), div. (0->0), fcn. (4247->12), ass. (0->142)
t134 = sin(pkin(8));
t135 = cos(pkin(8));
t175 = t134 ^ 2 + t135 ^ 2;
t160 = t175 * mrSges(4,3);
t137 = -pkin(1) - qJ(3);
t223 = -qJD(1) * qJD(3) + qJDD(1) * t137;
t138 = sin(qJ(5));
t141 = cos(qJ(5));
t139 = sin(qJ(4));
t142 = cos(qJ(4));
t178 = t135 * t142;
t149 = t134 * t139 - t178;
t90 = -t142 * t134 - t139 * t135;
t59 = -t138 * t149 - t141 * t90;
t172 = qJD(4) * t142;
t173 = qJD(4) * t139;
t84 = -t134 * t173 + t135 * t172;
t85 = t90 * qJD(4);
t145 = -qJD(5) * t59 - t138 * t84 + t141 * t85;
t211 = t138 * t90 - t141 * t149;
t26 = qJD(5) * t211 + t138 * t85 + t141 * t84;
t101 = qJD(1) * t137 + qJD(2);
t158 = -pkin(6) * qJD(1) + t101;
t77 = t158 * t134;
t78 = t158 * t135;
t43 = t139 * t78 + t142 * t77;
t82 = t90 * qJD(1);
t34 = pkin(7) * t82 + t43;
t184 = t138 * t34;
t42 = -t139 * t77 + t142 * t78;
t174 = qJD(1) * t134;
t83 = qJD(1) * t178 - t139 * t174;
t33 = -pkin(7) * t83 + t42;
t32 = qJD(4) * pkin(4) + t33;
t8 = t141 * t32 - t184;
t181 = t141 * t34;
t9 = t138 * t32 + t181;
t222 = -t145 * t8 - t26 * t9;
t125 = qJDD(4) + qJDD(5);
t130 = qJD(4) + qJD(5);
t162 = -t138 * t83 + t141 * t82;
t55 = qJD(1) * t85 - qJDD(1) * t149;
t56 = qJD(1) * qJD(4) * t149 + qJDD(1) * t90;
t17 = qJD(5) * t162 + t138 * t56 + t141 * t55;
t52 = t138 * t82 + t141 * t83;
t18 = -qJD(5) * t52 - t138 * t55 + t141 * t56;
t195 = Ifges(6,4) * t52;
t93 = qJDD(2) + t223;
t157 = -pkin(6) * qJDD(1) + t93;
t68 = t157 * t134;
t69 = t157 * t135;
t20 = -qJD(4) * t43 - t139 * t68 + t142 * t69;
t6 = qJDD(4) * pkin(4) - pkin(7) * t55 + t20;
t19 = t139 * t69 + t142 * t68 + t78 * t172 - t173 * t77;
t7 = pkin(7) * t56 + t19;
t2 = qJD(5) * t8 + t138 * t6 + t141 * t7;
t44 = Ifges(6,4) * t162;
t24 = Ifges(6,1) * t52 + Ifges(6,5) * t130 + t44;
t3 = -qJD(5) * t9 - t138 * t7 + t141 * t6;
t116 = qJD(1) * qJ(2) + qJD(3);
t97 = pkin(3) * t174 + t116;
t63 = -pkin(4) * t82 + t97;
t221 = t3 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,5) * t17 + Ifges(6,6) * t18 + Ifges(6,3) * t125 - (Ifges(6,5) * t162 - Ifges(6,6) * t52) * t130 / 0.2e1 + (t162 * t8 + t52 * t9) * mrSges(6,3) - (-Ifges(6,2) * t52 + t24 + t44) * t162 / 0.2e1 - t63 * (mrSges(6,1) * t52 + mrSges(6,2) * t162) - (Ifges(6,1) * t162 - t195) * t52 / 0.2e1;
t129 = pkin(8) + qJ(4);
t118 = sin(t129);
t119 = cos(t129);
t120 = qJ(5) + t129;
t111 = sin(t120);
t112 = cos(t120);
t152 = t111 * mrSges(6,1) + t112 * mrSges(6,2);
t220 = t118 * mrSges(5,1) + t119 * mrSges(5,2) + t152;
t23 = Ifges(6,2) * t162 + Ifges(6,6) * t130 + t195;
t218 = t23 / 0.2e1;
t171 = -m(4) - m(5) - m(6);
t215 = m(3) - t171;
t132 = qJD(1) * qJD(2);
t102 = -qJDD(1) * qJ(2) - t132;
t140 = sin(qJ(1));
t143 = cos(qJ(1));
t210 = -g(1) * t140 + g(2) * t143;
t206 = -t149 * t20 - t19 * t90 + t42 * t85 + t43 * t84;
t205 = mrSges(2,1) + mrSges(4,3) + mrSges(5,3) + mrSges(6,3) - mrSges(3,2);
t121 = t134 * pkin(3);
t155 = mrSges(4,1) * t134 + mrSges(4,2) * t135;
t204 = -m(5) * t121 + mrSges(2,2) - mrSges(3,3) - m(6) * (pkin(4) * t118 + t121) - t155 - t220;
t203 = m(4) * t116 + m(5) * t97 - mrSges(5,1) * t82 + mrSges(5,2) * t83 + t155 * qJD(1);
t200 = t52 / 0.2e1;
t198 = t83 / 0.2e1;
t196 = Ifges(5,4) * t83;
t136 = -pkin(6) - qJ(3);
t187 = -pkin(6) + t137;
t94 = t187 * t134;
t95 = t187 * t135;
t62 = t139 * t95 + t142 * t94;
t186 = mrSges(6,1) * t112;
t185 = mrSges(6,2) * t111;
t180 = qJDD(1) * pkin(1);
t106 = qJ(2) + t121;
t168 = qJDD(1) * t135;
t169 = qJDD(1) * t134;
t177 = mrSges(4,1) * t169 + mrSges(4,2) * t168;
t164 = -t56 * mrSges(5,1) + t55 * mrSges(5,2);
t163 = -t18 * mrSges(6,1) + t17 * mrSges(6,2);
t100 = qJDD(3) - t102;
t61 = -t139 * t94 + t142 * t95;
t159 = t175 * t93;
t156 = t175 * t101;
t86 = pkin(3) * t169 + t100;
t154 = mrSges(5,1) * t119 - mrSges(5,2) * t118;
t40 = pkin(7) * t149 + t61;
t41 = pkin(7) * t90 + t62;
t21 = -t138 * t41 + t141 * t40;
t22 = t138 * t40 + t141 * t41;
t37 = qJD(3) * t149 - qJD(4) * t62;
t36 = qJD(3) * t90 + t95 * t172 - t173 * t94;
t144 = qJD(1) ^ 2;
t126 = -pkin(7) + t136;
t117 = qJDD(2) - t180;
t99 = t143 * t185;
t98 = t140 * t186;
t79 = Ifges(5,4) * t82;
t74 = pkin(4) * t84 + qJD(2);
t73 = -pkin(4) * t90 + t106;
t72 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t83;
t71 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t82;
t48 = t83 * Ifges(5,1) + Ifges(5,5) * qJD(4) + t79;
t47 = t82 * Ifges(5,2) + Ifges(5,6) * qJD(4) + t196;
t46 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t56;
t45 = qJDD(4) * mrSges(5,1) - mrSges(5,3) * t55;
t39 = mrSges(6,1) * t130 - mrSges(6,3) * t52;
t38 = -mrSges(6,2) * t130 + mrSges(6,3) * t162;
t35 = -pkin(4) * t56 + t86;
t31 = -pkin(7) * t85 + t37;
t30 = -pkin(7) * t84 + t36;
t25 = -mrSges(6,1) * t162 + mrSges(6,2) * t52;
t13 = -mrSges(6,2) * t125 + mrSges(6,3) * t18;
t12 = mrSges(6,1) * t125 - mrSges(6,3) * t17;
t11 = t141 * t33 - t184;
t10 = -t138 * t33 - t181;
t5 = -qJD(5) * t22 - t138 * t30 + t141 * t31;
t4 = qJD(5) * t21 + t138 * t31 + t141 * t30;
t1 = [m(3) * (-pkin(1) * t117 + (-t102 + t132) * qJ(2)) + qJ(2) * t177 - (Ifges(4,4) * t135 - Ifges(4,2) * t134) * t169 + t73 * t163 + t106 * t164 - t206 * mrSges(5,3) + (Ifges(5,1) * t85 - Ifges(5,4) * t84) * t198 + (Ifges(4,1) * t135 - Ifges(4,4) * t134) * t168 + t100 * t155 + t97 * (mrSges(5,1) * t84 + mrSges(5,2) * t85) - 0.2e1 * t102 * mrSges(3,3) + m(5) * (t106 * t86 + t19 * t62 + t20 * t61 + t36 * t43 + t37 * t42) - t84 * t47 / 0.2e1 + t82 * (Ifges(5,4) * t85 - Ifges(5,2) * t84) / 0.2e1 + qJD(4) * (Ifges(5,5) * t85 - Ifges(5,6) * t84) / 0.2e1 + t85 * t48 / 0.2e1 + t36 * t71 + t37 * t72 + t74 * t25 + t61 * t45 + t62 * t46 + t4 * t38 + t5 * t39 + t21 * t12 + t22 * t13 - (mrSges(5,2) * t86 + Ifges(5,1) * t55 + Ifges(5,4) * t56 + Ifges(5,5) * qJDD(4)) * t149 + (mrSges(6,2) * t35 - mrSges(6,3) * t3 + Ifges(6,1) * t17 + Ifges(6,4) * t18 + Ifges(6,5) * t125) * t211 + (-t215 * (t143 * pkin(1) + t140 * qJ(2)) + (-m(4) * qJ(3) + m(5) * t136 + m(6) * t126 - t205) * t143 + t204 * t140) * g(2) + ((-m(4) * t137 - m(6) * (-pkin(1) + t126) - m(5) * (-pkin(1) + t136) + m(3) * pkin(1) + t205) * t140 + (-qJ(2) * t215 + t204) * t143) * g(1) + (Ifges(6,1) * t145 - Ifges(6,4) * t26) * t200 + t130 * (Ifges(6,5) * t145 - Ifges(6,6) * t26) / 0.2e1 + t63 * (mrSges(6,1) * t26 + mrSges(6,2) * t145) + t162 * (Ifges(6,4) * t145 - Ifges(6,2) * t26) / 0.2e1 - t26 * t218 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1) + t203 * qJD(2) + m(6) * (t2 * t22 + t21 * t3 + t35 * t73 + t4 * t9 + t5 * t8 + t63 * t74) + (-mrSges(5,1) * t86 + Ifges(5,4) * t55 + Ifges(5,2) * t56 + Ifges(5,6) * qJDD(4)) * t90 + t222 * mrSges(6,3) - (-mrSges(6,1) * t35 + mrSges(6,3) * t2 + Ifges(6,4) * t17 + Ifges(6,2) * t18 + Ifges(6,6) * t125) * t59 + t145 * t24 / 0.2e1 + m(4) * (qJ(2) * t100 - qJD(3) * t156 + t137 * t159) + (-t93 - t223) * t160 + (-t180 + t117) * mrSges(3,2); t211 * t12 + t59 * t13 + t26 * t38 + t145 * t39 - t149 * t45 - t90 * t46 + t84 * t71 + t85 * t72 + (-m(3) * qJ(2) - mrSges(3,3)) * t144 + (mrSges(3,2) - t160) * qJDD(1) + m(3) * t117 + m(5) * t206 + m(6) * (t2 * t59 + t211 * t3 - t222) + m(4) * t159 + (-m(6) * t63 - t203 - t25) * qJD(1) + t210 * t215; -t144 * t160 - t162 * t38 + t52 * t39 - t82 * t71 + t83 * t72 + t163 + t164 + t177 + (g(1) * t143 + g(2) * t140) * t171 + (-t162 * t9 + t52 * t8 + t35) * m(6) + (t42 * t83 - t43 * t82 + t86) * m(5) + (qJD(1) * t156 + t100) * m(4); -t83 * (Ifges(5,1) * t82 - t196) / 0.2e1 + (-t98 + (-t154 + t185) * t140) * g(1) + (-t99 + (t154 + t186) * t143) * g(2) + (t141 * t12 + t138 * t13 - t83 * t25 + (g(3) * t118 + t210 * t119 + t138 * t2 + t141 * t3 - t63 * t83) * m(6) + (-t138 * t39 + t141 * t38 + (-t138 * t8 + t141 * t9) * m(6)) * qJD(5)) * pkin(4) + t47 * t198 - t97 * (mrSges(5,1) * t83 + mrSges(5,2) * t82) - qJD(4) * (Ifges(5,5) * t82 - Ifges(5,6) * t83) / 0.2e1 - t42 * t71 + t43 * t72 + Ifges(5,5) * t55 + Ifges(5,6) * t56 - t11 * t38 - t10 * t39 - t19 * mrSges(5,2) + t20 * mrSges(5,1) - (-Ifges(5,2) * t83 + t48 + t79) * t82 / 0.2e1 + t52 * t218 - m(6) * (t10 * t8 + t11 * t9) + t221 + (t42 * t82 + t43 * t83) * mrSges(5,3) + t220 * g(3) + Ifges(5,3) * qJDD(4); t23 * t200 - t8 * t38 + t9 * t39 - g(1) * (-t140 * t185 + t98) - g(2) * (-t143 * t186 + t99) + g(3) * t152 + t221;];
tau = t1;
