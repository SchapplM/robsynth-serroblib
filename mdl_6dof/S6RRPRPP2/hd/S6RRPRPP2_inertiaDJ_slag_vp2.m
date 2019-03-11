% Calculate time derivative of joint inertia matrix for
% S6RRPRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRPP2_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP2_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP2_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP2_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP2_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPP2_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPP2_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:49:50
% EndTime: 2019-03-09 09:49:59
% DurationCPUTime: 5.24s
% Computational Cost: add. (2780->367), mult. (6336->504), div. (0->0), fcn. (5472->6), ass. (0->160)
t246 = Ifges(7,4) + Ifges(6,5);
t245 = Ifges(7,2) + Ifges(6,3);
t239 = Ifges(5,1) + Ifges(6,1) + Ifges(7,1);
t234 = Ifges(6,4) + Ifges(5,5);
t166 = -Ifges(7,5) + t234;
t121 = cos(qJ(4));
t244 = t246 * t121;
t119 = sin(qJ(4));
t243 = t246 * t119;
t241 = Ifges(6,6) - Ifges(7,6);
t240 = t245 * t119 + t244;
t193 = Ifges(5,4) * t119;
t238 = t239 * t121 - t193 + t243;
t237 = Ifges(6,2) + Ifges(5,3);
t236 = -t245 * t121 + t243;
t168 = qJD(4) * t121;
t118 = sin(pkin(9));
t120 = sin(qJ(2));
t122 = cos(qJ(2));
t180 = cos(pkin(9));
t88 = t118 * t122 + t120 * t180;
t158 = t88 * t168;
t126 = -t118 * t120 + t122 * t180;
t81 = t126 * qJD(2);
t184 = t119 * t81;
t129 = t158 + t184;
t169 = qJD(4) * t119;
t159 = t88 * t169;
t182 = t121 * t81;
t128 = t159 - t182;
t179 = qJ(5) * t121;
t211 = pkin(4) + pkin(5);
t235 = -t119 * t211 + t179;
t218 = -2 * mrSges(7,3);
t80 = t88 * qJD(2);
t233 = -t128 * t246 + t245 * t129 + t241 * t80;
t27 = mrSges(5,1) * t80 + mrSges(5,3) * t128;
t28 = -t80 * mrSges(6,1) - mrSges(6,2) * t128;
t232 = -t27 + t28;
t231 = -t126 * t241 + t240 * t88;
t181 = t121 * t88;
t53 = -mrSges(5,1) * t126 - mrSges(5,3) * t181;
t54 = mrSges(6,1) * t126 + mrSges(6,2) * t181;
t230 = -t53 + t54;
t229 = t80 * qJ(5) - qJD(5) * t126;
t228 = t240 * qJD(4);
t227 = Ifges(6,6) * t169 + t168 * t234;
t196 = -qJ(3) - pkin(7);
t146 = qJD(2) * t196;
t124 = -qJD(3) * t120 + t122 * t146;
t79 = qJD(3) * t122 + t120 * t146;
t42 = t118 * t124 + t180 * t79;
t171 = qJD(2) * t120;
t161 = pkin(2) * t171;
t43 = pkin(3) * t80 - pkin(8) * t81 + t161;
t112 = -pkin(2) * t122 - pkin(1);
t49 = -pkin(3) * t126 - pkin(8) * t88 + t112;
t163 = t119 * t43 + t121 * t42 + t49 * t168;
t101 = t196 * t120;
t102 = t196 * t122;
t60 = t118 * t101 - t102 * t180;
t6 = -t169 * t60 + t163;
t164 = t119 * t42 + t60 * t168 + t49 * t169;
t7 = t121 * t43 - t164;
t225 = -t119 * t7 + t121 * t6;
t3 = t6 + t229;
t4 = -pkin(4) * t80 - t7;
t224 = t119 * t4 + t121 * t3;
t223 = t166 * t80 + (-Ifges(5,4) + t246) * t129 - t239 * t128;
t199 = t126 * Ifges(7,5);
t222 = -t126 * t234 + t238 * t88 + t199;
t221 = t238 * qJD(4);
t192 = Ifges(5,4) * t121;
t220 = t119 * t239 + t192 - t244;
t219 = 2 * m(7);
t217 = -2 * Ifges(4,4);
t41 = t118 * t79 - t180 * t124;
t216 = 0.2e1 * t41;
t59 = -t180 * t101 - t102 * t118;
t215 = 0.2e1 * t59;
t212 = m(6) + m(7);
t207 = pkin(2) * t118;
t204 = t121 * pkin(4);
t201 = t41 * t59;
t200 = t80 * mrSges(4,3);
t198 = -mrSges(6,2) + mrSges(7,3);
t197 = -Ifges(5,6) - Ifges(7,6);
t29 = mrSges(7,2) * t80 + mrSges(7,3) * t129;
t31 = -mrSges(6,2) * t129 + mrSges(6,3) * t80;
t195 = t29 + t31;
t183 = t119 * t88;
t50 = -mrSges(7,2) * t126 + mrSges(7,3) * t183;
t55 = -mrSges(6,2) * t183 - mrSges(6,3) * t126;
t194 = t50 + t55;
t24 = t119 * t49 + t121 * t60;
t187 = Ifges(7,5) * t121;
t186 = qJ(5) * t81;
t185 = t119 * mrSges(7,1);
t178 = qJD(4) * t88;
t175 = t119 * qJ(5);
t110 = pkin(8) + t207;
t173 = qJ(6) - t110;
t170 = qJD(2) * t122;
t167 = t119 * qJD(5);
t17 = -qJ(5) * t126 + t24;
t160 = qJ(5) * t178;
t153 = t180 * pkin(2);
t150 = t80 * mrSges(4,1) + t81 * mrSges(4,2);
t147 = t168 / 0.2e1;
t57 = t119 * t60;
t23 = t121 * t49 - t57;
t145 = 0.2e1 * t161;
t86 = t173 * t121;
t144 = -pkin(4) * t169 + t167;
t111 = -t153 - pkin(3);
t113 = mrSges(7,2) * t168;
t91 = -mrSges(7,1) * t169 + t113;
t141 = -t121 * mrSges(5,1) + t119 * mrSges(5,2);
t140 = t119 * mrSges(5,1) + t121 * mrSges(5,2);
t99 = -t121 * mrSges(6,1) - t119 * mrSges(6,3);
t139 = t119 * mrSges(6,1) - t121 * mrSges(6,3);
t135 = -Ifges(5,2) * t119 + t192;
t132 = -qJ(6) * t81 - qJD(6) * t88;
t26 = -t80 * mrSges(7,1) + mrSges(7,3) * t128;
t131 = t129 * Ifges(6,6) + t182 * t234 + t237 * t80;
t130 = t121 * (m(6) * t110 - t198);
t34 = -Ifges(5,6) * t126 + t135 * t88;
t127 = -t126 * t197 + t231 - t34;
t125 = t111 - t175;
t21 = -mrSges(7,1) * t129 - mrSges(7,2) * t128;
t105 = Ifges(5,2) * t121 + t193;
t100 = mrSges(7,1) * t121 + mrSges(7,2) * t119;
t95 = t135 * qJD(4);
t92 = t140 * qJD(4);
t90 = t139 * qJD(4);
t85 = t173 * t119;
t84 = t125 - t204;
t78 = qJ(5) * t168 + t144;
t70 = t121 * t211 - t125;
t63 = -qJD(4) * t86 - qJD(6) * t119;
t62 = -qJD(6) * t121 + t169 * t173;
t61 = (-pkin(5) * t119 + t179) * qJD(4) + t144;
t52 = mrSges(7,1) * t126 - mrSges(7,3) * t181;
t51 = mrSges(5,2) * t126 - mrSges(5,3) * t183;
t48 = (mrSges(7,2) * t121 - t185) * t88;
t47 = t139 * t88;
t30 = -mrSges(5,2) * t80 - mrSges(5,3) * t129;
t25 = (pkin(4) * t119 - t179) * t88 + t59;
t22 = mrSges(5,1) * t129 - mrSges(5,2) * t128;
t20 = mrSges(6,1) * t129 + mrSges(6,3) * t128;
t19 = t235 * t88 - t59;
t18 = pkin(4) * t126 - t23;
t13 = -Ifges(5,4) * t128 - Ifges(5,2) * t129 + t80 * Ifges(5,6);
t10 = qJ(6) * t183 + t17;
t9 = t57 + t211 * t126 + (-qJ(6) * t88 - t49) * t121;
t8 = (pkin(4) * t81 + t160) * t119 + (-t186 + (pkin(4) * qJD(4) - qJD(5)) * t88) * t121 + t41;
t5 = (-t211 * t81 - t160) * t119 + (t186 + (-qJD(4) * t211 + qJD(5)) * t88) * t121 - t41;
t2 = qJ(6) * t158 + (-qJD(4) * t60 - t132) * t119 + t163 + t229;
t1 = qJ(6) * t159 - t211 * t80 + (t132 - t43) * t121 + t164;
t11 = [0.2e1 * m(6) * (t17 * t3 + t18 * t4 + t25 * t8) + (mrSges(4,3) * t215 - t126 * t217 + (t222 + t199) * t121 + t127 * t119) * t81 + (mrSges(4,2) * t145 + mrSges(4,3) * t216 + 0.2e1 * Ifges(4,1) * t81 + (mrSges(5,2) * t216 + t223) * t121 + (mrSges(5,1) * t216 - t13 + t233) * t119 + (t217 + t166 * t121 + (Ifges(6,6) + t197) * t119) * t80 + (t127 * t121 + (t126 * t166 - t222) * t119) * qJD(4)) * t88 + 0.2e1 * t112 * t150 + 0.2e1 * (-pkin(1) * (mrSges(3,1) * t120 + mrSges(3,2) * t122) + (-Ifges(3,2) + Ifges(3,1)) * t120 * t122 + (-t120 ^ 2 + t122 ^ 2) * Ifges(3,4)) * qJD(2) - (mrSges(4,1) * t145 - 0.2e1 * t42 * mrSges(4,3) + ((2 * Ifges(4,2)) + (2 * Ifges(7,3)) + t237) * t80 + t131) * t126 + t22 * t215 + (t1 * t9 + t10 * t2 + t19 * t5) * t219 - 0.2e1 * t60 * t200 + 0.2e1 * m(5) * (t23 * t7 + t24 * t6 + t201) + 0.2e1 * m(4) * (t112 * t161 + t42 * t60 + t201) + 0.2e1 * t19 * t21 + 0.2e1 * t25 * t20 + 0.2e1 * t9 * t26 + 0.2e1 * t23 * t27 + 0.2e1 * t18 * t28 + 0.2e1 * t10 * t29 + 0.2e1 * t24 * t30 + 0.2e1 * t17 * t31 + 0.2e1 * t8 * t47 + 0.2e1 * t5 * t48 + 0.2e1 * t2 * t50 + 0.2e1 * t6 * t51 + 0.2e1 * t1 * t52 + 0.2e1 * t7 * t53 + 0.2e1 * t4 * t54 + 0.2e1 * t3 * t55; (-mrSges(3,1) * t170 + mrSges(3,2) * t171) * pkin(7) + t41 * t141 + (-t1 * t119 + t10 * t169 - t121 * t2 - t168 * t9) * mrSges(7,3) - (t220 * t88 + t34) * t169 / 0.2e1 - (-Ifges(5,6) * t169 + t227) * t126 / 0.2e1 + m(4) * (t118 * t42 - t180 * t41) * pkin(2) + t126 * (Ifges(7,6) * t119 + t187) * qJD(4) / 0.2e1 - t200 * t207 + m(7) * (-t1 * t85 + t10 * t62 + t19 * t61 - t2 * t86 + t5 * t70 + t63 * t9) - Ifges(3,6) * t171 + ((t31 + t30) * t121 + (-t51 - t55) * t169 + t230 * t168 + t232 * t119 + m(6) * ((-t119 * t17 + t121 * t18) * qJD(4) + t224) + m(5) * ((-t119 * t24 - t121 * t23) * qJD(4) + t225)) * t110 + Ifges(3,5) * t170 + t236 * (t88 * t147 + t184 / 0.2e1) + (m(5) * t41 + t22) * t111 + m(6) * (-t25 * t78 + t8 * t84) - t41 * mrSges(4,1) - t42 * mrSges(4,2) + t61 * t48 + t62 * t50 + t63 * t52 + t70 * t21 - t78 * t47 - Ifges(4,6) * t80 + t84 * t20 - t85 * t26 - t86 * t29 + t25 * t90 + t19 * t91 + t59 * t92 + t8 * t99 + t5 * t100 + t121 * t13 / 0.2e1 - t80 * (Ifges(7,5) * t119 - Ifges(7,6) * t121) / 0.2e1 + (-t95 / 0.2e1 + t228 / 0.2e1) * t183 + t220 * t182 / 0.2e1 + t221 * t181 / 0.2e1 + t222 * t147 + t223 * t119 / 0.2e1 + (t168 * t18 - t169 * t17 + t224) * mrSges(6,2) + (-t168 * t23 - t169 * t24 + t225) * mrSges(5,3) - t129 * t105 / 0.2e1 + t231 * t169 / 0.2e1 - t233 * t121 / 0.2e1 + ((Ifges(5,6) - Ifges(6,6)) * t121 + t234 * t119) * t80 / 0.2e1 + (-mrSges(4,3) * t153 + Ifges(4,5)) * t81; (t61 * t70 - t62 * t86 - t63 * t85) * t219 + 0.2e1 * t111 * t92 + 0.2e1 * t84 * t90 + 0.2e1 * t61 * t100 + 0.2e1 * t70 * t91 + 0.2e1 * (-m(6) * t84 - t99) * t78 + (t218 * t62 - t228 + t95) * t121 + (t218 * t63 + t221) * t119 + ((-t218 * t85 + t220) * t121 + (t218 * t86 - t105 + t236) * t119) * qJD(4); m(4) * t161 + (-t26 - t232) * t121 + (t30 + t195) * t119 + ((t51 + t194) * t121 + (t52 + t230) * t119) * qJD(4) + m(6) * (t119 * t3 - t121 * t4 + (t119 * t18 + t121 * t17) * qJD(4)) + m(7) * (-t1 * t121 + t119 * t2 + (t10 * t121 + t119 * t9) * qJD(4)) + m(5) * (t119 * t6 + t121 * t7 + (-t119 * t23 + t121 * t24) * qJD(4)) + t150; m(7) * (t119 * t62 - t121 * t63 + (-t119 * t85 - t121 * t86) * qJD(4)); 0; (t119 * t197 - t187) * t81 + t194 * qJD(5) + t195 * qJ(5) + m(7) * (qJ(5) * t2 + qJD(5) * t10 - t1 * t211) + m(6) * (-pkin(4) * t4 + qJ(5) * t3 + qJD(5) * t17) + t131 + (-t119 * t166 + t121 * t197) * t178 - t1 * mrSges(7,1) + t2 * mrSges(7,2) + t3 * mrSges(6,3) - t4 * mrSges(6,1) - t6 * mrSges(5,2) + t7 * mrSges(5,1) - pkin(4) * t28 + Ifges(7,3) * t80 - t211 * t26; m(7) * (qJ(5) * t62 - qJD(5) * t86 - t211 * t63) - t63 * mrSges(7,1) + t62 * mrSges(7,2) + qJD(5) * t130 + ((-mrSges(6,2) * pkin(4) + mrSges(7,3) * t211 - Ifges(7,5)) * t121 + (qJ(5) * t198 + t197) * t119 + (m(6) * (-t175 - t204) + t99 + t141) * t110) * qJD(4) + t227; t113 + m(7) * t167 + m(6) * t78 + (m(7) * t235 - t139 - t140 - t185) * qJD(4); 0.2e1 * (qJ(5) * t212 + mrSges(7,2) + mrSges(6,3)) * qJD(5); m(6) * t4 + m(7) * t1 + t26 + t28; m(7) * t63 + qJD(4) * t130; t212 * t169; 0; 0; m(7) * t5 + t21; m(7) * t61 + t91; 0; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t11(1) t11(2) t11(4) t11(7) t11(11) t11(16); t11(2) t11(3) t11(5) t11(8) t11(12) t11(17); t11(4) t11(5) t11(6) t11(9) t11(13) t11(18); t11(7) t11(8) t11(9) t11(10) t11(14) t11(19); t11(11) t11(12) t11(13) t11(14) t11(15) t11(20); t11(16) t11(17) t11(18) t11(19) t11(20) t11(21);];
Mq  = res;
