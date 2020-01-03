% Calculate vector of inverse dynamics joint torques for
% S5RPPRR7
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
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
% Datum: 2019-12-31 18:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPRR7_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR7_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR7_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR7_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR7_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR7_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR7_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR7_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR7_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:59:31
% EndTime: 2019-12-31 17:59:38
% DurationCPUTime: 4.01s
% Computational Cost: add. (1868->331), mult. (3473->466), div. (0->0), fcn. (1826->10), ass. (0->156)
t82 = sin(qJ(4));
t204 = -t82 / 0.2e1;
t84 = cos(qJ(5));
t147 = qJD(4) * t84;
t85 = cos(qJ(4));
t151 = qJD(1) * t85;
t81 = sin(qJ(5));
t55 = -t151 * t81 + t147;
t203 = -t55 / 0.2e1;
t149 = qJD(4) * t81;
t56 = t151 * t84 + t149;
t202 = -t56 / 0.2e1;
t152 = qJD(1) * t82;
t67 = qJD(5) + t152;
t201 = -t67 / 0.2e1;
t169 = Ifges(5,4) * t85;
t200 = Ifges(5,2) * t204 + t169 / 0.2e1;
t78 = qJ(1) + pkin(8);
t74 = sin(t78);
t75 = cos(t78);
t193 = -g(1) * t74 + g(2) * t75;
t178 = t85 / 0.2e1;
t199 = -m(5) - m(6);
t140 = qJD(1) * qJD(4);
t59 = qJDD(1) * t85 - t140 * t82;
t24 = qJD(5) * t55 + qJDD(4) * t81 + t59 * t84;
t25 = -qJD(5) * t56 + qJDD(4) * t84 - t59 * t81;
t5 = -mrSges(6,1) * t25 + mrSges(6,2) * t24;
t171 = -qJDD(4) * mrSges(5,1) + mrSges(5,3) * t59 + t5;
t134 = mrSges(5,3) * t151;
t198 = -qJD(4) * mrSges(5,1) - mrSges(6,1) * t55 + mrSges(6,2) * t56 + t134;
t170 = Ifges(5,4) * t82;
t112 = t85 * Ifges(5,1) - t170;
t50 = Ifges(6,4) * t55;
t21 = Ifges(6,1) * t56 + Ifges(6,5) * t67 + t50;
t197 = Ifges(5,5) * qJD(4) + qJD(1) * t112 + t84 * t21;
t80 = cos(pkin(8));
t73 = -pkin(1) * t80 - pkin(2);
t66 = -pkin(6) + t73;
t51 = qJD(1) * t66 + qJD(3);
t34 = qJD(2) * t85 + t51 * t82;
t150 = qJD(4) * t34;
t48 = qJDD(1) * t66 + qJDD(3);
t15 = -qJDD(2) * t82 + t48 * t85 - t150;
t144 = qJD(5) * t84;
t145 = qJD(5) * t81;
t30 = qJD(4) * pkin(7) + t34;
t120 = pkin(4) * t82 - pkin(7) * t85;
t79 = sin(pkin(8));
t70 = pkin(1) * t79 + qJ(3);
t47 = t120 + t70;
t35 = t47 * qJD(1);
t8 = -t30 * t81 + t35 * t84;
t9 = t30 * t84 + t35 * t81;
t196 = -t8 * t144 - t9 * t145;
t60 = -qJDD(1) * t82 - t140 * t85;
t53 = qJDD(5) - t60;
t10 = mrSges(6,1) * t53 - mrSges(6,3) * t24;
t11 = -mrSges(6,2) * t53 + mrSges(6,3) * t25;
t195 = -t81 * t10 + t84 * t11;
t146 = qJD(4) * t85;
t148 = qJD(4) * t82;
t14 = -qJD(2) * t148 + t85 * qJDD(2) + t51 * t146 + t82 * t48;
t194 = t14 * t82 + t15 * t85;
t192 = m(4) - t199;
t191 = -mrSges(3,1) + mrSges(4,2) - mrSges(5,3);
t190 = Ifges(5,6) * qJD(4) / 0.2e1 + qJD(1) * t200 + Ifges(6,5) * t202 + Ifges(6,6) * t203 + Ifges(6,3) * t201;
t113 = mrSges(6,1) * t81 + mrSges(6,2) * t84;
t179 = -t81 / 0.2e1;
t168 = Ifges(6,4) * t56;
t20 = Ifges(6,2) * t55 + Ifges(6,6) * t67 + t168;
t33 = -qJD(2) * t82 + t51 * t85;
t29 = -qJD(4) * pkin(4) - t33;
t189 = t113 * t29 + t20 * t179;
t115 = mrSges(5,1) * t82 + mrSges(5,2) * t85;
t156 = t85 * mrSges(6,3);
t188 = -m(6) * t120 + mrSges(3,2) - mrSges(4,3) - t115 + t156;
t13 = -qJDD(4) * pkin(4) - t15;
t31 = -mrSges(6,2) * t67 + mrSges(6,3) * t55;
t32 = mrSges(6,1) * t67 - mrSges(6,3) * t56;
t135 = mrSges(5,3) * t152;
t63 = -qJD(4) * mrSges(5,2) - t135;
t186 = m(5) * (t15 + t150) - m(6) * (-t147 * t9 + t149 * t8 + t13) + qJD(4) * (t31 * t84 - t32 * t81 + t63) - t171;
t184 = t24 / 0.2e1;
t183 = t25 / 0.2e1;
t182 = t53 / 0.2e1;
t180 = t56 / 0.2e1;
t177 = -m(5) - m(4);
t83 = sin(qJ(1));
t175 = pkin(1) * t83;
t12 = qJDD(4) * pkin(7) + t14;
t141 = qJD(1) * qJD(3);
t52 = qJDD(1) * t70 + t141;
t18 = -pkin(4) * t60 - pkin(7) * t59 + t52;
t2 = -qJD(5) * t9 - t12 * t81 + t18 * t84;
t172 = t2 * t81;
t86 = cos(qJ(1));
t77 = t86 * pkin(1);
t167 = Ifges(6,4) * t81;
t166 = Ifges(6,4) * t84;
t161 = t81 * t82;
t160 = t82 * t84;
t159 = t84 * mrSges(6,3);
t143 = qJD(5) * t85;
t142 = qJDD(1) * mrSges(4,2);
t139 = Ifges(6,5) * t24 + Ifges(6,6) * t25 + Ifges(6,3) * t53;
t138 = t75 * pkin(2) + t74 * qJ(3) + t77;
t132 = t66 * t146;
t131 = t81 * t148;
t129 = m(6) * pkin(7) + mrSges(6,3);
t124 = -t140 / 0.2e1;
t121 = pkin(4) * t85 + pkin(7) * t82;
t1 = qJD(5) * t8 + t12 * t84 + t18 * t81;
t118 = t1 * t84 - t172;
t117 = t8 * t84 + t81 * t9;
t116 = mrSges(5,1) * t85 - mrSges(5,2) * t82;
t114 = -mrSges(6,1) * t84 + mrSges(6,2) * t81;
t111 = Ifges(6,1) * t84 - t167;
t110 = Ifges(6,1) * t81 + t166;
t108 = -Ifges(6,2) * t81 + t166;
t107 = Ifges(6,2) * t84 + t167;
t106 = -Ifges(5,5) * t82 - Ifges(5,6) * t85;
t105 = Ifges(6,5) * t84 - Ifges(6,6) * t81;
t104 = Ifges(6,5) * t81 + Ifges(6,6) * t84;
t103 = -t81 * t31 - t84 * t32;
t62 = t70 * qJD(1);
t102 = qJD(3) * t62 + t52 * t70;
t27 = t160 * t66 + t47 * t81;
t26 = -t161 * t66 + t47 * t84;
t99 = t62 * t116;
t98 = t82 * (-Ifges(5,2) * t85 - t170);
t97 = t85 * (-Ifges(5,1) * t82 - t169);
t96 = m(6) * pkin(4) - t114;
t95 = t143 * t81 + t147 * t82;
t94 = -t143 * t84 + t131;
t91 = Ifges(6,5) * t85 - t111 * t82;
t90 = Ifges(6,6) * t85 - t108 * t82;
t89 = Ifges(6,3) * t85 - t105 * t82;
t43 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t60;
t88 = t43 + t103 * qJD(5) + t198 * qJD(4) + m(5) * (-qJD(4) * t33 + t14) + m(6) * (qJD(4) * t29 + t118 + t196) + t195;
t61 = qJDD(1) * t73 + qJDD(3);
t58 = t121 * qJD(1);
t57 = t115 * qJD(1);
t54 = qJD(4) * t121 + qJD(3);
t49 = t113 * t85;
t39 = t160 * t75 - t74 * t81;
t38 = t161 * t75 + t74 * t84;
t37 = t160 * t74 + t75 * t81;
t36 = -t161 * t74 + t75 * t84;
t17 = t33 * t84 + t58 * t81;
t16 = -t33 * t81 + t58 * t84;
t7 = -qJD(5) * t27 - t132 * t81 + t54 * t84;
t6 = qJD(5) * t26 + t132 * t84 + t54 * t81;
t4 = t24 * Ifges(6,1) + t25 * Ifges(6,4) + t53 * Ifges(6,5);
t3 = t24 * Ifges(6,4) + t25 * Ifges(6,2) + t53 * Ifges(6,6);
t19 = [t20 * t131 / 0.2e1 + t2 * (mrSges(6,1) * t82 - t156 * t84) + t1 * (-mrSges(6,2) * t82 - t156 * t81) + (t52 + t141) * mrSges(4,3) + t55 * (qJD(4) * t90 - t107 * t143) / 0.2e1 + t67 * (qJD(4) * t89 - t104 * t143) / 0.2e1 + t97 * t140 / 0.2e1 + (m(6) * (-t13 * t85 + t148 * t29) + t198 * t148 - t171 * t85 + m(5) * ((-t33 * t82 + t34 * t85) * qJD(4) + t194) + t82 * t43) * t66 + t29 * (-mrSges(6,1) * t94 - mrSges(6,2) * t95) + t85 * t3 * t179 + t60 * t200 + m(5) * t102 + (-t146 * t34 + t148 * t33 - t194) * mrSges(5,3) + (t8 * mrSges(6,1) - t9 * mrSges(6,2) - t190) * t146 - (t84 * t20 + t81 * t21) * t143 / 0.2e1 - t197 * t148 / 0.2e1 + t70 * (-mrSges(5,1) * t60 + mrSges(5,2) * t59) + t61 * mrSges(4,2) + qJD(3) * t57 + t13 * t49 + t6 * t31 + t7 * t32 + m(6) * (t1 * t27 + t2 * t26 + t6 * t9 + t7 * t8) + t26 * t10 + t27 * t11 + t59 * t112 / 0.2e1 + qJD(4) ^ 2 * t106 / 0.2e1 + t63 * t132 + (m(3) * t175 + mrSges(2,1) * t83 - t39 * mrSges(6,1) + mrSges(2,2) * t86 + t38 * mrSges(6,2) - t192 * (t75 * qJ(3) - t175) + (m(4) * pkin(2) + t199 * (-pkin(2) - pkin(6)) - t191) * t74 + t188 * t75) * g(1) + (-m(3) * t77 - m(4) * t138 - mrSges(2,1) * t86 - t37 * mrSges(6,1) + mrSges(2,2) * t83 - t36 * mrSges(6,2) + t199 * (t75 * pkin(6) + t138) + t191 * t75 + t188 * t74) * g(2) + (0.2e1 * Ifges(5,5) * t178 - Ifges(5,6) * t82) * qJDD(4) + t82 * t139 / 0.2e1 + (Ifges(5,1) * t59 + Ifges(5,4) * t60 + t84 * t4) * t178 + (qJD(4) * t91 - t110 * t143) * t180 + (Ifges(6,3) * t82 + t105 * t85) * t182 + (Ifges(6,6) * t82 + t108 * t85) * t183 + t73 * t142 + qJD(4) * t99 + (Ifges(5,4) * t59 + Ifges(5,2) * t60) * t204 + (t8 * t95 + t9 * t94) * mrSges(6,3) + m(4) * (t61 * t73 + t102) + (Ifges(6,5) * t82 + t111 * t85) * t184 + t98 * t124 + t52 * t115 + (t70 * mrSges(4,3) + Ifges(4,1) + Ifges(2,3) + Ifges(3,3) + (0.2e1 * mrSges(3,1) * t80 - 0.2e1 * mrSges(3,2) * t79 + m(3) * (t79 ^ 2 + t80 ^ 2) * pkin(1)) * pkin(1)) * qJDD(1); (m(3) + m(4)) * qJDD(2) + (-m(3) - m(6) + t177) * g(3) - t186 * t82 + t88 * t85; t142 + t61 * m(4) + t186 * t85 + t88 * t82 + (-m(6) * t117 - mrSges(4,3) * qJD(1) + t177 * t62 + t103 - t57) * qJD(1) + t192 * t193; (-t135 - t63) * t33 + t21 * t144 / 0.2e1 + (-m(6) * t29 + t134 - t198) * t34 + (-t32 * t144 - t31 * t145 + m(6) * (-qJD(5) * t117 + t118) + t195) * pkin(7) + (-t172 + t196) * mrSges(6,3) + t189 * qJD(5) + t190 * t151 + (t98 / 0.2e1 - t97 / 0.2e1) * qJD(1) ^ 2 + (t105 * t67 + t108 * t55 + t111 * t56) * qJD(5) / 0.2e1 - (t55 * t90 + t56 * t91 + t67 * t89) * qJD(1) / 0.2e1 + (-t99 - t8 * (mrSges(6,1) * t85 + t159 * t82) - t9 * (-mrSges(6,2) * t85 + mrSges(6,3) * t161)) * qJD(1) + t84 * t3 / 0.2e1 + t81 * t4 / 0.2e1 + Ifges(5,5) * t59 + Ifges(5,6) * t60 - t17 * t31 - t16 * t32 - t14 * mrSges(5,2) + t15 * mrSges(5,1) - pkin(4) * t5 + t13 * t114 + (-pkin(4) * t13 - t16 * t8 - t17 * t9) * m(6) + t104 * t182 + t1 * t159 + (t129 * t82 + t85 * t96 + t116) * t193 + t107 * t183 + t110 * t184 + t106 * t124 + (-t129 * t85 + t82 * t96 + t115) * g(3) + (t189 + t197 / 0.2e1) * t152 + Ifges(5,3) * qJDD(4); -t1 * mrSges(6,2) + t2 * mrSges(6,1) - t29 * (mrSges(6,1) * t56 + mrSges(6,2) * t55) + (Ifges(6,1) * t55 - t168) * t202 + t20 * t180 + (Ifges(6,5) * t55 - Ifges(6,6) * t56) * t201 - t8 * t31 + t9 * t32 - g(1) * (mrSges(6,1) * t36 - mrSges(6,2) * t37) - g(2) * (mrSges(6,1) * t38 + mrSges(6,2) * t39) + g(3) * t49 + (t55 * t8 + t56 * t9) * mrSges(6,3) + t139 + (-Ifges(6,2) * t56 + t21 + t50) * t203;];
tau = t19;
