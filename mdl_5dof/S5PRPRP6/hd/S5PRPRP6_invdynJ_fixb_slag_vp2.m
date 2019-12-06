% Calculate vector of inverse dynamics joint torques for
% S5PRPRP6
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
%   pkin=[a2,a3,a4,a5,d2,d4,theta1]';
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
% Datum: 2019-12-05 15:41
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRPRP6_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP6_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP6_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRP6_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP6_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRPRP6_invdynJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP6_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRP6_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRP6_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:40:21
% EndTime: 2019-12-05 15:40:31
% DurationCPUTime: 4.00s
% Computational Cost: add. (881->256), mult. (1718->310), div. (0->0), fcn. (821->6), ass. (0->116)
t135 = mrSges(3,1) - mrSges(4,2);
t186 = -mrSges(5,3) - t135;
t185 = Ifges(5,1) + Ifges(6,1);
t57 = sin(qJ(4));
t59 = cos(qJ(4));
t76 = pkin(4) * t57 - qJ(5) * t59;
t83 = t57 * mrSges(6,1) - t59 * mrSges(6,3);
t85 = mrSges(5,1) * t57 + mrSges(5,2) * t59;
t184 = -m(6) * t76 + mrSges(3,2) - t83 - t85;
t116 = t59 * qJD(2);
t171 = mrSges(5,1) + mrSges(6,1);
t99 = mrSges(6,2) * t116;
t129 = -mrSges(5,3) * t116 + qJD(4) * t171 - t99;
t117 = t57 * qJD(2);
t105 = mrSges(6,2) * t117;
t43 = qJD(4) * mrSges(6,3) - t105;
t130 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t117 + t43;
t111 = qJD(2) * qJD(4);
t34 = qJDD(2) * t59 - t111 * t57;
t16 = -qJDD(4) * mrSges(6,1) + t34 * mrSges(6,2);
t35 = qJDD(2) * t57 + t111 * t59;
t18 = -mrSges(6,2) * t35 + qJDD(4) * mrSges(6,3);
t175 = (-qJDD(4) * mrSges(5,2) - mrSges(5,3) * t35 + t18) * t57 + (qJDD(4) * mrSges(5,1) - mrSges(5,3) * t34 - t16) * t59 - (t129 * t57 - t130 * t59) * qJD(4);
t61 = -pkin(2) - pkin(6);
t60 = cos(qJ(2));
t124 = qJD(1) * t60;
t92 = qJD(3) - t124;
t30 = qJD(2) * t61 + t92;
t142 = t30 * t57;
t11 = qJD(4) * qJ(5) + t142;
t141 = t30 * t59;
t9 = -qJD(4) * pkin(4) + qJD(5) - t141;
t88 = t11 * t59 + t9 * t57;
t112 = qJD(1) * qJD(2);
t58 = sin(qJ(2));
t48 = t58 * t112;
t36 = qJDD(1) * t60 - t48;
t75 = qJDD(3) - t36;
t12 = qJDD(2) * t61 + t75;
t122 = qJD(4) * t57;
t4 = t12 * t59 - t122 * t30;
t121 = qJD(4) * t59;
t5 = t57 * t12 + t30 * t121;
t89 = t4 * t59 + t5 * t57;
t2 = qJDD(4) * qJ(5) + qJD(4) * qJD(5) + t5;
t3 = -qJDD(4) * pkin(4) + qJDD(5) - t4;
t90 = t2 * t57 - t3 * t59;
t183 = m(5) * t89 + m(6) * (qJD(4) * t88 + t90) + t175;
t172 = -m(5) - m(6);
t169 = Ifges(6,4) + Ifges(5,5);
t168 = -Ifges(5,6) + Ifges(6,6);
t144 = Ifges(6,5) * t57;
t146 = Ifges(5,4) * t57;
t181 = t185 * t59 + t144 - t146;
t55 = sin(pkin(7));
t56 = cos(pkin(7));
t162 = g(1) * t56 + g(2) * t55;
t119 = t11 * qJD(4);
t178 = -t59 * t119 - t9 * t122 - t90;
t173 = t129 * t59 + t130 * t57;
t151 = g(3) * t60;
t170 = mrSges(5,2) - mrSges(6,3);
t167 = qJD(2) * t181 + qJD(4) * t169;
t32 = t83 * qJD(2);
t33 = t85 * qJD(2);
t131 = t32 + t33;
t145 = Ifges(5,4) * t59;
t166 = t57 * (Ifges(6,3) * t59 - t144) + t59 * (-Ifges(5,1) * t57 - t145);
t125 = qJD(1) * t58;
t38 = qJ(3) + t76;
t10 = qJD(2) * t38 + t125;
t45 = qJD(2) * qJ(3) + t125;
t84 = t59 * mrSges(6,1) + t57 * mrSges(6,3);
t86 = mrSges(5,1) * t59 - mrSges(5,2) * t57;
t165 = t10 * t84 + t45 * t86;
t164 = t168 * t59 - t169 * t57;
t160 = (t57 ^ 2 + t59 ^ 2) * t30;
t62 = qJD(2) ^ 2;
t155 = t59 / 0.2e1;
t154 = -m(4) - m(5);
t143 = Ifges(6,5) * t59;
t140 = t45 * t60;
t139 = t57 * t58;
t137 = t58 * t59;
t134 = mrSges(4,3) - mrSges(3,2);
t128 = t60 * pkin(2) + t58 * qJ(3);
t118 = t45 * qJD(2);
t115 = qJDD(1) * t58;
t114 = qJDD(2) * mrSges(4,2);
t113 = m(6) - t154;
t110 = qJDD(2) * qJ(3);
t98 = t60 * t112;
t94 = -t111 / 0.2e1;
t91 = -qJD(5) * t59 + qJD(3);
t80 = -t57 * Ifges(5,2) + t145;
t77 = pkin(4) * t59 + qJ(5) * t57;
t73 = t110 + t115;
t13 = (qJD(3) + t124) * qJD(2) + t73;
t74 = qJ(3) * t13 + qJD(3) * t45;
t70 = t57 * (-Ifges(5,2) * t59 - t146);
t49 = Ifges(6,5) * t116;
t39 = -qJD(2) * pkin(2) + t92;
t37 = t98 + t115;
t31 = t77 * qJD(2);
t26 = Ifges(5,6) * qJD(4) + qJD(2) * t80;
t25 = Ifges(6,6) * qJD(4) + Ifges(6,3) * t117 + t49;
t24 = -t139 * t55 + t56 * t59;
t23 = t137 * t55 + t56 * t57;
t22 = t139 * t56 + t55 * t59;
t21 = -t137 * t56 + t55 * t57;
t20 = -qJDD(2) * pkin(2) + t75;
t14 = qJD(4) * t77 + t91;
t7 = mrSges(5,1) * t35 + mrSges(5,2) * t34;
t6 = mrSges(6,1) * t35 - mrSges(6,3) * t34;
t1 = pkin(4) * t35 - qJ(5) * t34 + (t91 + t124) * qJD(2) + t73;
t8 = [m(2) * qJDD(1) + (-m(2) - m(3) - t113) * g(3) + (t6 + t7 - t135 * t62 + t134 * qJDD(2) + t173 * qJD(2) + m(3) * t37 + m(4) * (qJD(2) * t39 + t13) + m(6) * (t11 * t117 - t116 * t9 + t1) + m(5) * (qJD(2) * t160 + t13)) * t58 + (t134 * t62 + t135 * qJDD(2) + t131 * qJD(2) + m(3) * t36 + m(4) * (-t20 + t118) + m(6) * (qJD(2) * t10 + t178) + m(5) * (-t89 + t118) - t175) * t60; t181 * t34 / 0.2e1 + ((m(4) * pkin(2) + t172 * t61 + mrSges(6,2) - t186) * t58 + (qJ(3) * (-m(4) + t172) - mrSges(4,3) + t184) * t60) * t162 + (-m(4) * t128 + t172 * (t60 * pkin(6) + t128) + t186 * t60 + t184 * t58) * g(3) + (t169 * qJDD(4) + t185 * t34) * t155 + (t168 * t57 + t169 * t59) * qJDD(4) / 0.2e1 + (t165 + t164 * qJD(4) / 0.2e1) * qJD(4) + (-t26 / 0.2e1 + t25 / 0.2e1) * t121 + t57 * (Ifges(6,5) * t34 + Ifges(6,6) * qJDD(4)) / 0.2e1 + (-t151 + t178) * mrSges(6,2) + (-t48 + t20) * mrSges(4,2) - t57 * (Ifges(5,4) * t34 + Ifges(5,6) * qJDD(4)) / 0.2e1 - pkin(2) * t114 + (t98 - t37) * mrSges(3,2) - t173 * t125 + (t59 * (-Ifges(6,1) * t57 + t143) + t166) * t111 / 0.2e1 - t131 * t124 - t167 * t122 / 0.2e1 + t1 * t83 + t13 * t85 + m(4) * (-pkin(2) * t20 + t74) + (t48 + t36) * mrSges(3,1) - t89 * mrSges(5,3) + (-m(6) * (t10 * t60 + (t11 * t57 - t59 * t9) * t58) - m(5) * (t160 * t58 + t140) - m(4) * (t39 * t58 + t140)) * qJD(1) + m(6) * (t1 * t38 + t10 * t14) + (t143 / 0.2e1 - t80 / 0.2e1 + (Ifges(6,3) + Ifges(5,2) / 0.2e1) * t57 + (Ifges(6,5) - Ifges(5,4)) * t155) * t35 + t183 * t61 + (Ifges(4,1) + Ifges(3,3)) * qJDD(2) + (-g(3) * t58 + qJD(2) * qJD(3) + t110 + t13 - t98) * mrSges(4,3) + t38 * t6 + t14 * t32 + qJD(3) * t33 + qJ(3) * t7 + m(5) * t74 + t70 * t94; t114 - t62 * mrSges(4,3) + m(4) * t20 + (-m(6) * t10 + t154 * t45 - t131) * qJD(2) + (-t162 * t58 + t151) * t113 + t183; (t170 * t22 + t171 * t21) * g(1) + (-t170 * t24 - t171 * t23) * g(2) + (Ifges(6,2) + Ifges(5,3)) * qJDD(4) + t26 * t116 / 0.2e1 + t164 * t94 - t165 * qJD(2) + t129 * t142 - t130 * t141 + t167 * t117 / 0.2e1 + t168 * t35 + t169 * t34 + t9 * t105 + t11 * t99 + qJD(5) * t43 - t31 * t32 + (-pkin(4) * t3 + qJ(5) * t2 + qJD(5) * t11 + t77 * t151 - t10 * t31 - t88 * t30 - g(1) * (-pkin(4) * t21 + qJ(5) * t22) - g(2) * (pkin(4) * t23 - qJ(5) * t24)) * m(6) - pkin(4) * t16 + qJ(5) * t18 - t5 * mrSges(5,2) + t2 * mrSges(6,3) - t3 * mrSges(6,1) + t4 * mrSges(5,1) + (t70 / 0.2e1 - t166 / 0.2e1) * t62 + (t84 + t86) * t151 - (-Ifges(6,1) * t117 + t25 + t49) * t116 / 0.2e1; t32 * t116 - qJD(4) * t43 + (-g(1) * t21 + g(2) * t23 + t10 * t116 - t151 * t59 - t119 + t3) * m(6) + t16;];
tau = t8;
