% Calculate vector of inverse dynamics joint torques for
% S5RPPPR1
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
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
% Datum: 2022-01-20 09:13
% Revision: 008671b0a00594318b890887636eaaff83cd5e2f (2021-12-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPPR1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR1_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR1_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR1_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR1_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPPR1_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 09:12:14
% EndTime: 2022-01-20 09:12:23
% DurationCPUTime: 4.02s
% Computational Cost: add. (1742->280), mult. (3851->402), div. (0->0), fcn. (2564->14), ass. (0->137)
t109 = sin(pkin(7));
t93 = pkin(1) * t109 + qJ(3);
t78 = qJD(1) * qJD(3) + qJDD(1) * t93;
t169 = m(5) + m(6);
t147 = m(4) + t169;
t108 = sin(pkin(8));
t111 = cos(pkin(8));
t186 = t108 ^ 2 + t111 ^ 2;
t107 = sin(pkin(9));
t114 = sin(qJ(5));
t110 = cos(pkin(9));
t116 = cos(qJ(5));
t153 = t110 * t116;
t127 = t107 * t114 - t153;
t81 = t107 * t116 + t110 * t114;
t74 = t81 * qJD(5);
t26 = (-qJD(1) * t74 - qJDD(1) * t127) * t108;
t185 = Ifges(6,5) * t26;
t177 = t127 * qJD(5);
t27 = (qJD(1) * t177 - qJDD(1) * t81) * t108;
t184 = Ifges(6,6) * t27;
t144 = qJDD(1) * t111;
t89 = qJDD(5) - t144;
t183 = Ifges(6,3) * t89;
t151 = qJD(1) * t111;
t182 = t127 * t151 - t177;
t122 = qJD(1) * t81;
t181 = t111 * t122 - t74;
t48 = qJDD(2) * t111 - t108 * t78;
t49 = qJDD(2) * t108 + t111 * t78;
t179 = -t108 * t48 + t111 * t49;
t130 = mrSges(5,1) * t107 + mrSges(5,2) * t110;
t123 = t130 * t108;
t152 = qJD(1) * t108;
t138 = t107 * t152;
t85 = t93 * qJD(1);
t59 = qJD(2) * t111 - t108 * t85;
t58 = qJD(4) - t59;
t38 = pkin(4) * t138 + t58;
t52 = t108 * t122;
t54 = -t114 * t138 + t152 * t153;
t175 = m(6) * t38 + mrSges(6,1) * t52 + mrSges(6,2) * t54 + qJD(1) * t123;
t155 = t108 * t110;
t142 = pkin(6) * t155;
t126 = -pkin(4) * t111 - t142;
t148 = qJD(4) * t108;
t134 = -qJ(4) * t108 - pkin(2);
t112 = cos(pkin(7));
t166 = pkin(1) * t112;
t77 = -pkin(3) * t111 + t134 - t166;
t37 = -qJD(1) * t148 + qJDD(1) * t77 + qJDD(3);
t13 = -t107 * t49 + t110 * t37;
t10 = qJDD(1) * t126 + t13;
t145 = qJDD(1) * t108;
t136 = t107 * t145;
t14 = t107 * t37 + t110 * t49;
t11 = -pkin(6) * t136 + t14;
t50 = qJD(1) * t77 + qJD(3);
t60 = qJD(2) * t108 + t111 * t85;
t20 = -t107 * t60 + t110 * t50;
t12 = qJD(1) * t126 + t20;
t21 = t107 * t50 + t110 * t60;
t15 = -pkin(6) * t138 + t21;
t5 = -t114 * t15 + t116 * t12;
t1 = qJD(5) * t5 + t10 * t114 + t11 * t116;
t6 = t114 * t12 + t116 * t15;
t2 = -qJD(5) * t6 + t10 * t116 - t11 * t114;
t174 = -t2 * mrSges(6,1) + t1 * mrSges(6,2);
t164 = pkin(4) * t107;
t173 = -m(6) * t164 + mrSges(3,2) - mrSges(4,3) - t130;
t131 = mrSges(4,1) * t111 - mrSges(4,2) * t108;
t172 = -(m(5) * pkin(3) + t110 * mrSges(5,1) - t107 * mrSges(5,2)) * t111 - mrSges(3,1) - t131 - t108 * mrSges(6,3);
t170 = t54 / 0.2e1;
t135 = t110 * t145;
t61 = mrSges(5,1) * t136 + mrSges(5,2) * t135;
t7 = -t27 * mrSges(6,1) + t26 * mrSges(6,2);
t168 = -t61 - t7;
t167 = Ifges(6,4) * t54;
t115 = sin(qJ(1));
t165 = pkin(1) * t115;
t117 = cos(qJ(1));
t102 = t117 * pkin(1);
t154 = t110 * t111;
t36 = t107 * t77 + t93 * t154;
t106 = qJ(1) + pkin(7);
t99 = sin(t106);
t160 = t111 * t99;
t101 = cos(t106);
t158 = t101 * t111;
t157 = t107 * t108;
t156 = t107 * t111;
t149 = qJD(3) * t111;
t143 = t183 + t184 + t185;
t97 = -pkin(2) - t166;
t132 = -mrSges(4,1) * t144 + mrSges(4,2) * t145;
t129 = -t108 * t59 + t111 * t60;
t65 = t110 * t77;
t28 = -t142 + t65 + (-t107 * t93 - pkin(4)) * t111;
t29 = -pkin(6) * t157 + t36;
t8 = -t114 * t29 + t116 * t28;
t9 = t114 * t28 + t116 * t29;
t128 = -t108 * (-pkin(6) - qJ(4)) + t111 * (pkin(4) * t110 + pkin(3));
t42 = qJDD(4) - t48;
t125 = -mrSges(5,1) * t111 - mrSges(5,3) * t155;
t124 = mrSges(5,2) * t111 - mrSges(5,3) * t157;
t105 = pkin(9) + qJ(5);
t100 = cos(t105);
t98 = sin(t105);
t90 = qJD(5) - t151;
t83 = qJDD(1) * t97 + qJDD(3);
t76 = t125 * qJD(1);
t75 = t124 * qJD(1);
t72 = t125 * qJDD(1);
t71 = t124 * qJDD(1);
t70 = -t107 * t148 + t110 * t149;
t69 = -t107 * t149 - t110 * t148;
t67 = (t93 + t164) * t108;
t63 = t127 * t108;
t62 = t81 * t108;
t57 = t108 * t177;
t56 = t108 * t74;
t47 = Ifges(6,4) * t52;
t46 = t100 * t158 + t98 * t99;
t45 = t100 * t99 - t158 * t98;
t44 = -t100 * t160 + t101 * t98;
t43 = t100 * t101 + t160 * t98;
t35 = -t156 * t93 + t65;
t34 = mrSges(6,1) * t90 - mrSges(6,3) * t54;
t33 = -mrSges(6,2) * t90 - mrSges(6,3) * t52;
t32 = pkin(4) * t136 + t42;
t19 = Ifges(6,1) * t54 + Ifges(6,5) * t90 - t47;
t18 = -Ifges(6,2) * t52 + Ifges(6,6) * t90 + t167;
t17 = -mrSges(6,2) * t89 + mrSges(6,3) * t27;
t16 = mrSges(6,1) * t89 - mrSges(6,3) * t26;
t4 = -qJD(5) * t9 - t114 * t70 + t116 * t69;
t3 = qJD(5) * t8 + t114 * t69 + t116 * t70;
t22 = [((Ifges(5,1) * t110 - Ifges(5,4) * t107) * t135 - (Ifges(5,4) * t110 - Ifges(5,2) * t107) * t136 + m(5) * (qJD(3) * t58 + t42 * t93) + t93 * t61 + Ifges(4,1) * t145 + (-Ifges(5,5) * t110 + Ifges(5,6) * t107 + Ifges(4,4)) * t144) * t108 + (-m(3) * t102 - mrSges(2,1) * t117 - t46 * mrSges(6,1) + mrSges(2,2) * t115 - t45 * mrSges(6,2) - t147 * (t101 * pkin(2) + t99 * qJ(3) + t102) + t173 * t99 + (-(m(5) * qJ(4) + mrSges(5,3)) * t108 - m(6) * t128 + t172) * t101) * g(2) + (m(3) * t165 + mrSges(2,1) * t115 - t44 * mrSges(6,1) + mrSges(2,2) * t117 - t43 * mrSges(6,2) - t147 * (t101 * qJ(3) - t165) + t173 * t101 + (-m(5) * t134 + t108 * mrSges(5,3) + m(4) * pkin(2) - m(6) * (-pkin(2) - t128) - t172) * t99) * g(1) + (t186 * t78 + t179) * mrSges(4,3) + m(4) * (t129 * qJD(3) + t179 * t93 + t83 * t97) + (Ifges(3,3) + Ifges(2,3) + (0.2e1 * t112 * mrSges(3,1) - 0.2e1 * t109 * mrSges(3,2) + m(3) * (t109 ^ 2 + t112 ^ 2) * pkin(1)) * pkin(1)) * qJDD(1) - (mrSges(6,2) * t32 - mrSges(6,3) * t2 + Ifges(6,1) * t26 + Ifges(6,4) * t27 + Ifges(6,5) * t89) * t63 + t97 * t132 - t83 * t131 + t14 * t124 + t13 * t125 - (-mrSges(6,1) * t32 + mrSges(6,3) * t1 + Ifges(6,4) * t26 + Ifges(6,2) * t27 + Ifges(6,6) * t89) * t62 + t67 * t7 + t36 * t71 + t35 * t72 + t70 * t75 + t69 * t76 - t56 * t19 / 0.2e1 + t57 * t18 / 0.2e1 + m(5) * (t13 * t35 + t14 * t36 + t20 * t69 + t21 * t70) + t3 * t33 + t4 * t34 + t8 * t16 + t9 * t17 + m(6) * (t1 * t9 + t2 * t8 + t3 * t6 + t32 * t67 + t4 * t5) + t175 * qJD(3) * t108 + (t5 * t56 + t6 * t57) * mrSges(6,3) + t90 * (-Ifges(6,5) * t56 + Ifges(6,6) * t57) / 0.2e1 + t38 * (-mrSges(6,1) * t57 - mrSges(6,2) * t56) - t52 * (-Ifges(6,4) * t56 + Ifges(6,2) * t57) / 0.2e1 + (-Ifges(6,1) * t56 + Ifges(6,4) * t57) * t170 + (-Ifges(5,5) * t135 - t143 / 0.2e1 + Ifges(5,6) * t136 - t183 / 0.2e1 - t184 / 0.2e1 - t185 / 0.2e1 + Ifges(4,4) * t145 + (Ifges(5,3) + Ifges(4,2)) * t144 + t174) * t111 + t42 * t123; m(3) * qJDD(2) - t62 * t16 - t63 * t17 - t56 * t33 + t57 * t34 + t168 * t111 + (-t107 * t72 + t110 * t71) * t108 + m(4) * (t108 * t49 + t111 * t48) + m(5) * (-t111 * t42 + (-t107 * t13 + t110 * t14) * t108) + m(6) * (-t1 * t63 - t111 * t32 - t2 * t62 + t5 * t57 - t56 * t6) + (-m(3) - t147) * g(3); t107 * t71 + t110 * t72 - t127 * t16 + t81 * t17 + t181 * t34 + t182 * t33 + m(5) * (t14 * t107 + t13 * t110) + m(4) * t83 + t132 + (-g(1) * t99 + g(2) * t101) * t147 + (t1 * t81 - t127 * t2 + t181 * t5 + t182 * t6) * m(6) + ((t107 * t76 - t110 * t75) * t111 - m(5) * (t154 * t21 - t156 * t20) - m(4) * t129 - t186 * mrSges(4,3) * qJD(1) + (-m(5) * t58 - t175) * t108) * qJD(1); t52 * t33 + t54 * t34 + t169 * t111 * g(3) + m(5) * t42 + ((t107 * t75 + t110 * t76 - m(5) * (-t107 * t21 - t110 * t20)) * qJD(1) + t169 * (-g(1) * t101 - g(2) * t99)) * t108 - t168 + (t5 * t54 + t52 * t6 + t32) * m(6); -t38 * (mrSges(6,1) * t54 - mrSges(6,2) * t52) - t54 * (-Ifges(6,1) * t52 - t167) / 0.2e1 + t18 * t170 - t90 * (-Ifges(6,5) * t52 - Ifges(6,6) * t54) / 0.2e1 - t5 * t33 + t6 * t34 - g(1) * (mrSges(6,1) * t45 - mrSges(6,2) * t46) - g(2) * (-mrSges(6,1) * t43 + mrSges(6,2) * t44) - g(3) * (-mrSges(6,1) * t98 - mrSges(6,2) * t100) * t108 + (-t5 * t52 + t54 * t6) * mrSges(6,3) + t143 + (-Ifges(6,2) * t54 + t19 - t47) * t52 / 0.2e1 - t174;];
tau = t22;
