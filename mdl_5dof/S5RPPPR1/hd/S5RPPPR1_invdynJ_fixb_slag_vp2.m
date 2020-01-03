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
% Datum: 2020-01-03 11:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:20:05
% EndTime: 2020-01-03 11:20:17
% DurationCPUTime: 3.73s
% Computational Cost: add. (1742->279), mult. (3851->398), div. (0->0), fcn. (2564->14), ass. (0->136)
t110 = sin(pkin(7));
t92 = pkin(1) * t110 + qJ(3);
t78 = qJD(1) * qJD(3) + qJDD(1) * t92;
t109 = sin(pkin(8));
t112 = cos(pkin(8));
t184 = t109 ^ 2 + t112 ^ 2;
t167 = m(5) + m(6);
t183 = -m(6) - m(4);
t108 = sin(pkin(9));
t115 = sin(qJ(5));
t111 = cos(pkin(9));
t117 = cos(qJ(5));
t153 = t111 * t117;
t129 = t108 * t115 - t153;
t81 = t108 * t117 + t111 * t115;
t74 = t81 * qJD(5);
t26 = (-qJD(1) * t74 - qJDD(1) * t129) * t109;
t182 = Ifges(6,5) * t26;
t174 = t129 * qJD(5);
t27 = (qJD(1) * t174 - qJDD(1) * t81) * t109;
t181 = Ifges(6,6) * t27;
t144 = qJDD(1) * t112;
t89 = qJDD(5) - t144;
t180 = Ifges(6,3) * t89;
t151 = qJD(1) * t112;
t179 = t129 * t151 - t174;
t124 = qJD(1) * t81;
t178 = t112 * t124 - t74;
t48 = qJDD(2) * t112 - t109 * t78;
t49 = qJDD(2) * t109 + t112 * t78;
t176 = -t109 * t48 + t112 * t49;
t132 = mrSges(4,1) * t112 - mrSges(4,2) * t109;
t173 = -mrSges(3,1) - t132 + (-m(6) * (pkin(4) * t111 + pkin(3)) - m(5) * pkin(3) - t111 * mrSges(5,1) + t108 * mrSges(5,2)) * t112 + (-m(6) * (pkin(6) + qJ(4)) - mrSges(6,3) - m(5) * qJ(4) - mrSges(5,3)) * t109;
t131 = mrSges(5,1) * t108 + mrSges(5,2) * t111;
t125 = t131 * t109;
t152 = qJD(1) * t109;
t138 = t108 * t152;
t85 = t92 * qJD(1);
t59 = qJD(2) * t112 - t109 * t85;
t58 = qJD(4) - t59;
t38 = pkin(4) * t138 + t58;
t52 = t109 * t124;
t54 = -t115 * t138 + t152 * t153;
t172 = m(6) * t38 + mrSges(6,1) * t52 + mrSges(6,2) * t54 + qJD(1) * t125;
t155 = t109 * t111;
t142 = pkin(6) * t155;
t128 = -pkin(4) * t112 - t142;
t148 = qJD(4) * t109;
t113 = cos(pkin(7));
t97 = -pkin(1) * t113 - pkin(2);
t77 = -pkin(3) * t112 - qJ(4) * t109 + t97;
t37 = -qJD(1) * t148 + qJDD(1) * t77 + qJDD(3);
t13 = -t108 * t49 + t111 * t37;
t10 = qJDD(1) * t128 + t13;
t145 = qJDD(1) * t109;
t137 = t108 * t145;
t14 = t108 * t37 + t111 * t49;
t11 = -pkin(6) * t137 + t14;
t50 = qJD(1) * t77 + qJD(3);
t60 = qJD(2) * t109 + t112 * t85;
t20 = -t108 * t60 + t111 * t50;
t12 = qJD(1) * t128 + t20;
t21 = t108 * t50 + t111 * t60;
t15 = -pkin(6) * t138 + t21;
t5 = -t115 * t15 + t117 * t12;
t1 = qJD(5) * t5 + t10 * t115 + t11 * t117;
t6 = t115 * t12 + t117 * t15;
t2 = -qJD(5) * t6 + t10 * t117 - t11 * t115;
t171 = -t2 * mrSges(6,1) + t1 * mrSges(6,2);
t164 = pkin(4) * t108;
t170 = -m(6) * t164 + mrSges(3,2) - mrSges(4,3) - t131;
t168 = t54 / 0.2e1;
t136 = t111 * t145;
t61 = mrSges(5,1) * t137 + mrSges(5,2) * t136;
t7 = -t27 * mrSges(6,1) + mrSges(6,2) * t26;
t166 = -t61 - t7;
t165 = Ifges(6,4) * t54;
t116 = sin(qJ(1));
t102 = t116 * pkin(1);
t118 = cos(qJ(1));
t103 = t118 * pkin(1);
t154 = t111 * t112;
t36 = t108 * t77 + t154 * t92;
t107 = qJ(1) + pkin(7);
t99 = sin(t107);
t161 = t112 * t99;
t160 = pkin(2) * t99 + t102;
t101 = cos(t107);
t158 = t101 * t112;
t157 = t108 * t109;
t156 = t108 * t112;
t149 = qJD(3) * t112;
t147 = m(4) + t167;
t143 = t180 + t181 + t182;
t133 = -mrSges(4,1) * t144 + mrSges(4,2) * t145;
t130 = -t109 * t59 + t112 * t60;
t65 = t111 * t77;
t28 = -t142 + t65 + (-t108 * t92 - pkin(4)) * t112;
t29 = -pkin(6) * t157 + t36;
t8 = -t115 * t29 + t117 * t28;
t9 = t115 * t28 + t117 * t29;
t42 = qJDD(4) - t48;
t127 = -mrSges(5,1) * t112 - mrSges(5,3) * t155;
t126 = mrSges(5,2) * t112 - mrSges(5,3) * t157;
t106 = pkin(9) + qJ(5);
t100 = cos(t106);
t98 = sin(t106);
t90 = qJD(5) - t151;
t83 = qJDD(1) * t97 + qJDD(3);
t76 = t127 * qJD(1);
t75 = t126 * qJD(1);
t72 = t127 * qJDD(1);
t71 = t126 * qJDD(1);
t70 = -t108 * t148 + t111 * t149;
t69 = -t108 * t149 - t111 * t148;
t67 = (t92 + t164) * t109;
t63 = t129 * t109;
t62 = t81 * t109;
t57 = t109 * t174;
t56 = t109 * t74;
t47 = Ifges(6,4) * t52;
t46 = t100 * t158 + t98 * t99;
t45 = -t100 * t99 + t158 * t98;
t44 = t100 * t161 - t101 * t98;
t43 = -t100 * t101 - t161 * t98;
t35 = -t156 * t92 + t65;
t34 = mrSges(6,1) * t90 - mrSges(6,3) * t54;
t33 = -mrSges(6,2) * t90 - mrSges(6,3) * t52;
t32 = pkin(4) * t137 + t42;
t19 = Ifges(6,1) * t54 + Ifges(6,5) * t90 - t47;
t18 = -Ifges(6,2) * t52 + Ifges(6,6) * t90 + t165;
t17 = -mrSges(6,2) * t89 + mrSges(6,3) * t27;
t16 = mrSges(6,1) * t89 - mrSges(6,3) * t26;
t4 = -qJD(5) * t9 - t115 * t70 + t117 * t69;
t3 = qJD(5) * t8 + t115 * t69 + t117 * t70;
t22 = [(Ifges(3,3) + Ifges(2,3) + (0.2e1 * t113 * mrSges(3,1) - 0.2e1 * t110 * mrSges(3,2) + m(3) * (t110 ^ 2 + t113 ^ 2) * pkin(1)) * pkin(1)) * qJDD(1) + (-Ifges(6,1) * t56 + Ifges(6,4) * t57) * t168 + t90 * (-Ifges(6,5) * t56 + Ifges(6,6) * t57) / 0.2e1 - t52 * (-Ifges(6,4) * t56 + Ifges(6,2) * t57) / 0.2e1 + t38 * (-mrSges(6,1) * t57 - mrSges(6,2) * t56) + (t5 * t56 + t57 * t6) * mrSges(6,3) + (-(Ifges(5,4) * t111 - Ifges(5,2) * t108) * t137 + Ifges(4,1) * t145 + m(5) * (qJD(3) * t58 + t42 * t92) + (Ifges(5,1) * t111 - Ifges(5,4) * t108) * t136 + t92 * t61 + (-Ifges(5,5) * t111 + Ifges(5,6) * t108 + Ifges(4,4)) * t144 + t172 * qJD(3)) * t109 + m(4) * (t130 * qJD(3) + t176 * t92 + t83 * t97) + t97 * t133 - t83 * t132 + t14 * t126 + t13 * t127 + t42 * t125 - (-mrSges(6,1) * t32 + mrSges(6,3) * t1 + Ifges(6,4) * t26 + Ifges(6,2) * t27 + Ifges(6,6) * t89) * t62 + t67 * t7 + t36 * t71 + t35 * t72 + t70 * t75 + t69 * t76 - t56 * t19 / 0.2e1 + t57 * t18 / 0.2e1 + t3 * t33 + t4 * t34 + t8 * t16 + t9 * t17 - (mrSges(6,2) * t32 - mrSges(6,3) * t2 + Ifges(6,1) * t26 + Ifges(6,4) * t27 + Ifges(6,5) * t89) * t63 + (-m(3) * t103 - mrSges(2,1) * t118 - t46 * mrSges(6,1) + mrSges(2,2) * t116 + t45 * mrSges(6,2) + (-m(5) + t183) * (pkin(2) * t101 + qJ(3) * t99 + t103) + t170 * t99 + t173 * t101) * g(2) + (-m(3) * t102 - m(5) * t160 - mrSges(2,1) * t116 - t44 * mrSges(6,1) - mrSges(2,2) * t118 - t43 * mrSges(6,2) + t183 * (-qJ(3) * t101 + t160) + (m(5) * qJ(3) - t170) * t101 + t173 * t99) * g(3) + (-t143 / 0.2e1 + Ifges(5,6) * t137 + Ifges(4,4) * t145 - t180 / 0.2e1 - t182 / 0.2e1 - t181 / 0.2e1 - Ifges(5,5) * t136 + (Ifges(5,3) + Ifges(4,2)) * t144 + t171) * t112 + (t184 * t78 + t176) * mrSges(4,3) + m(6) * (t1 * t9 + t2 * t8 + t3 * t6 + t32 * t67 + t4 * t5) + m(5) * (t13 * t35 + t14 * t36 + t20 * t69 + t21 * t70); m(3) * qJDD(2) - t62 * t16 - t63 * t17 - t56 * t33 + t57 * t34 + t166 * t112 + (-t108 * t72 + t111 * t71) * t109 + m(4) * (t109 * t49 + t112 * t48) + m(5) * (-t112 * t42 + (-t108 * t13 + t111 * t14) * t109) + m(6) * (-t1 * t63 - t112 * t32 - t2 * t62 + t5 * t57 - t56 * t6) + (-m(3) - t147) * g(1); t108 * t71 + t111 * t72 - t129 * t16 + t81 * t17 + t178 * t34 + t179 * t33 + m(5) * (t108 * t14 + t111 * t13) + m(4) * t83 + t133 + (g(2) * t101 + g(3) * t99) * t147 + (t1 * t81 - t129 * t2 + t178 * t5 + t179 * t6) * m(6) + ((t108 * t76 - t111 * t75) * t112 - m(5) * (t154 * t21 - t156 * t20) - m(4) * t130 - t184 * mrSges(4,3) * qJD(1) + (-m(5) * t58 - t172) * t109) * qJD(1); t52 * t33 + t54 * t34 + t167 * t112 * g(1) + m(5) * t42 + ((t108 * t75 + t111 * t76 - m(5) * (-t108 * t21 - t111 * t20)) * qJD(1) + t167 * (-g(2) * t99 + g(3) * t101)) * t109 - t166 + (t5 * t54 + t52 * t6 + t32) * m(6); -t38 * (mrSges(6,1) * t54 - mrSges(6,2) * t52) - t54 * (-Ifges(6,1) * t52 - t165) / 0.2e1 + t18 * t168 - t90 * (-Ifges(6,5) * t52 - Ifges(6,6) * t54) / 0.2e1 - t5 * t33 + t6 * t34 - g(2) * (mrSges(6,1) * t43 - mrSges(6,2) * t44) - g(3) * (mrSges(6,1) * t45 + mrSges(6,2) * t46) - g(1) * (-mrSges(6,1) * t98 - mrSges(6,2) * t100) * t109 + (-t5 * t52 + t54 * t6) * mrSges(6,3) + t143 + (-Ifges(6,2) * t54 + t19 - t47) * t52 / 0.2e1 - t171;];
tau = t22;
