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
% Datum: 2019-12-05 17:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 17:28:46
% EndTime: 2019-12-05 17:28:53
% DurationCPUTime: 3.62s
% Computational Cost: add. (1742->280), mult. (3851->397), div. (0->0), fcn. (2564->14), ass. (0->137)
t106 = sin(pkin(7));
t92 = pkin(1) * t106 + qJ(3);
t78 = qJD(1) * qJD(3) + qJDD(1) * t92;
t105 = sin(pkin(8));
t108 = cos(pkin(8));
t184 = t105 ^ 2 + t108 ^ 2;
t167 = m(5) + m(6);
t183 = -m(6) - m(4);
t104 = sin(pkin(9));
t111 = sin(qJ(5));
t107 = cos(pkin(9));
t113 = cos(qJ(5));
t151 = t107 * t113;
t126 = t104 * t111 - t151;
t81 = t104 * t113 + t107 * t111;
t74 = t81 * qJD(5);
t26 = (-qJD(1) * t74 - qJDD(1) * t126) * t105;
t182 = Ifges(6,5) * t26;
t174 = t126 * qJD(5);
t27 = (qJD(1) * t174 - qJDD(1) * t81) * t105;
t181 = Ifges(6,6) * t27;
t142 = qJDD(1) * t108;
t89 = qJDD(5) - t142;
t180 = Ifges(6,3) * t89;
t149 = qJD(1) * t108;
t179 = t126 * t149 - t174;
t121 = qJD(1) * t81;
t178 = t108 * t121 - t74;
t48 = qJDD(2) * t108 - t105 * t78;
t49 = qJDD(2) * t105 + t108 * t78;
t176 = -t105 * t48 + t108 * t49;
t129 = -mrSges(4,1) * t108 + mrSges(4,2) * t105;
t133 = -qJ(4) * t105 - pkin(2);
t173 = -m(6) * (-t108 * (pkin(4) * t107 + pkin(3)) - pkin(2)) + mrSges(3,1) + m(4) * pkin(2) - t129 - m(5) * t133 - (-m(5) * pkin(3) - t107 * mrSges(5,1) + t104 * mrSges(5,2)) * t108 + (-m(6) * (-pkin(6) - qJ(4)) + mrSges(6,3) + mrSges(5,3)) * t105;
t128 = mrSges(5,1) * t104 + mrSges(5,2) * t107;
t122 = t128 * t105;
t150 = qJD(1) * t105;
t137 = t104 * t150;
t85 = t92 * qJD(1);
t59 = qJD(2) * t108 - t105 * t85;
t58 = qJD(4) - t59;
t38 = pkin(4) * t137 + t58;
t52 = t105 * t121;
t54 = -t111 * t137 + t150 * t151;
t172 = m(6) * t38 + mrSges(6,1) * t52 + mrSges(6,2) * t54 + qJD(1) * t122;
t153 = t105 * t107;
t140 = pkin(6) * t153;
t125 = -pkin(4) * t108 - t140;
t146 = qJD(4) * t105;
t109 = cos(pkin(7));
t164 = pkin(1) * t109;
t77 = -pkin(3) * t108 + t133 - t164;
t37 = -qJD(1) * t146 + qJDD(1) * t77 + qJDD(3);
t13 = -t104 * t49 + t107 * t37;
t10 = qJDD(1) * t125 + t13;
t143 = qJDD(1) * t105;
t135 = t104 * t143;
t14 = t104 * t37 + t107 * t49;
t11 = -pkin(6) * t135 + t14;
t50 = qJD(1) * t77 + qJD(3);
t60 = qJD(2) * t105 + t108 * t85;
t20 = -t104 * t60 + t107 * t50;
t12 = qJD(1) * t125 + t20;
t21 = t104 * t50 + t107 * t60;
t15 = -pkin(6) * t137 + t21;
t5 = -t111 * t15 + t113 * t12;
t1 = qJD(5) * t5 + t10 * t111 + t11 * t113;
t6 = t111 * t12 + t113 * t15;
t2 = -qJD(5) * t6 + t10 * t113 - t11 * t111;
t171 = -t2 * mrSges(6,1) + t1 * mrSges(6,2);
t161 = pkin(4) * t104;
t170 = -m(6) * t161 + mrSges(3,2) - mrSges(4,3) - t128;
t168 = t54 / 0.2e1;
t134 = t107 * t143;
t61 = mrSges(5,1) * t135 + mrSges(5,2) * t134;
t7 = -t27 * mrSges(6,1) + mrSges(6,2) * t26;
t166 = -t61 - t7;
t165 = Ifges(6,4) * t54;
t112 = sin(qJ(1));
t163 = pkin(1) * t112;
t114 = cos(qJ(1));
t162 = pkin(1) * t114;
t152 = t107 * t108;
t36 = t104 * t77 + t152 * t92;
t103 = qJ(1) + pkin(7);
t97 = sin(t103);
t158 = t108 * t97;
t99 = cos(t103);
t157 = t108 * t99;
t155 = t104 * t105;
t154 = t104 * t108;
t147 = qJD(3) * t108;
t145 = -m(4) - t167;
t141 = t180 + t181 + t182;
t95 = -pkin(2) - t164;
t131 = -mrSges(4,1) * t142 + mrSges(4,2) * t143;
t127 = -t105 * t59 + t108 * t60;
t65 = t107 * t77;
t28 = -t140 + t65 + (-t104 * t92 - pkin(4)) * t108;
t29 = -pkin(6) * t155 + t36;
t8 = -t111 * t29 + t113 * t28;
t9 = t111 * t28 + t113 * t29;
t42 = qJDD(4) - t48;
t124 = -mrSges(5,1) * t108 - mrSges(5,3) * t153;
t123 = mrSges(5,2) * t108 - mrSges(5,3) * t155;
t102 = pkin(9) + qJ(5);
t98 = cos(t102);
t96 = sin(t102);
t90 = qJD(5) - t149;
t83 = qJDD(1) * t95 + qJDD(3);
t76 = t124 * qJD(1);
t75 = t123 * qJD(1);
t72 = t124 * qJDD(1);
t71 = t123 * qJDD(1);
t70 = -t104 * t146 + t107 * t147;
t69 = -t104 * t147 - t107 * t146;
t67 = (t92 + t161) * t105;
t63 = t126 * t105;
t62 = t81 * t105;
t57 = t105 * t174;
t56 = t105 * t74;
t47 = Ifges(6,4) * t52;
t46 = -t157 * t98 - t96 * t97;
t45 = t157 * t96 - t97 * t98;
t44 = t158 * t98 - t96 * t99;
t43 = t158 * t96 + t98 * t99;
t35 = -t154 * t92 + t65;
t34 = mrSges(6,1) * t90 - mrSges(6,3) * t54;
t33 = -mrSges(6,2) * t90 - mrSges(6,3) * t52;
t32 = pkin(4) * t135 + t42;
t19 = Ifges(6,1) * t54 + Ifges(6,5) * t90 - t47;
t18 = -Ifges(6,2) * t52 + Ifges(6,6) * t90 + t165;
t17 = -mrSges(6,2) * t89 + mrSges(6,3) * t27;
t16 = mrSges(6,1) * t89 - mrSges(6,3) * t26;
t4 = -qJD(5) * t9 - t111 * t70 + t113 * t69;
t3 = qJD(5) * t8 + t111 * t69 + t113 * t70;
t22 = [(Ifges(3,3) + Ifges(2,3) + (0.2e1 * t109 * mrSges(3,1) - 0.2e1 * t106 * mrSges(3,2) + m(3) * (t106 ^ 2 + t109 ^ 2) * pkin(1)) * pkin(1)) * qJDD(1) + (-Ifges(6,1) * t56 + Ifges(6,4) * t57) * t168 + t90 * (-Ifges(6,5) * t56 + Ifges(6,6) * t57) / 0.2e1 - t52 * (-Ifges(6,4) * t56 + Ifges(6,2) * t57) / 0.2e1 + t38 * (-mrSges(6,1) * t57 - mrSges(6,2) * t56) + (t5 * t56 + t57 * t6) * mrSges(6,3) + (m(3) * t163 + mrSges(2,1) * t112 + t44 * mrSges(6,1) + mrSges(2,2) * t114 - t43 * mrSges(6,2) + (-m(5) + t183) * (qJ(3) * t99 - t163) + t170 * t99 + t173 * t97) * g(3) + (mrSges(2,1) * t114 - t46 * mrSges(6,1) - mrSges(2,2) * t112 - t45 * mrSges(6,2) + (m(5) + m(3)) * t162 + t183 * (-qJ(3) * t97 - t162) + (m(5) * qJ(3) - t170) * t97 + t173 * t99) * g(2) + (Ifges(4,4) * t143 - t141 / 0.2e1 - t180 / 0.2e1 - t182 / 0.2e1 - t181 / 0.2e1 + Ifges(5,6) * t135 - Ifges(5,5) * t134 + (Ifges(4,2) + Ifges(5,3)) * t142 + t171) * t108 - (mrSges(6,2) * t32 - mrSges(6,3) * t2 + Ifges(6,1) * t26 + Ifges(6,4) * t27 + Ifges(6,5) * t89) * t63 - (-mrSges(6,1) * t32 + mrSges(6,3) * t1 + Ifges(6,4) * t26 + Ifges(6,2) * t27 + Ifges(6,6) * t89) * t62 + m(5) * (t13 * t35 + t14 * t36 + t20 * t69 + t21 * t70) + m(6) * (t1 * t9 + t2 * t8 + t3 * t6 + t32 * t67 + t4 * t5) + t42 * t122 + t36 * t71 + t35 * t72 + t70 * t75 + t69 * t76 + t67 * t7 + m(4) * (t127 * qJD(3) + t176 * t92 + t83 * t95) - t56 * t19 / 0.2e1 + t57 * t18 / 0.2e1 + t3 * t33 + t4 * t34 + t8 * t16 + t9 * t17 + (Ifges(4,1) * t143 + m(5) * (qJD(3) * t58 + t42 * t92) + t92 * t61 - (Ifges(5,4) * t107 - Ifges(5,2) * t104) * t135 + (Ifges(5,1) * t107 - Ifges(5,4) * t104) * t134 + (-Ifges(5,5) * t107 + Ifges(5,6) * t104 + Ifges(4,4)) * t142 + t172 * qJD(3)) * t105 + t95 * t131 + t83 * t129 + t14 * t123 + t13 * t124 + (t184 * t78 + t176) * mrSges(4,3); m(3) * qJDD(2) - t62 * t16 - t63 * t17 - t56 * t33 + t57 * t34 + t166 * t108 + (-t104 * t72 + t107 * t71) * t105 + m(4) * (t105 * t49 + t108 * t48) + m(5) * (-t108 * t42 + (-t104 * t13 + t107 * t14) * t105) + m(6) * (-t1 * t63 - t108 * t32 - t2 * t62 + t5 * t57 - t56 * t6) + (-m(3) + t145) * g(1); t104 * t71 + t107 * t72 - t126 * t16 + t81 * t17 + t178 * t34 + t179 * t33 + m(5) * (t14 * t104 + t13 * t107) + m(4) * t83 + t131 + (g(2) * t99 + g(3) * t97) * t145 + (t1 * t81 - t126 * t2 + t178 * t5 + t179 * t6) * m(6) + ((t104 * t76 - t107 * t75) * t108 - m(5) * (t152 * t21 - t154 * t20) - m(4) * t127 - t184 * mrSges(4,3) * qJD(1) + (-m(5) * t58 - t172) * t105) * qJD(1); t52 * t33 + t54 * t34 + t167 * t108 * g(1) + m(5) * t42 + ((t104 * t75 + t107 * t76 - m(5) * (-t104 * t21 - t107 * t20)) * qJD(1) + t167 * (g(2) * t97 - g(3) * t99)) * t105 - t166 + (t5 * t54 + t52 * t6 + t32) * m(6); -t38 * (mrSges(6,1) * t54 - mrSges(6,2) * t52) - t54 * (-Ifges(6,1) * t52 - t165) / 0.2e1 + t18 * t168 - t90 * (-Ifges(6,5) * t52 - Ifges(6,6) * t54) / 0.2e1 - t5 * t33 + t6 * t34 - g(2) * (mrSges(6,1) * t43 + mrSges(6,2) * t44) - g(3) * (-mrSges(6,1) * t45 + mrSges(6,2) * t46) - g(1) * (-mrSges(6,1) * t96 - mrSges(6,2) * t98) * t105 + (-t5 * t52 + t54 * t6) * mrSges(6,3) + t141 + (-Ifges(6,2) * t54 + t19 - t47) * t52 / 0.2e1 - t171;];
tau = t22;
