% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
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
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:55
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRPRR5_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR5_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR5_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR5_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR5_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR5_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR5_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:53:26
% EndTime: 2019-12-05 15:53:36
% DurationCPUTime: 2.53s
% Computational Cost: add. (2294->256), mult. (6067->373), div. (0->0), fcn. (4455->8), ass. (0->116)
t102 = sin(pkin(9));
t105 = sin(qJ(4));
t103 = cos(pkin(9));
t108 = cos(qJ(4));
t128 = t103 * t108;
t112 = t102 * t105 - t128;
t124 = qJD(4) * t108;
t109 = cos(qJ(2));
t126 = qJD(1) * t109;
t134 = pkin(6) + qJ(3);
t93 = t134 * t102;
t94 = t134 * t103;
t151 = -t93 * t124 + qJD(3) * t128 + (-qJD(3) * t102 - qJD(4) * t94) * t105 + t112 * t126;
t90 = t102 * t108 + t103 * t105;
t111 = t90 * t109;
t61 = -t105 * t93 + t108 * t94;
t150 = qJD(1) * t111 - qJD(3) * t90 - qJD(4) * t61;
t83 = t90 * qJD(4);
t159 = -pkin(7) * t83 + t151;
t82 = t112 * qJD(4);
t158 = pkin(7) * t82 + t150;
t104 = sin(qJ(5));
t107 = cos(qJ(5));
t106 = sin(qJ(2));
t127 = qJD(1) * t106;
t96 = qJD(2) * qJ(3) + t127;
t118 = pkin(6) * qJD(2) + t96;
t72 = t118 * t102;
t73 = t118 * t103;
t41 = -t105 * t72 + t108 * t73;
t80 = t112 * qJD(2);
t33 = -pkin(7) * t80 + t41;
t130 = t107 * t33;
t40 = -t105 * t73 - t108 * t72;
t81 = t90 * qJD(2);
t32 = -pkin(7) * t81 + t40;
t29 = qJD(4) * pkin(4) + t32;
t10 = t104 * t29 + t130;
t101 = qJD(4) + qJD(5);
t117 = -t104 * t81 - t107 * t80;
t55 = -t104 * t80 + t107 * t81;
t135 = Ifges(6,4) * t55;
t49 = Ifges(6,4) * t117;
t18 = Ifges(6,1) * t55 + Ifges(6,5) * t101 + t49;
t92 = (qJD(3) + t126) * qJD(2);
t25 = -t72 * t124 + t92 * t128 + (-qJD(4) * t73 - t102 * t92) * t105;
t75 = qJD(2) * t83;
t13 = -pkin(7) * t75 + t25;
t26 = -qJD(4) * t41 - t90 * t92;
t74 = qJD(2) * t82;
t14 = pkin(7) * t74 + t26;
t133 = t104 * t33;
t9 = t107 * t29 - t133;
t2 = qJD(5) * t9 + t104 * t14 + t107 * t13;
t22 = t117 * qJD(5) - t104 * t75 - t107 * t74;
t23 = -qJD(5) * t55 + t104 * t74 - t107 * t75;
t3 = -qJD(5) * t10 - t104 * t13 + t107 * t14;
t115 = qJD(3) - t126;
t98 = -pkin(3) * t103 - pkin(2);
t84 = qJD(2) * t98 + t115;
t59 = pkin(4) * t80 + t84;
t157 = t3 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,5) * t22 + Ifges(6,6) * t23 - (Ifges(6,5) * t117 - Ifges(6,6) * t55) * t101 / 0.2e1 + (t10 * t55 + t117 * t9) * mrSges(6,3) - (-Ifges(6,2) * t55 + t18 + t49) * t117 / 0.2e1 - t59 * (mrSges(6,1) * t55 + mrSges(6,2) * t117) - (Ifges(6,1) * t117 - t135) * t55 / 0.2e1;
t17 = Ifges(6,2) * t117 + Ifges(6,6) * t101 + t135;
t155 = t17 / 0.2e1;
t60 = -t105 * t94 - t108 * t93;
t44 = -pkin(7) * t90 + t60;
t45 = -pkin(7) * t112 + t61;
t16 = t104 * t44 + t107 * t45;
t154 = -qJD(5) * t16 - t159 * t104 + t158 * t107;
t15 = -t104 * t45 + t107 * t44;
t153 = qJD(5) * t15 + t158 * t104 + t159 * t107;
t129 = t102 ^ 2 + t103 ^ 2;
t122 = t129 * mrSges(4,3);
t145 = t117 / 0.2e1;
t143 = t55 / 0.2e1;
t141 = -t82 / 0.2e1;
t140 = -t83 / 0.2e1;
t139 = pkin(4) * t83;
t48 = t75 * mrSges(5,1) - t74 * mrSges(5,2);
t6 = -t23 * mrSges(6,1) + t22 * mrSges(6,2);
t137 = -t48 - t6;
t136 = Ifges(5,4) * t81;
t131 = t106 * (-qJD(2) * pkin(2) + t115);
t125 = qJD(2) * t106;
t113 = -mrSges(4,1) * t103 + mrSges(4,2) * t102;
t24 = -mrSges(6,1) * t117 + mrSges(6,2) * t55;
t123 = mrSges(5,1) * t80 + mrSges(5,2) * t81 + qJD(2) * t113 + t24;
t121 = t129 * t92;
t120 = t129 * t96;
t119 = qJD(1) * t125;
t116 = t129 * qJD(3);
t76 = t90 * t106;
t77 = t112 * t106;
t42 = t104 * t77 - t107 * t76;
t43 = -t104 * t76 - t107 * t77;
t57 = -t104 * t90 - t107 * t112;
t58 = -t104 * t112 + t107 * t90;
t110 = qJD(2) ^ 2;
t78 = Ifges(5,4) * t80;
t69 = pkin(4) * t112 + t98;
t66 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t81;
t65 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t80;
t62 = pkin(4) * t75 + t119;
t51 = t81 * Ifges(5,1) + Ifges(5,5) * qJD(4) - t78;
t50 = -t80 * Ifges(5,2) + Ifges(5,6) * qJD(4) + t136;
t47 = -qJD(2) * t111 + t106 * t82;
t46 = -t106 * t83 - t109 * t80;
t37 = mrSges(6,1) * t101 - mrSges(6,3) * t55;
t36 = -mrSges(6,2) * t101 + mrSges(6,3) * t117;
t28 = -qJD(5) * t58 + t104 * t82 - t107 * t83;
t27 = qJD(5) * t57 - t104 * t83 - t107 * t82;
t12 = t107 * t32 - t133;
t11 = -t104 * t32 - t130;
t8 = -qJD(5) * t43 - t104 * t46 + t107 * t47;
t7 = qJD(5) * t42 + t104 * t47 + t107 * t46;
t1 = [t7 * t36 + t8 * t37 + t46 * t65 + t47 * t66 + t137 * t109 + (-t22 * t42 + t23 * t43) * mrSges(6,3) + (-t74 * t76 + t75 * t77) * mrSges(5,3) + t123 * t125 + m(5) * (-t25 * t77 - t26 * t76 + t40 * t47 + t41 * t46 + (t84 - t126) * t125) + m(4) * (t106 * t121 + (t131 + (t120 - t127) * t109) * qJD(2)) + m(6) * (t10 * t7 - t109 * t62 + t125 * t59 + t2 * t43 + t3 * t42 + t8 * t9) + (-t106 * mrSges(3,1) + (-mrSges(3,2) + t122) * t109) * t110; t153 * t36 + (t15 * t3 + t16 * t2 + t62 * t69 + t154 * t9 + (-t127 + t139) * t59 + t153 * t10) * m(6) + t154 * t37 + t98 * t48 + t101 * (Ifges(6,5) * t27 + Ifges(6,6) * t28) / 0.2e1 + t59 * (-mrSges(6,1) * t28 + mrSges(6,2) * t27) + t62 * (-mrSges(6,1) * t57 + mrSges(6,2) * t58) + t69 * t6 + t150 * t66 + (t119 * t98 - t127 * t84 + t150 * t40 + t151 * t41 + t25 * t61 + t26 * t60) * m(5) + t151 * t65 + t84 * (mrSges(5,1) * t83 - mrSges(5,2) * t82) + qJD(4) * (-Ifges(5,5) * t82 - Ifges(5,6) * t83) / 0.2e1 + t28 * t155 + t27 * t18 / 0.2e1 + (-(t109 * t120 + t131) * qJD(1) - pkin(2) * t119 + qJ(3) * t121 + t96 * t116) * m(4) + (t10 * t28 - t15 * t22 + t16 * t23 + t2 * t57 - t27 * t9 - t3 * t58) * mrSges(6,3) + (-t109 * qJD(2) * t122 + ((mrSges(5,1) * t112 + mrSges(5,2) * t90 + t113) * qJD(2) - t123) * t106) * qJD(1) + (t141 * t81 - t74 * t90) * Ifges(5,1) + (t112 * t75 - t140 * t80) * Ifges(5,2) + (-t112 * t25 - t26 * t90 + t40 * t82 - t41 * t83 + t60 * t74 - t61 * t75) * mrSges(5,3) + (t112 * t74 + t140 * t81 - t141 * t80 - t75 * t90) * Ifges(5,4) + (qJD(2) * t116 + t121) * mrSges(4,3) + t24 * t139 + (t143 * t27 + t58 * t22) * Ifges(6,1) + (t145 * t28 + t57 * t23) * Ifges(6,2) + (t143 * t28 + t145 * t27 + t57 * t22 + t58 * t23) * Ifges(6,4) + t50 * t140 + t51 * t141; -t117 * t36 + t55 * t37 + t80 * t65 + t81 * t66 + (m(4) + m(5)) * t119 - t110 * t122 - m(5) * (-t40 * t81 - t41 * t80) - m(4) * qJD(2) * t120 - t137 + (-t10 * t117 + t55 * t9 + t62) * m(6); t81 * t50 / 0.2e1 - Ifges(5,5) * t74 - Ifges(5,6) * t75 - t40 * t65 + t41 * t66 - t12 * t36 - t11 * t37 + (-t40 * t80 + t41 * t81) * mrSges(5,3) - t84 * (mrSges(5,1) * t81 - mrSges(5,2) * t80) - qJD(4) * (-Ifges(5,5) * t80 - Ifges(5,6) * t81) / 0.2e1 - t81 * (-Ifges(5,1) * t80 - t136) / 0.2e1 - t25 * mrSges(5,2) + t26 * mrSges(5,1) - m(6) * (t10 * t12 + t11 * t9) + t55 * t155 + (-Ifges(5,2) * t81 + t51 - t78) * t80 / 0.2e1 + (-t81 * t24 + (-t104 * t37 + t107 * t36) * qJD(5) + (t104 * t23 - t107 * t22) * mrSges(6,3) + (t104 * t2 + t107 * t3 - t59 * t81 + (t10 * t107 - t104 * t9) * qJD(5)) * m(6)) * pkin(4) + t157; t10 * t37 + t17 * t143 - t9 * t36 + t157;];
tauc = t1(:);
