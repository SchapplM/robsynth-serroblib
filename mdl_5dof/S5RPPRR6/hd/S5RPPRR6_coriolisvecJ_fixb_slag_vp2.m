% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
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
% Datum: 2019-12-31 17:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPPRR6_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR6_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR6_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR6_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR6_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR6_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR6_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:57:41
% EndTime: 2019-12-31 17:57:45
% DurationCPUTime: 1.62s
% Computational Cost: add. (2278->226), mult. (5669->320), div. (0->0), fcn. (3886->8), ass. (0->112)
t82 = cos(qJ(5));
t125 = Ifges(6,4) * t82;
t80 = sin(qJ(5));
t101 = -Ifges(6,2) * t80 + t125;
t126 = Ifges(6,4) * t80;
t103 = Ifges(6,1) * t82 - t126;
t104 = mrSges(6,1) * t80 + mrSges(6,2) * t82;
t115 = pkin(6) * qJD(1);
t70 = sin(pkin(8)) * pkin(1) + qJ(3);
t69 = t70 * qJD(1);
t78 = cos(pkin(9));
t73 = t78 * qJD(2);
t76 = sin(pkin(9));
t44 = t73 + (-t69 - t115) * t76;
t53 = t76 * qJD(2) + t78 * t69;
t45 = t78 * t115 + t53;
t81 = sin(qJ(4));
t83 = cos(qJ(4));
t24 = t44 * t81 + t45 * t83;
t22 = qJD(4) * pkin(7) + t24;
t92 = -cos(pkin(8)) * pkin(1) - pkin(3) * t78 - pkin(2);
t59 = qJD(1) * t92 + qJD(3);
t67 = t76 * t81 - t83 * t78;
t60 = t67 * qJD(1);
t68 = t76 * t83 + t78 * t81;
t61 = t68 * qJD(1);
t28 = pkin(4) * t60 - pkin(7) * t61 + t59;
t7 = -t22 * t80 + t28 * t82;
t8 = t22 * t82 + t28 * t80;
t107 = t7 * t82 + t8 * t80;
t136 = t82 / 0.2e1;
t137 = -t80 / 0.2e1;
t47 = qJD(4) * t80 + t61 * t82;
t140 = t47 / 0.2e1;
t127 = Ifges(6,4) * t47;
t46 = qJD(4) * t82 - t61 * t80;
t58 = qJD(5) + t60;
t17 = Ifges(6,2) * t46 + Ifges(6,6) * t58 + t127;
t43 = Ifges(6,4) * t46;
t18 = Ifges(6,1) * t47 + Ifges(6,5) * t58 + t43;
t23 = t44 * t83 - t45 * t81;
t21 = -qJD(4) * pkin(4) - t23;
t99 = Ifges(6,5) * t82 - Ifges(6,6) * t80;
t153 = -t107 * mrSges(6,3) + t58 * t99 / 0.2e1 + t46 * t101 / 0.2e1 + t103 * t140 + t21 * t104 + t17 * t137 + t18 * t136;
t57 = Ifges(5,4) * t60;
t145 = t61 * Ifges(5,1) / 0.2e1 - t57 / 0.2e1;
t152 = t59 * mrSges(5,2) + Ifges(5,5) * qJD(4) + t145 + t153;
t111 = qJD(1) * (t76 ^ 2 + t78 ^ 2);
t151 = mrSges(4,3) * t111;
t88 = t67 * qJD(3);
t14 = -qJD(1) * t88 + qJD(4) * t23;
t62 = t67 * qJD(4);
t54 = qJD(1) * t62;
t63 = t68 * qJD(4);
t55 = qJD(1) * t63;
t35 = pkin(4) * t55 + pkin(7) * t54;
t1 = qJD(5) * t7 + t14 * t82 + t35 * t80;
t2 = -qJD(5) * t8 - t14 * t80 + t35 * t82;
t109 = t1 * t82 - t2 * t80;
t31 = qJD(5) * t46 - t54 * t82;
t19 = mrSges(6,1) * t55 - mrSges(6,3) * t31;
t32 = -qJD(5) * t47 + t54 * t80;
t20 = -mrSges(6,2) * t55 + mrSges(6,3) * t32;
t33 = -mrSges(6,2) * t58 + mrSges(6,3) * t46;
t34 = mrSges(6,1) * t58 - mrSges(6,3) * t47;
t150 = -t80 * t19 + t82 * t20 + m(6) * t109 + (-m(6) * t107 - t80 * t33 - t82 * t34) * qJD(5);
t147 = -t7 * t80 + t8 * t82;
t146 = t2 * mrSges(6,1) - t1 * mrSges(6,2) + Ifges(6,5) * t31 + Ifges(6,6) * t32;
t144 = t31 / 0.2e1;
t143 = t32 / 0.2e1;
t142 = -t46 / 0.2e1;
t141 = -t47 / 0.2e1;
t139 = t55 / 0.2e1;
t138 = -t58 / 0.2e1;
t129 = pkin(6) + t70;
t89 = t68 * qJD(3);
t15 = qJD(1) * t89 + qJD(4) * t24;
t64 = t129 * t76;
t65 = t129 * t78;
t93 = -t83 * t64 - t65 * t81;
t124 = t15 * t93;
t123 = t15 * t67;
t119 = t60 * Ifges(5,2);
t117 = -qJD(4) * mrSges(5,1) - mrSges(6,1) * t46 + mrSges(6,2) * t47 + t61 * mrSges(5,3);
t112 = t55 * mrSges(5,1) - t54 * mrSges(5,2);
t108 = -t1 * t80 - t2 * t82;
t105 = t82 * mrSges(6,1) - t80 * mrSges(6,2);
t102 = Ifges(6,1) * t80 + t125;
t100 = Ifges(6,2) * t82 + t126;
t98 = Ifges(6,5) * t80 + Ifges(6,6) * t82;
t96 = t82 * t33 - t80 * t34;
t36 = pkin(4) * t67 - pkin(7) * t68 + t92;
t40 = -t64 * t81 + t65 * t83;
t12 = t36 * t82 - t40 * t80;
t13 = t36 * t80 + t40 * t82;
t94 = -(-t69 * t76 + t73) * t76 + t53 * t78;
t48 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t60;
t91 = t48 + t96;
t87 = t8 * mrSges(6,2) - t58 * Ifges(6,3) - t47 * Ifges(6,5) - t46 * Ifges(6,6) + Ifges(5,6) * qJD(4) - t119 / 0.2e1 + Ifges(5,4) * t61 - t59 * mrSges(5,1) - t7 * mrSges(6,1);
t51 = Ifges(6,3) * t55;
t42 = pkin(4) * t63 + pkin(7) * t62;
t41 = pkin(4) * t61 + pkin(7) * t60;
t26 = qJD(4) * t40 + t89;
t25 = qJD(4) * t93 - t88;
t11 = t23 * t82 + t41 * t80;
t10 = -t23 * t80 + t41 * t82;
t9 = -mrSges(6,1) * t32 + mrSges(6,2) * t31;
t6 = t31 * Ifges(6,1) + t32 * Ifges(6,4) + t55 * Ifges(6,5);
t5 = t31 * Ifges(6,4) + t32 * Ifges(6,2) + t55 * Ifges(6,6);
t4 = -qJD(5) * t13 - t25 * t80 + t42 * t82;
t3 = qJD(5) * t12 + t25 * t82 + t42 * t80;
t16 = [t92 * t112 + t25 * t48 + t3 * t33 + t4 * t34 - t93 * t9 + t12 * t19 + t13 * t20 + (t119 / 0.2e1 - t87) * t63 - (t145 + t152) * t62 + t117 * t26 + m(6) * (t1 * t13 + t12 * t2 + t21 * t26 + t3 * t8 + t4 * t7 - t124) + m(5) * (t14 * t40 - t23 * t26 + t24 * t25 - t124) + (t23 * t62 - t24 * t63 - t40 * t55 + t54 * t93) * mrSges(5,3) + (m(4) * (t111 * t70 + t94) + 0.2e1 * t151) * qJD(3) + (t51 / 0.2e1 + Ifges(5,4) * t54 - t14 * mrSges(5,3) + (Ifges(6,3) / 0.2e1 + Ifges(5,2)) * t55 + t146) * t67 + (t99 * t139 + t101 * t143 + t103 * t144 - Ifges(5,1) * t54 - Ifges(5,4) * t55 + t6 * t136 + t5 * t137 + (mrSges(5,3) + t104) * t15 + t108 * mrSges(6,3) + (t21 * t105 + t98 * t138 + t100 * t142 + t102 * t141 + t18 * t137 - t82 * t17 / 0.2e1 - t147 * mrSges(6,3)) * qJD(5)) * t68; (-t54 * mrSges(5,3) + t9) * t67 + t117 * t63 - t91 * t62 + m(5) * (-t23 * t63 - t24 * t62 + t123) + m(6) * (-t147 * t62 + t21 * t63 + t123) + (m(5) * t14 - t55 * mrSges(5,3) + t150) * t68; t82 * t19 + t80 * t20 - t117 * t61 + t96 * qJD(5) + t91 * t60 - m(5) * (-t23 * t61 - t24 * t60) + t112 + (-m(4) * t94 - t151) * qJD(1) + (t147 * t58 - t21 * t61 - t108) * m(6); t98 * t139 + t100 * t143 + t102 * t144 + t5 * t136 + t80 * t6 / 0.2e1 - t23 * t48 - Ifges(5,5) * t54 - Ifges(5,6) * t55 - t11 * t33 - t10 * t34 - pkin(4) * t9 - t14 * mrSges(5,2) - t117 * t24 + (-mrSges(5,1) - t105) * t15 + t109 * mrSges(6,3) + (t24 * mrSges(5,3) + t87) * t61 - (t57 / 0.2e1 + t23 * mrSges(5,3) + (-Ifges(5,1) / 0.2e1 + Ifges(5,2) / 0.2e1) * t61 - t152) * t60 + t153 * qJD(5) + (-pkin(4) * t15 - t10 * t7 - t11 * t8 - t21 * t24) * m(6) + t150 * pkin(7); t51 - t21 * (mrSges(6,1) * t47 + mrSges(6,2) * t46) + (Ifges(6,1) * t46 - t127) * t141 + t17 * t140 + (Ifges(6,5) * t46 - Ifges(6,6) * t47) * t138 - t7 * t33 + t8 * t34 + (t46 * t7 + t47 * t8) * mrSges(6,3) + (-Ifges(6,2) * t47 + t18 + t43) * t142 + t146;];
tauc = t16(:);
