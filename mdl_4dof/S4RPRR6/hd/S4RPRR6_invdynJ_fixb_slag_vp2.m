% Calculate vector of inverse dynamics joint torques for
% S4RPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RPRR6_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR6_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR6_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR6_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR6_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR6_invdynJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR6_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR6_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR6_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:52:29
% EndTime: 2019-12-31 16:52:36
% DurationCPUTime: 4.23s
% Computational Cost: add. (2025->265), mult. (5141->367), div. (0->0), fcn. (3710->12), ass. (0->116)
t136 = qJD(1) * qJD(2);
t87 = qJDD(1) * qJ(2) + t136;
t104 = sin(pkin(7));
t105 = cos(pkin(7));
t139 = t104 ^ 2 + t105 ^ 2;
t102 = pkin(7) + qJ(3);
t97 = qJ(4) + t102;
t90 = sin(t97);
t91 = cos(t97);
t125 = t91 * mrSges(5,1) - t90 * mrSges(5,2);
t95 = sin(t102);
t96 = cos(t102);
t171 = -t96 * mrSges(4,1) + t95 * mrSges(4,2) - t125;
t103 = qJD(3) + qJD(4);
t107 = sin(qJ(4));
t110 = cos(qJ(4));
t108 = sin(qJ(3));
t111 = cos(qJ(3));
t140 = t105 * t111;
t77 = -t104 * t108 + t140;
t69 = t77 * qJD(1);
t78 = t104 * t111 + t105 * t108;
t70 = t78 * qJD(1);
t124 = -t107 * t70 + t110 * t69;
t71 = t77 * qJD(3);
t48 = qJD(1) * t71 + qJDD(1) * t78;
t72 = t78 * qJD(3);
t49 = -qJD(1) * t72 + qJDD(1) * t77;
t12 = qJD(4) * t124 + t107 * t49 + t110 * t48;
t46 = t107 * t69 + t110 * t70;
t13 = -qJD(4) * t46 - t107 * t48 + t110 * t49;
t145 = pkin(5) + qJ(2);
t84 = t145 * t104;
t79 = qJD(1) * t84;
t85 = t145 * t105;
t80 = qJD(1) * t85;
t53 = -t108 * t79 + t111 * t80;
t32 = pkin(6) * t69 + t53;
t144 = t107 * t32;
t142 = t108 * t80;
t52 = -t111 * t79 - t142;
t31 = -pkin(6) * t70 + t52;
t30 = qJD(3) * pkin(3) + t31;
t15 = t110 * t30 - t144;
t150 = Ifges(5,4) * t46;
t141 = t110 * t32;
t16 = t107 * t30 + t141;
t123 = pkin(5) * qJDD(1) + t87;
t64 = t123 * t104;
t65 = t123 * t105;
t25 = -qJD(3) * t53 - t108 * t65 - t111 * t64;
t11 = qJDD(3) * pkin(3) - pkin(6) * t48 + t25;
t138 = qJD(3) * t111;
t24 = -qJD(3) * t142 - t108 * t64 + t111 * t65 - t79 * t138;
t14 = pkin(6) * t49 + t24;
t2 = qJD(4) * t15 + t107 * t11 + t110 * t14;
t40 = Ifges(5,4) * t124;
t22 = Ifges(5,1) * t46 + Ifges(5,5) * t103 + t40;
t3 = -qJD(4) * t16 - t107 * t14 + t11 * t110;
t92 = pkin(2) * t105 + pkin(1);
t83 = -qJD(1) * t92 + qJD(2);
t54 = -pkin(3) * t69 + t83;
t98 = qJDD(3) + qJDD(4);
t170 = t3 * mrSges(5,1) - t2 * mrSges(5,2) + Ifges(5,5) * t12 + Ifges(5,6) * t13 + Ifges(5,3) * t98 - (Ifges(5,5) * t124 - Ifges(5,6) * t46) * t103 / 0.2e1 + (t124 * t15 + t16 * t46) * mrSges(5,3) - (-Ifges(5,2) * t46 + t22 + t40) * t124 / 0.2e1 - t54 * (mrSges(5,1) * t46 + mrSges(5,2) * t124) - (Ifges(5,1) * t124 - t150) * t46 / 0.2e1;
t109 = sin(qJ(1));
t112 = cos(qJ(1));
t160 = g(1) * t112 + g(2) * t109;
t119 = -mrSges(3,1) * t105 + mrSges(3,2) * t104;
t168 = m(3) * pkin(1) + mrSges(2,1) + m(5) * (pkin(3) * t96 + t92) + m(4) * t92 - t119 - t171;
t128 = m(3) * qJ(2) + mrSges(3,3);
t167 = -m(4) * t145 + mrSges(2,2) - mrSges(4,3) - t128 + m(5) * (-pkin(6) - t145) - mrSges(5,3);
t21 = Ifges(5,2) * t124 + Ifges(5,6) * t103 + t150;
t165 = t21 / 0.2e1;
t157 = t46 / 0.2e1;
t155 = t70 / 0.2e1;
t154 = pkin(3) * t72;
t151 = Ifges(4,4) * t70;
t56 = -t108 * t84 + t111 * t85;
t134 = qJDD(1) * t104;
t133 = qJDD(1) * t105;
t131 = -t13 * mrSges(5,1) + t12 * mrSges(5,2);
t126 = -t49 * mrSges(4,1) + t48 * mrSges(4,2);
t55 = -t108 * t85 - t111 * t84;
t122 = -mrSges(3,1) * t133 + mrSges(3,2) * t134;
t120 = mrSges(5,1) * t90 + mrSges(5,2) * t91;
t38 = -pkin(6) * t78 + t55;
t39 = pkin(6) * t77 + t56;
t19 = -t107 * t39 + t110 * t38;
t20 = t107 * t38 + t110 * t39;
t50 = -t107 * t78 + t110 * t77;
t51 = t107 * t77 + t110 * t78;
t82 = -qJDD(1) * t92 + qJDD(2);
t36 = -t84 * t138 + qJD(2) * t140 + (-qJD(2) * t104 - qJD(3) * t85) * t108;
t37 = -qJD(2) * t78 - qJD(3) * t56;
t94 = -qJDD(1) * pkin(1) + qJDD(2);
t66 = Ifges(4,4) * t69;
t62 = -pkin(3) * t77 - t92;
t58 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t70;
t57 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t69;
t42 = t70 * Ifges(4,1) + Ifges(4,5) * qJD(3) + t66;
t41 = t69 * Ifges(4,2) + Ifges(4,6) * qJD(3) + t151;
t35 = mrSges(5,1) * t103 - mrSges(5,3) * t46;
t34 = -mrSges(5,2) * t103 + mrSges(5,3) * t124;
t33 = -pkin(3) * t49 + t82;
t29 = -pkin(6) * t71 + t37;
t28 = -pkin(6) * t72 + t36;
t27 = -qJD(4) * t51 - t107 * t71 - t110 * t72;
t26 = qJD(4) * t50 - t107 * t72 + t110 * t71;
t23 = -mrSges(5,1) * t124 + mrSges(5,2) * t46;
t18 = t110 * t31 - t144;
t17 = -t107 * t31 - t141;
t7 = -mrSges(5,2) * t98 + mrSges(5,3) * t13;
t6 = mrSges(5,1) * t98 - mrSges(5,3) * t12;
t5 = -qJD(4) * t20 - t107 * t28 + t110 * t29;
t4 = qJD(4) * t19 + t107 * t29 + t110 * t28;
t1 = [(t167 * t109 - t168 * t112) * g(2) + (-mrSges(5,1) * t33 + mrSges(5,3) * t2 + Ifges(5,4) * t12 + Ifges(5,2) * t13 + Ifges(5,6) * t98) * t50 + (Ifges(5,1) * t26 + Ifges(5,4) * t27) * t157 + (Ifges(3,4) * t104 + Ifges(3,2) * t105) * t133 + (Ifges(3,1) * t104 + Ifges(3,4) * t105) * t134 + m(4) * (t24 * t56 + t25 * t55 + t36 * t53 + t37 * t52 - t82 * t92) + (mrSges(5,2) * t33 - mrSges(5,3) * t3 + Ifges(5,1) * t12 + Ifges(5,4) * t13 + Ifges(5,5) * t98) * t51 + m(3) * (-pkin(1) * t94 + (t136 + t87) * qJ(2) * t139) + t23 * t154 + (t55 * mrSges(4,1) - t56 * mrSges(4,2)) * qJDD(3) + t62 * t131 - t92 * t126 - pkin(1) * t122 + t94 * t119 + (t168 * t109 + t167 * t112) * g(1) + (Ifges(4,1) * t71 - Ifges(4,4) * t72) * t155 + t83 * (mrSges(4,1) * t72 + mrSges(4,2) * t71) + t69 * (Ifges(4,4) * t71 - Ifges(4,2) * t72) / 0.2e1 + qJD(3) * (Ifges(4,5) * t71 - Ifges(4,6) * t72) / 0.2e1 + m(5) * (t15 * t5 + t154 * t54 + t16 * t4 + t19 * t3 + t2 * t20 + t33 * t62) + t124 * (Ifges(5,4) * t26 + Ifges(5,2) * t27) / 0.2e1 + (-t55 * t48 + t56 * t49 - t52 * t71 - t53 * t72) * mrSges(4,3) + t103 * (Ifges(5,5) * t26 + Ifges(5,6) * t27) / 0.2e1 - t72 * t41 / 0.2e1 + t71 * t42 / 0.2e1 + t37 * t58 + t54 * (-mrSges(5,1) * t27 + mrSges(5,2) * t26) + t36 * t57 + t4 * t34 + t5 * t35 + t27 * t165 + t19 * t6 + t20 * t7 + t26 * t22 / 0.2e1 + (-t15 * t26 + t16 * t27) * mrSges(5,3) + 0.2e1 * t139 * t87 * mrSges(3,3) + (-mrSges(4,1) * t82 + mrSges(4,3) * t24 + Ifges(4,4) * t48 + Ifges(4,2) * t49 + Ifges(4,6) * qJDD(3)) * t77 + (t82 * mrSges(4,2) - t25 * mrSges(4,3) + Ifges(4,1) * t48 + Ifges(4,4) * t49 + Ifges(4,5) * qJDD(3)) * t78 + Ifges(2,3) * qJDD(1); m(3) * t94 - t124 * t34 + t46 * t35 - t69 * t57 + t70 * t58 + t122 + t126 + t131 + (-g(1) * t109 + g(2) * t112) * (m(3) + m(4) + m(5)) - t128 * t139 * qJD(1) ^ 2 + (-t124 * t16 + t15 * t46 + t33) * m(5) + (t52 * t70 - t53 * t69 + t82) * m(4); t160 * (mrSges(4,1) * t95 + mrSges(4,2) * t96 + t120) + t41 * t155 - m(5) * (t15 * t17 + t16 * t18) - t70 * (Ifges(4,1) * t69 - t151) / 0.2e1 - (-Ifges(4,2) * t70 + t42 + t66) * t69 / 0.2e1 + t46 * t165 + t170 + (t107 * t7 + t110 * t6 - t70 * t23 + (-g(3) * t96 + t107 * t2 + t110 * t3 + t160 * t95 - t54 * t70) * m(5) + (-t107 * t35 + t110 * t34 + (-t107 * t15 + t110 * t16) * m(5)) * qJD(4)) * pkin(3) + t171 * g(3) - t83 * (mrSges(4,1) * t70 + mrSges(4,2) * t69) - qJD(3) * (Ifges(4,5) * t69 - Ifges(4,6) * t70) / 0.2e1 + t53 * t58 - t52 * t57 + Ifges(4,5) * t48 + Ifges(4,6) * t49 - t18 * t34 - t17 * t35 - t24 * mrSges(4,2) + t25 * mrSges(4,1) + (t52 * t69 + t53 * t70) * mrSges(4,3) + Ifges(4,3) * qJDD(3); -g(3) * t125 + t160 * t120 - t15 * t34 + t21 * t157 + t16 * t35 + t170;];
tau = t1;
