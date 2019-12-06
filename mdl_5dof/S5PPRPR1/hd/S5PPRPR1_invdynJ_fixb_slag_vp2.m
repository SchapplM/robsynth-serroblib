% Calculate vector of inverse dynamics joint torques for
% S5PPRPR1
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
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
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
% Datum: 2019-12-05 15:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PPRPR1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR1_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR1_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRPR1_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR1_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR1_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRPR1_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRPR1_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:00:53
% EndTime: 2019-12-05 15:00:59
% DurationCPUTime: 2.30s
% Computational Cost: add. (1167->212), mult. (2736->296), div. (0->0), fcn. (2029->14), ass. (0->107)
t83 = sin(pkin(9));
t86 = cos(pkin(9));
t104 = -mrSges(5,1) * t86 + mrSges(5,2) * t83;
t156 = -mrSges(4,1) + t104;
t87 = cos(pkin(8));
t93 = cos(qJ(3));
t129 = t93 * t87;
t84 = sin(pkin(8));
t91 = sin(qJ(3));
t130 = t84 * t91;
t67 = qJD(1) * t130;
t46 = qJD(1) * t129 - t67;
t100 = qJD(4) - t46;
t70 = pkin(4) * t86 + pkin(3);
t81 = pkin(9) + qJ(5);
t75 = sin(t81);
t77 = cos(t81);
t155 = m(5) * pkin(3) + m(6) * t70 + t77 * mrSges(6,1) - t75 * mrSges(6,2) - t156;
t128 = pkin(6) + qJ(4);
t154 = -m(5) * qJ(4) - m(6) * t128 + mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t153 = t83 ^ 2 + t86 ^ 2;
t117 = qJD(1) * qJD(3);
t152 = qJDD(1) * t91 + t93 * t117;
t85 = sin(pkin(7));
t88 = cos(pkin(7));
t151 = g(1) * t88 + g(2) * t85;
t118 = qJDD(3) * t86;
t119 = qJDD(3) * t83;
t53 = -mrSges(5,1) * t118 + mrSges(5,2) * t119;
t90 = sin(qJ(5));
t92 = cos(qJ(5));
t101 = t83 * t90 - t86 * t92;
t49 = t101 * qJD(5);
t56 = t83 * t92 + t86 * t90;
t30 = -qJD(3) * t49 + qJDD(3) * t56;
t50 = t56 * qJD(5);
t31 = -qJD(3) * t50 - qJDD(3) * t101;
t7 = -t31 * mrSges(6,1) + t30 * mrSges(6,2);
t150 = t53 + t7;
t59 = t128 * t83;
t60 = t128 * t86;
t36 = -t59 * t92 - t60 * t90;
t149 = qJD(5) * t36 - t100 * t101;
t37 = -t59 * t90 + t60 * t92;
t148 = -qJD(5) * t37 - t100 * t56;
t147 = mrSges(5,3) * t153;
t82 = pkin(8) + qJ(3);
t76 = sin(t82);
t146 = t151 * t76;
t121 = qJDD(1) * t93;
t115 = t84 * t121 + t152 * t87;
t116 = qJDD(3) * qJ(4);
t18 = t116 + (qJD(4) - t67) * qJD(3) + t115;
t72 = t86 * qJDD(2);
t10 = -t18 * t83 + t72;
t11 = t83 * qJDD(2) + t86 * t18;
t103 = -t10 * t83 + t11 * t86;
t45 = t101 * qJD(3);
t47 = t56 * qJD(3);
t144 = mrSges(6,1) * t45 + mrSges(6,2) * t47 + t156 * qJD(3);
t140 = t47 / 0.2e1;
t139 = m(3) + m(4);
t138 = m(5) + m(6);
t135 = Ifges(6,4) * t47;
t78 = cos(t82);
t132 = t78 * t85;
t131 = t78 * t88;
t57 = t84 * t93 + t87 * t91;
t48 = t57 * qJD(1);
t43 = qJD(3) * qJ(4) + t48;
t35 = t83 * qJD(2) + t86 * t43;
t126 = pkin(6) * qJD(3);
t123 = t46 * qJD(3);
t113 = t138 + t139;
t109 = t91 * t117;
t74 = t86 * qJD(2);
t28 = t74 + (-t43 - t126) * t83;
t29 = t126 * t86 + t35;
t5 = t28 * t92 - t29 * t90;
t6 = t28 * t90 + t29 * t92;
t102 = -(-t43 * t83 + t74) * t83 + t35 * t86;
t55 = -t129 + t130;
t33 = (-t109 + t121) * t87 - t152 * t84;
t97 = qJDD(4) - t33;
t52 = t57 * qJD(3);
t51 = t55 * qJD(3);
t44 = Ifges(6,4) * t45;
t42 = -qJD(3) * pkin(3) + t100;
t41 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t47;
t40 = -qJD(5) * mrSges(6,2) - mrSges(6,3) * t45;
t38 = -qJD(3) * t70 + t100;
t32 = -t109 * t84 + t115;
t25 = -qJDD(3) * pkin(3) + t97;
t24 = Ifges(6,1) * t47 + Ifges(6,5) * qJD(5) - t44;
t23 = -Ifges(6,2) * t45 + Ifges(6,6) * qJD(5) + t135;
t22 = t101 * t57;
t21 = t56 * t57;
t20 = -qJDD(5) * mrSges(6,2) + mrSges(6,3) * t31;
t19 = qJDD(5) * mrSges(6,1) - mrSges(6,3) * t30;
t14 = -qJDD(3) * t70 + t97;
t9 = pkin(6) * t118 + t11;
t8 = t72 + (-pkin(6) * qJDD(3) - t18) * t83;
t4 = t49 * t57 + t51 * t56;
t3 = t101 * t51 - t50 * t57;
t2 = -qJD(5) * t6 + t8 * t92 - t9 * t90;
t1 = qJD(5) * t5 + t8 * t90 + t9 * t92;
t12 = [-t21 * t19 - t22 * t20 + t3 * t40 + t4 * t41 + (-qJDD(3) * mrSges(4,1) + t150) * t55 + t144 * t52 + (-m(2) - t113) * g(3) + m(5) * (-t102 * t51 + t103 * t57 + t25 * t55 + t42 * t52) + m(4) * (t32 * t57 - t33 * t55 - t46 * t52 - t48 * t51) + m(6) * (-t1 * t22 + t14 * t55 - t2 * t21 + t3 * t6 + t38 * t52 + t4 * t5) + (m(2) + m(3) * (t84 ^ 2 + t87 ^ 2)) * qJDD(1) + (-mrSges(4,2) + t147) * (-qJD(3) * t51 + qJDD(3) * t57); -t101 * t19 + t56 * t20 - t49 * t40 - t50 * t41 + t139 * qJDD(2) + m(5) * (t10 * t86 + t11 * t83) + m(6) * (t1 * t56 - t101 * t2 - t49 * t6 - t5 * t50) + (-g(1) * t85 + g(2) * t88) * t113; (mrSges(6,2) * t14 - mrSges(6,3) * t2 + Ifges(6,1) * t30 + Ifges(6,4) * t31 + Ifges(6,5) * qJDD(5)) * t56 - (-mrSges(6,1) * t14 + mrSges(6,3) * t1 + Ifges(6,4) * t30 + Ifges(6,2) * t31 + Ifges(6,6) * qJDD(5)) * t101 + t25 * t104 + (Ifges(5,4) * t83 + Ifges(5,2) * t86) * t118 + (Ifges(5,1) * t83 + Ifges(5,4) * t86) * t119 + (-Ifges(6,1) * t49 - Ifges(6,4) * t50) * t140 + (t5 * t49 - t6 * t50) * mrSges(6,3) + t38 * (mrSges(6,1) * t50 - mrSges(6,2) * t49) + qJD(5) * (-Ifges(6,5) * t49 - Ifges(6,6) * t50) / 0.2e1 - t45 * (-Ifges(6,4) * t49 - Ifges(6,2) * t50) / 0.2e1 + (t103 + t153 * (qJD(3) * qJD(4) + t116 - t123)) * mrSges(5,3) + (g(3) * t154 + t151 * t155) * t76 + (-g(3) * t155 + t151 * t154) * t78 - t70 * t7 + (t123 - t32) * mrSges(4,2) + t148 * t41 + t149 * t40 + (t1 * t37 - t14 * t70 + t148 * t5 + t149 * t6 + t2 * t36) * m(6) + (-pkin(3) * t25 + qJ(4) * t103 + t100 * t102) * m(5) + (-m(5) * t42 - m(6) * t38 - t144) * t48 - t50 * t23 / 0.2e1 - pkin(3) * t53 - t49 * t24 / 0.2e1 + t33 * mrSges(4,1) + t36 * t19 + t37 * t20 + Ifges(4,3) * qJDD(3); -qJD(3) ^ 2 * t147 + t138 * t78 * g(3) + t45 * t40 + t47 * t41 + (t6 * t45 + t5 * t47 + t14 - t146) * m(6) + (-qJD(3) * t102 - t146 + t25) * m(5) + t150; Ifges(6,5) * t30 + Ifges(6,6) * t31 + Ifges(6,3) * qJDD(5) - t1 * mrSges(6,2) + t2 * mrSges(6,1) - t38 * (mrSges(6,1) * t47 - mrSges(6,2) * t45) - t47 * (-Ifges(6,1) * t45 - t135) / 0.2e1 + t23 * t140 - qJD(5) * (-Ifges(6,5) * t45 - Ifges(6,6) * t47) / 0.2e1 - t5 * t40 + t6 * t41 - g(1) * ((-t131 * t75 + t77 * t85) * mrSges(6,1) + (-t131 * t77 - t75 * t85) * mrSges(6,2)) - g(2) * ((-t132 * t75 - t77 * t88) * mrSges(6,1) + (-t132 * t77 + t75 * t88) * mrSges(6,2)) - g(3) * (-mrSges(6,1) * t75 - mrSges(6,2) * t77) * t76 + (-t45 * t5 + t47 * t6) * mrSges(6,3) + (-Ifges(6,2) * t47 + t24 - t44) * t45 / 0.2e1;];
tau = t12;
