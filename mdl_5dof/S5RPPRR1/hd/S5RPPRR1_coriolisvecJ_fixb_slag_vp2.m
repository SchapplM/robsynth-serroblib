% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
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
% Datum: 2019-12-05 17:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPPRR1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR1_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR1_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR1_coriolisvecJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR1_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:37:59
% EndTime: 2019-12-05 17:38:03
% DurationCPUTime: 0.93s
% Computational Cost: add. (1009->192), mult. (2166->279), div. (0->0), fcn. (1129->4), ass. (0->96)
t75 = pkin(1) + qJ(3);
t121 = qJD(1) * t75;
t70 = qJD(4) + qJD(5);
t69 = qJD(1) * qJ(2) + qJD(3);
t55 = -qJD(1) * pkin(6) + t69;
t77 = sin(qJ(4));
t79 = cos(qJ(4));
t120 = (t77 ^ 2 + t79 ^ 2) * t55;
t100 = qJD(1) * t77;
t76 = sin(qJ(5));
t78 = cos(qJ(5));
t99 = qJD(1) * t79;
t42 = -t78 * t100 - t76 * t99;
t119 = -t42 / 0.2e1;
t118 = t42 / 0.2e1;
t43 = -t76 * t100 + t78 * t99;
t117 = -t43 / 0.2e1;
t116 = t43 / 0.2e1;
t115 = t77 / 0.2e1;
t114 = -t79 / 0.2e1;
t74 = qJ(2) - pkin(6);
t113 = pkin(7) - t74;
t112 = mrSges(6,3) * t42;
t111 = Ifges(5,4) * t77;
t110 = Ifges(5,4) * t79;
t109 = Ifges(6,4) * t43;
t37 = -pkin(7) * t100 + t55 * t77;
t108 = t37 * t76;
t107 = t37 * t78;
t106 = t43 * mrSges(6,3);
t103 = t78 * t79;
t83 = t76 * t77 - t103;
t21 = t70 * t83 * qJD(1);
t46 = -t76 * t79 - t78 * t77;
t105 = t46 * t21;
t25 = t70 * t46;
t20 = t25 * qJD(1);
t104 = t83 * t20;
t92 = qJD(1) * qJD(2);
t95 = qJD(4) * t79;
t34 = t55 * t95 + t77 * t92;
t102 = Ifges(5,5) * qJD(4);
t101 = Ifges(5,6) * qJD(4);
t56 = -qJD(2) + t121;
t97 = qJD(3) * t56;
t96 = qJD(4) * t77;
t94 = qJD(5) * t76;
t93 = qJD(5) * t78;
t91 = qJD(1) * qJD(3);
t90 = pkin(7) * t99;
t52 = t113 * t79;
t89 = m(5) * t74 - mrSges(5,3);
t23 = -mrSges(6,1) * t42 + mrSges(6,2) * t43;
t63 = pkin(4) * t77 + t75;
t44 = t63 * qJD(1) - qJD(2);
t88 = -m(6) * t44 - t23;
t87 = m(3) * qJ(2) + mrSges(4,2) + mrSges(3,3);
t38 = t79 * t55 - t90;
t57 = pkin(4) * t95 + qJD(3);
t86 = mrSges(5,1) * t79 - mrSges(5,2) * t77;
t85 = -t21 * mrSges(6,1) + t20 * mrSges(6,2);
t32 = qJD(4) * pkin(4) + t38;
t11 = t32 * t78 - t108;
t12 = t32 * t76 + t107;
t51 = t113 * t77;
t27 = -t51 * t78 - t52 * t76;
t26 = t51 * t76 - t52 * t78;
t53 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t100;
t54 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t99;
t84 = t53 * t79 - t54 * t77;
t65 = t79 * t92;
t28 = t65 + (pkin(7) * qJD(1) - t55) * t96;
t29 = -qJD(4) * t90 + t34;
t2 = qJD(5) * t11 + t28 * t76 + t29 * t78;
t24 = t70 * t103 - t76 * t96 - t77 * t94;
t3 = -qJD(5) * t12 + t28 * t78 - t29 * t76;
t82 = t11 * t25 + t12 * t24 - t2 * t46 - t3 * t83;
t16 = Ifges(6,2) * t42 + Ifges(6,6) * t70 + t109;
t39 = Ifges(6,4) * t42;
t17 = Ifges(6,1) * t43 + Ifges(6,5) * t70 + t39;
t81 = t3 * mrSges(6,1) - t2 * mrSges(6,2) + t11 * t112 + t16 * t116 - t44 * (mrSges(6,1) * t43 + mrSges(6,2) * t42) + (Ifges(6,1) * t42 - t109) * t117 - t70 * (Ifges(6,5) * t42 - Ifges(6,6) * t43) / 0.2e1 + Ifges(6,6) * t21 + Ifges(6,5) * t20 + (-Ifges(6,2) * t43 + t17 + t39) * t119;
t80 = qJD(1) ^ 2;
t50 = t57 * qJD(1);
t48 = qJD(1) * (mrSges(5,1) * t77 + mrSges(5,2) * t79);
t41 = t102 + (Ifges(5,1) * t79 - t111) * qJD(1);
t40 = t101 + (-Ifges(5,2) * t77 + t110) * qJD(1);
t36 = qJD(2) * t77 - qJD(4) * t52;
t35 = qJD(2) * t79 + t113 * t96;
t33 = -t55 * t96 + t65;
t31 = mrSges(6,1) * t70 - t106;
t30 = -mrSges(6,2) * t70 + t112;
t14 = t38 * t78 - t108;
t13 = -t38 * t76 - t107;
t5 = -t27 * qJD(5) + t35 * t78 - t36 * t76;
t4 = t26 * qJD(5) + t35 * t76 + t36 * t78;
t1 = [t57 * t23 + t63 * t85 + t70 * (Ifges(6,5) * t25 - Ifges(6,6) * t24) / 0.2e1 + t44 * (t24 * mrSges(6,1) + t25 * mrSges(6,2)) + qJD(3) * t48 + t50 * (-t46 * mrSges(6,1) - mrSges(6,2) * t83) + t25 * t17 / 0.2e1 + t4 * t30 + t5 * t31 - t24 * t16 / 0.2e1 + (t24 * t119 + t105) * Ifges(6,2) + (t25 * t116 - t104) * Ifges(6,1) + 0.2e1 * (qJD(3) * mrSges(4,3) + qJD(2) * t87) * qJD(1) + m(5) * (qJD(2) * t120 + t75 * t91 + t97) + m(4) * (qJD(2) * t69 + t97 + (qJ(2) * qJD(2) + qJD(3) * t75) * qJD(1)) + m(6) * (t11 * t5 + t12 * t4 + t2 * t27 + t26 * t3 + t44 * t57 + t50 * t63) + (-t26 * t20 + t27 * t21 - t82) * mrSges(6,3) + (mrSges(5,2) * t91 + qJD(2) * t54 + t89 * t33 + (-t40 / 0.2e1 + t74 * t53 - 0.3e1 / 0.2e1 * Ifges(5,4) * t99 - t101 / 0.2e1 + (t56 + t121) * mrSges(5,1)) * qJD(4)) * t79 + (t24 * t117 + t25 * t118 + t46 * t20 - t21 * t83) * Ifges(6,4) + (mrSges(5,1) * t91 + qJD(2) * t53 + t89 * t34 + (-t56 * mrSges(5,2) - t41 / 0.2e1 - t74 * t54 - t102 / 0.2e1 + (-t75 * mrSges(5,2) + 0.3e1 / 0.2e1 * t111 + (-0.3e1 / 0.2e1 * Ifges(5,1) + 0.3e1 / 0.2e1 * Ifges(5,2)) * t79) * qJD(1)) * qJD(4)) * t77; t42 * t30 - t43 * t31 + 0.2e1 * (-t50 / 0.2e1 + t11 * t117 + t12 * t118) * m(6) - t87 * t80 + (-t77 * t53 - t79 * t54 - t86 * qJD(4) + (-qJD(3) - t120) * m(5) + (-qJD(3) - t69) * m(4)) * qJD(1) - t85; -t80 * mrSges(4,3) + t24 * t30 + t25 * t31 + t84 * qJD(4) + (t104 - t105) * mrSges(6,3) + m(5) * (t33 * t79 + t34 * t77) + m(6) * t82 + (-m(5) * t56 - t48 + (qJD(2) - t56) * m(4) + t88) * qJD(1); (m(6) * (-t11 * t94 + t12 * t93 + t2 * t76 + t3 * t78) + t30 * t93 - t31 * t94 + (-t20 * t78 + t21 * t76) * mrSges(6,3)) * pkin(4) + (-t56 * t86 + t79 * t40 / 0.2e1 + t41 * t115 + ((-Ifges(5,1) * t77 - t110) * t114 + (-Ifges(5,2) * t79 - t111) * t115) * qJD(1) + t88 * t79 * pkin(4) + (-Ifges(5,5) * t77 / 0.2e1 + Ifges(5,6) * t114) * qJD(4)) * qJD(1) - m(6) * (t11 * t13 + t12 * t14) - t84 * t55 + t81 - t14 * t30 - t13 * t31 + t33 * mrSges(5,1) - t34 * mrSges(5,2) + t12 * t106; -t11 * t30 + (t31 + t106) * t12 + t81;];
tauc = t1(:);
