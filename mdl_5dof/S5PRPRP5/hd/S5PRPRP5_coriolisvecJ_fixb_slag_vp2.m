% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-05 15:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRPRP5_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP5_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP5_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP5_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP5_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRP5_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRP5_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:37:35
% EndTime: 2019-12-05 15:37:41
% DurationCPUTime: 1.69s
% Computational Cost: add. (1054->203), mult. (2872->265), div. (0->0), fcn. (1869->6), ass. (0->97)
t121 = Ifges(5,1) + Ifges(6,1);
t120 = Ifges(6,4) + Ifges(5,5);
t70 = cos(qJ(2));
t88 = qJD(1) * t70;
t76 = qJD(3) - t88;
t123 = mrSges(5,1) + mrSges(6,1);
t65 = sin(pkin(8));
t67 = sin(qJ(4));
t100 = t65 * t67;
t66 = cos(pkin(8));
t69 = cos(qJ(4));
t99 = t69 * t66;
t51 = -t99 + t100;
t72 = t51 * t70;
t94 = pkin(6) + qJ(3);
t56 = t94 * t65;
t57 = t94 * t66;
t73 = -t69 * t56 - t57 * t67;
t122 = qJD(1) * t72 - qJD(3) * t51 + qJD(4) * t73;
t89 = t65 ^ 2 + t66 ^ 2;
t82 = t89 * mrSges(4,3);
t96 = mrSges(6,2) + mrSges(5,3);
t95 = -Ifges(5,4) + Ifges(6,5);
t119 = Ifges(6,6) - Ifges(5,6);
t26 = -t56 * t67 + t57 * t69;
t52 = t65 * t69 + t66 * t67;
t118 = qJD(4) * t26 + t76 * t52;
t45 = t51 * qJD(2);
t105 = Ifges(6,5) * t45;
t43 = Ifges(5,4) * t45;
t46 = t52 * qJD(2);
t117 = t120 * qJD(4) + t121 * t46 + t105 - t43;
t103 = t45 * mrSges(5,3);
t31 = -t45 * mrSges(6,2) + qJD(4) * mrSges(6,3);
t91 = -qJD(4) * mrSges(5,2) - t103 + t31;
t68 = sin(qJ(2));
t84 = t68 * qJD(1);
t59 = qJD(2) * qJ(3) + t84;
t78 = pkin(6) * qJD(2) + t59;
t37 = t78 * t66;
t104 = t37 * t67;
t36 = t78 * t65;
t11 = -t36 * t69 - t104;
t116 = -t11 + qJD(5);
t115 = -t45 / 0.2e1;
t114 = t45 / 0.2e1;
t112 = t46 / 0.2e1;
t12 = -t36 * t67 + t37 * t69;
t54 = (qJD(3) + t88) * qJD(2);
t3 = qJD(4) * t12 + t52 * t54;
t110 = t73 * t3;
t40 = t52 * t68;
t109 = t3 * t40;
t108 = t3 * t52;
t106 = Ifges(5,4) * t46;
t102 = t46 * mrSges(5,3);
t101 = (-qJD(2) * pkin(2) + t76) * t68;
t98 = -qJD(4) / 0.2e1;
t97 = qJD(4) / 0.2e1;
t47 = t51 * qJD(4);
t38 = qJD(2) * t47;
t48 = t52 * qJD(4);
t39 = qJD(2) * t48;
t15 = t39 * mrSges(6,1) + t38 * mrSges(6,3);
t16 = t39 * mrSges(5,1) - t38 * mrSges(5,2);
t93 = -t15 - t16;
t85 = qJD(4) * t69;
t92 = -t36 * t85 + t54 * t99;
t90 = -t46 * mrSges(6,2) + qJD(4) * t123 - t102;
t87 = qJD(2) * t68;
t86 = qJD(4) * t68;
t22 = mrSges(6,1) * t45 - mrSges(6,3) * t46;
t75 = -mrSges(4,1) * t66 + mrSges(4,2) * t65;
t83 = mrSges(5,1) * t45 + mrSges(5,2) * t46 + qJD(2) * t75 + t22;
t62 = -pkin(3) * t66 - pkin(2);
t81 = t89 * t54;
t80 = t89 * t59;
t79 = qJD(2) * t84;
t77 = t89 * qJD(3);
t49 = qJD(2) * t62 + t76;
t71 = qJD(2) ^ 2;
t42 = Ifges(6,5) * t46;
t41 = t51 * t68;
t24 = pkin(4) * t51 - qJ(5) * t52 + t62;
t21 = pkin(4) * t46 + qJ(5) * t45;
t18 = -Ifges(5,2) * t45 + Ifges(5,6) * qJD(4) + t106;
t17 = Ifges(6,6) * qJD(4) + Ifges(6,3) * t45 + t42;
t14 = t66 * t68 * t85 - t100 * t86 + t46 * t70;
t13 = -qJD(2) * t72 - t52 * t86;
t8 = qJD(4) * qJ(5) + t12;
t7 = -qJD(4) * pkin(4) + t116;
t6 = pkin(4) * t45 - qJ(5) * t46 + t49;
t5 = pkin(4) * t48 + qJ(5) * t47 - qJD(5) * t52;
t4 = pkin(4) * t39 + qJ(5) * t38 - qJD(5) * t46 + t79;
t2 = (-qJD(4) * t37 - t54 * t65) * t67 + t92;
t1 = -t54 * t100 + (qJD(5) - t104) * qJD(4) + t92;
t9 = [t93 * t70 - t90 * t14 + t91 * t13 + t83 * t87 + m(4) * (t68 * t81 + (t101 + (t80 - t84) * t70) * qJD(2)) + m(5) * (-t11 * t14 + t12 * t13 - t2 * t41 + t109 + (t49 - t88) * t87) + m(6) * (-t1 * t41 + t13 * t8 + t14 * t7 - t4 * t70 + t6 * t87 + t109) + (-t68 * mrSges(3,1) + (-mrSges(3,2) + t82) * t70) * t71 + t96 * (-t38 * t40 + t39 * t41); -(t121 * t52 + t95 * t51 - t73 * t96) * t38 + (t95 * t52 + (Ifges(5,2) + Ifges(6,3)) * t51 - t96 * t26) * t39 + (-t70 * qJD(2) * t82 + ((mrSges(5,1) * t51 + mrSges(5,2) * t52 + t75) * qJD(2) - t83) * t68) * qJD(1) + (qJD(2) * t77 + t81) * mrSges(4,3) + t62 * t16 + t4 * (mrSges(6,1) * t51 - mrSges(6,3) * t52) + t5 * t22 + t24 * t15 + (-t2 * t51 + t108) * mrSges(5,3) + (-t1 * t51 + t108) * mrSges(6,2) - t118 * t90 + t91 * t122 + (t1 * t26 + t24 * t4 - t110 + t122 * t8 + t118 * t7 + (t5 - t84) * t6) * m(6) + (-t118 * t11 + t122 * t12 + t2 * t26 - t49 * t84 + t62 * t79 - t110) * m(5) + (-(t70 * t80 + t101) * qJD(1) - pkin(2) * t79 + qJ(3) * t81 + t59 * t77) * m(4) + (-t18 / 0.2e1 + t49 * mrSges(5,1) + t6 * mrSges(6,1) + t17 / 0.2e1 - t12 * mrSges(5,3) - t8 * mrSges(6,2) + Ifges(6,3) * t114 - Ifges(5,2) * t115 + t119 * t97 + t95 * t112) * t48 + (-t49 * mrSges(5,2) - t7 * mrSges(6,2) + t11 * mrSges(5,3) + t6 * mrSges(6,3) - Ifges(5,4) * t115 - Ifges(6,5) * t114 - t120 * t97 - t121 * t112 - t117 / 0.2e1) * t47; t90 * t46 + t91 * t45 + (m(4) + m(5)) * t79 - t71 * t82 - m(5) * (-t11 * t46 - t12 * t45) - m(4) * qJD(2) * t80 - t93 + (t8 * t45 - t7 * t46 + t4) * m(6); -t49 * (mrSges(5,1) * t46 - mrSges(5,2) * t45) + (Ifges(6,3) * t46 - t105) * t115 - t6 * (t46 * mrSges(6,1) + t45 * mrSges(6,3)) + t18 * t112 + qJD(5) * t31 - t21 * t22 + t1 * mrSges(6,3) - t2 * mrSges(5,2) + t119 * t39 - t120 * t38 - t123 * t3 + (t90 + t102) * t12 + (-t91 - t103) * t11 + (pkin(4) * t38 - qJ(5) * t39 + t45 * t7 + t46 * t8) * mrSges(6,2) + (t119 * t46 - t120 * t45) * t98 + (-t3 * pkin(4) + t1 * qJ(5) + t116 * t8 - t7 * t12 - t6 * t21) * m(6) + (-Ifges(5,2) * t46 + t117 - t43) * t114 - (-t121 * t45 - t106 + t17 + t42) * t46 / 0.2e1; -t38 * mrSges(6,2) - qJD(4) * t31 + t46 * t22 + 0.2e1 * (t3 / 0.2e1 + t8 * t98 + t6 * t112) * m(6);];
tauc = t9(:);
