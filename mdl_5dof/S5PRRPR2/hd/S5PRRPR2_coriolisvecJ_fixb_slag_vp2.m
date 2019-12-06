% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:18
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRRPR2_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR2_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR2_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR2_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR2_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR2_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR2_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:17:19
% EndTime: 2019-12-05 16:17:22
% DurationCPUTime: 0.85s
% Computational Cost: add. (1003->151), mult. (1836->245), div. (0->0), fcn. (972->6), ass. (0->95)
t59 = cos(pkin(9));
t103 = qJD(4) * t59;
t107 = pkin(2) * qJD(2);
t63 = cos(qJ(3));
t109 = t59 * t63;
t104 = qJ(4) * t59;
t58 = sin(pkin(9));
t49 = -t59 * pkin(4) - t58 * pkin(7) - pkin(3);
t60 = sin(qJ(5));
t62 = cos(qJ(5));
t34 = -t60 * t104 + t62 * t49;
t61 = sin(qJ(3));
t125 = t34 * qJD(5) + t62 * t103 - (t62 * t109 + t60 * t61) * t107;
t35 = t62 * t104 + t60 * t49;
t124 = -t35 * qJD(5) - t60 * t103 - (-t60 * t109 + t61 * t62) * t107;
t55 = t58 ^ 2;
t56 = t59 ^ 2;
t123 = (t55 + t56) * mrSges(5,3);
t80 = -t63 * t107 + qJD(4);
t57 = qJD(2) + qJD(3);
t106 = pkin(2) * qJD(3);
t89 = qJD(2) * t106;
t84 = t63 * t89;
t46 = t57 * qJD(4) + t84;
t117 = t46 * t55;
t90 = t55 * qJD(5) * t57;
t122 = t63 * pkin(2);
t121 = mrSges(6,3) * t58;
t120 = Ifges(6,4) * t60;
t119 = Ifges(6,4) * t62;
t48 = t57 * qJ(4) + t61 * t107;
t36 = -t59 * qJD(1) + t58 * t48;
t118 = t36 * t58;
t116 = t56 * t46;
t115 = t57 * t58;
t76 = mrSges(6,1) * t60 + mrSges(6,2) * t62;
t33 = t76 * t115;
t114 = t58 * t33;
t94 = t63 * t106;
t51 = qJD(4) + t94;
t113 = t58 * t51;
t112 = t59 * t57;
t111 = t59 * t60;
t110 = t59 * t62;
t105 = qJ(4) * t46;
t101 = qJD(5) * t58;
t100 = qJD(5) * t60;
t99 = qJD(5) * t62;
t97 = t60 * t121;
t96 = t62 * t121;
t95 = t61 * t106;
t92 = mrSges(6,3) * t101;
t91 = mrSges(6,3) * t99;
t88 = -t101 / 0.2e1;
t87 = t60 * t92;
t86 = t58 * t91;
t85 = t61 * t89;
t83 = t60 * t88;
t82 = t62 * t88;
t28 = t49 * t57 + t80;
t37 = t58 * qJD(1) + t59 * t48;
t6 = t62 * t28 - t60 * t37;
t7 = t60 * t28 + t62 * t37;
t81 = -t6 * t60 + t62 * t7;
t79 = mrSges(4,1) * t61 + mrSges(4,2) * t63;
t78 = -t59 * mrSges(5,1) + t58 * mrSges(5,2);
t77 = mrSges(6,1) * t62 - mrSges(6,2) * t60;
t75 = -Ifges(6,5) * t60 - Ifges(6,6) * t62;
t74 = t37 * t59 + t118;
t45 = t49 - t122;
t52 = t61 * pkin(2) + qJ(4);
t23 = t52 * t110 + t60 * t45;
t22 = -t52 * t111 + t62 * t45;
t73 = t60 * (-Ifges(6,2) * t62 - t120);
t72 = t62 * (-Ifges(6,1) * t60 - t119);
t71 = (t62 * Ifges(6,1) - t120) * t58;
t70 = (-t60 * Ifges(6,2) + t119) * t58;
t69 = t75 * qJD(5);
t68 = t77 * t101;
t50 = qJD(5) - t112;
t13 = Ifges(6,6) * t50 + t57 * t70;
t14 = Ifges(6,5) * t50 + t57 * t71;
t3 = qJD(5) * t6 + t46 * t110 + t60 * t85;
t4 = -qJD(5) * t7 - t46 * t111 + t62 * t85;
t64 = t3 * (t59 * mrSges(6,2) - t97) + t14 * t83 + t13 * t82 + t36 * t68 + t4 * (-t59 * mrSges(6,1) - t96) + t76 * t117 + t78 * t85 + t6 * t87 + (-t73 + t72) * t90 + (-t112 / 0.2e1 + t50 / 0.2e1) * t58 * t69 + ((-Ifges(6,5) * t59 + t71) * t83 + (-Ifges(6,6) * t59 + t70) * t82) * t57 + (t117 + t116) * mrSges(5,3);
t47 = -t57 * pkin(3) + t80;
t43 = t78 * t57;
t40 = t55 * t105;
t32 = t50 * mrSges(6,1) - t57 * t96;
t31 = -t50 * mrSges(6,2) - t57 * t97;
t30 = t52 * t117;
t29 = t57 * t68;
t9 = -t23 * qJD(5) - t51 * t111 + t62 * t95;
t8 = t22 * qJD(5) + t51 * t110 + t60 * t95;
t1 = [-t59 * t29 + (-t60 ^ 2 - t62 ^ 2) * mrSges(6,3) * t90 + (m(6) * (-t7 * t100 + t3 * t62 - t4 * t60 - t46 * t59 - t6 * t99) - t31 * t100 - t32 * t99) * t58; m(5) * (t52 * t116 + t30 + t74 * t51 + (t47 + (-pkin(3) - t122) * qJD(2)) * t95) + m(6) * (t36 * t113 + t4 * t22 + t3 * t23 + t6 * t9 + t7 * t8 + t30) + t64 - mrSges(4,1) * t85 - t7 * t86 - mrSges(4,2) * t84 + t33 * t113 + t58 * t52 * t29 + t8 * t31 + t9 * t32 + t43 * t95 + (-mrSges(4,1) * t95 - mrSges(4,2) * t94 + t51 * t123 + t22 * t87 - t23 * t86) * t57; (t79 * t107 + (t34 * t60 - t35 * t62) * t92 + t80 * t123) * t57 + (qJ(4) * t29 + qJD(4) * t33 - t7 * t91) * t58 + (-t79 * qJD(3) - t63 * t114 - t43 * t61) * t107 + t124 * t32 + t125 * t31 + t64 + (t80 * t118 + t124 * t6 + t125 * t7 + t3 * t35 + t4 * t34 + t40) * m(6) + (-pkin(3) * t85 + t74 * qJD(4) + t56 * t105 + t40 - (t47 * t61 + t74 * t63) * t107) * m(5); m(5) * t85 + m(6) * (t81 * qJD(5) + t3 * t60 + t4 * t62) + t31 * t99 - t32 * t100 + (-t114 + (-t62 * t31 + t60 * t32) * t59 - m(5) * t74 - m(6) * (t7 * t110 - t6 * t111 + t118) - t57 * t123) * t57; t4 * mrSges(6,1) - t3 * mrSges(6,2) - t6 * t31 + t7 * t32 + (-t36 * t77 + t60 * t14 / 0.2e1 + t62 * t13 / 0.2e1 - t50 * t75 / 0.2e1 + (-t72 / 0.2e1 + t73 / 0.2e1) * t115 + t69 + t81 * mrSges(6,3)) * t115;];
tauc = t1(:);
