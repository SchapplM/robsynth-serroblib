% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
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
% Datum: 2019-12-05 18:19
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPPR1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR1_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR1_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR1_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR1_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:18:10
% EndTime: 2019-12-05 18:18:13
% DurationCPUTime: 0.94s
% Computational Cost: add. (1537->162), mult. (3042->235), div. (0->0), fcn. (1864->8), ass. (0->101)
t91 = sin(pkin(9));
t93 = cos(pkin(9));
t120 = t91 ^ 2 + t93 ^ 2;
t141 = t120 * mrSges(5,3);
t90 = qJD(1) + qJD(2);
t118 = qJD(4) * t90;
t119 = qJD(1) * pkin(1);
t111 = qJD(2) * t119;
t98 = cos(qJ(2));
t105 = t98 * t111;
t96 = sin(qJ(2));
t106 = t96 * t111;
t92 = sin(pkin(8));
t94 = cos(pkin(8));
t57 = t94 * t105 - t92 * t106;
t45 = t57 + t118;
t110 = t45 * t120;
t147 = mrSges(5,1) * t93 - mrSges(5,2) * t91 + mrSges(4,1);
t114 = t98 * t119;
t115 = t96 * t119;
t78 = t92 * t115;
t64 = t94 * t114 - t78;
t146 = qJD(4) - t64;
t95 = sin(qJ(5));
t97 = cos(qJ(5));
t72 = -t91 * t95 + t93 * t97;
t67 = t72 * qJD(5);
t145 = t67 / 0.2e1;
t73 = t91 * t97 + t93 * t95;
t68 = t73 * qJD(5);
t144 = -t68 / 0.2e1;
t51 = t90 * t67;
t143 = t51 * t72;
t52 = t90 * t68;
t142 = t52 * t73;
t82 = pkin(2) * t92 + qJ(4);
t69 = (-pkin(7) - t82) * t91;
t86 = t93 * pkin(7);
t70 = t82 * t93 + t86;
t36 = t69 * t95 + t70 * t97;
t140 = -t36 * qJD(5) - t146 * t73;
t35 = t69 * t97 - t70 * t95;
t139 = t35 * qJD(5) + t146 * t72;
t132 = pkin(1) * t98;
t133 = pkin(1) * t96;
t138 = (mrSges(3,1) * t133 + mrSges(3,2) * t132) * t90;
t74 = t90 * pkin(2) + t114;
t50 = t94 * t115 + t92 * t74;
t42 = qJ(4) * t90 + t50;
t33 = t91 * qJD(3) + t93 * t42;
t85 = t93 * qJD(3);
t101 = -(-t42 * t91 + t85) * t91 + t33 * t93;
t137 = m(5) * t101 + t141 * t90;
t49 = t74 * t94 - t78;
t103 = qJD(4) - t49;
t130 = pkin(4) * t93;
t112 = -pkin(3) - t130;
t34 = t112 * t90 + t103;
t58 = t72 * t90;
t59 = t73 * t90;
t136 = m(4) * t49 - m(5) * t103 - m(6) * t34 + mrSges(6,1) * t58 - mrSges(6,2) * t59 + (m(5) * pkin(3) + t147) * t90;
t134 = t59 / 0.2e1;
t131 = pkin(2) * t94;
t127 = Ifges(6,4) * t59;
t126 = t51 * mrSges(6,3);
t125 = t52 * mrSges(6,3);
t123 = t90 * mrSges(4,2);
t122 = t94 * t96;
t83 = pkin(2) + t132;
t121 = pkin(1) * t122 + t92 * t83;
t80 = t92 * t133;
t23 = t52 * mrSges(6,1) + t51 * mrSges(6,2);
t109 = t83 * t94 - t80;
t104 = -pkin(3) - t109;
t28 = t85 + (-pkin(7) * t90 - t42) * t91;
t29 = t90 * t86 + t33;
t7 = t28 * t97 - t29 * t95;
t8 = t28 * t95 + t29 * t97;
t61 = qJ(4) + t121;
t43 = (-pkin(7) - t61) * t91;
t44 = t61 * t93 + t86;
t14 = t43 * t97 - t44 * t95;
t15 = t43 * t95 + t44 * t97;
t65 = (t132 * t94 - t80) * qJD(2);
t100 = pkin(1) * (t92 * t98 + t122);
t63 = qJD(2) * t100;
t2 = t7 * qJD(5) + t72 * t45;
t26 = Ifges(6,2) * t58 + Ifges(6,6) * qJD(5) + t127;
t55 = Ifges(6,4) * t58;
t27 = Ifges(6,1) * t59 + Ifges(6,5) * qJD(5) + t55;
t3 = -t8 * qJD(5) - t73 * t45;
t56 = qJD(1) * t63;
t99 = -mrSges(3,1) * t106 - mrSges(3,2) * t105 - t57 * mrSges(4,2) + t27 * t145 + t26 * t144 + qJD(5) * (Ifges(6,5) * t67 - Ifges(6,6) * t68) / 0.2e1 + t34 * (mrSges(6,1) * t68 + mrSges(6,2) * t67) + mrSges(5,3) * t110 + (t58 * t144 - t72 * t52) * Ifges(6,2) + (t67 * t134 + t73 * t51) * Ifges(6,1) + (-mrSges(6,1) * t72 + mrSges(6,2) * t73 - t147) * t56 + (t2 * t72 - t3 * t73 - t67 * t7 - t8 * t68) * mrSges(6,3) + (-t68 * t134 + t58 * t145 - t142 + t143) * Ifges(6,4);
t75 = t112 - t131;
t60 = qJD(4) + t65;
t53 = t104 - t130;
t47 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t59;
t46 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t58;
t6 = -t15 * qJD(5) - t73 * t60;
t5 = t14 * qJD(5) + t72 * t60;
t1 = [m(4) * (-t56 * t109 + t57 * t121 + t50 * t65) + m(6) * (t14 * t3 + t15 * t2 + t5 * t8 + t53 * t56 + t6 * t7) - t65 * t123 - t14 * t126 - t15 * t125 + t5 * t46 + t6 * t47 + t53 * t23 + t99 + m(5) * (t56 * t104 + t61 * t110) + t137 * t60 - t138 * qJD(2) - t136 * t63; t75 * t23 - t35 * t126 - t36 * t125 + t99 + m(4) * (-t56 * t94 + t57 * t92) * pkin(2) + m(5) * (t56 * (-pkin(3) - t131) + t82 * t110 + t101 * qJD(4)) + t140 * t47 + t139 * t46 + t118 * t141 + (t139 * t8 + t140 * t7 + t2 * t36 + t3 * t35 + t56 * t75) * m(6) + (-m(4) * t50 + t123 - t137) * t64 + (t136 * t100 + t138) * qJD(1); m(6) * (t2 * t73 + t3 * t72 + t67 * t8 - t68 * t7) + t67 * t46 - t68 * t47 + (-t142 - t143) * mrSges(6,3); -t90 ^ 2 * t141 - t58 * t46 + t59 * t47 + t23 + (-t8 * t58 + t7 * t59 + t56) * m(6) + (-t101 * t90 + t56) * m(5); Ifges(6,5) * t51 - Ifges(6,6) * t52 - t2 * mrSges(6,2) + t3 * mrSges(6,1) - t34 * (mrSges(6,1) * t59 + mrSges(6,2) * t58) - t59 * (Ifges(6,1) * t58 - t127) / 0.2e1 + t26 * t134 - qJD(5) * (Ifges(6,5) * t58 - Ifges(6,6) * t59) / 0.2e1 - t7 * t46 + t8 * t47 + (t58 * t7 + t59 * t8) * mrSges(6,3) - (-Ifges(6,2) * t59 + t27 + t55) * t58 / 0.2e1;];
tauc = t1(:);
