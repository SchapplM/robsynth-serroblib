% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRRRP1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP1_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP1_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP1_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP1_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRP1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:39:48
% EndTime: 2019-12-05 16:39:53
% DurationCPUTime: 1.16s
% Computational Cost: add. (879->198), mult. (1677->270), div. (0->0), fcn. (713->4), ass. (0->93)
t139 = Ifges(5,1) + Ifges(6,1);
t118 = pkin(2) * qJD(2);
t72 = qJD(2) + qJD(3);
t76 = sin(qJ(3));
t55 = pkin(7) * t72 + t118 * t76;
t101 = qJ(5) * t72 + t55;
t75 = sin(qJ(4));
t115 = qJD(1) * t75;
t77 = cos(qJ(4));
t12 = t101 * t77 + t115;
t128 = t72 * t75;
t51 = qJD(4) * mrSges(6,1) - mrSges(6,3) * t128;
t52 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t128;
t127 = t72 * t77;
t53 = -qJD(4) * mrSges(6,2) + mrSges(6,3) * t127;
t125 = t77 * mrSges(5,3);
t54 = -qJD(4) * mrSges(5,2) + t125 * t72;
t135 = (t53 + t54) * t77 - (t51 + t52) * t75;
t78 = cos(qJ(3));
t134 = pkin(2) * t78;
t29 = t55 * t77 + t115;
t117 = pkin(2) * qJD(3);
t104 = qJD(2) * t117;
t97 = t78 * t104;
t7 = -qJD(4) * t29 - t75 * t97;
t133 = t7 * t75;
t132 = Ifges(5,4) * t75;
t130 = Ifges(6,4) * t75;
t126 = t75 * t78;
t124 = t77 * t78;
t123 = -qJ(5) - pkin(7);
t70 = t77 * qJD(1);
t120 = qJD(4) * t70 + t77 * t97;
t113 = qJD(4) * t75;
t105 = t72 * t113;
t112 = qJD(4) * t77;
t34 = t72 * mrSges(6,2) * t112 + mrSges(6,1) * t105;
t66 = pkin(2) * t76 + pkin(7);
t116 = -qJ(5) - t66;
t114 = qJD(3) * t76;
t69 = t77 * qJD(5);
t111 = -qJD(2) - t72;
t110 = -qJD(3) + t72;
t108 = t78 * t118;
t107 = t78 * t117;
t106 = pkin(4) * t113;
t68 = -t77 * pkin(4) - pkin(3);
t100 = qJD(4) * t123;
t99 = qJD(4) * t116;
t98 = t76 * t104;
t6 = -t113 * t55 + t120;
t94 = t6 * t77 - t133;
t89 = t101 * t75;
t11 = t70 - t89;
t8 = qJD(4) * pkin(4) + t11;
t93 = -t12 * t75 - t8 * t77;
t92 = -mrSges(5,1) * t77 + mrSges(5,2) * t75;
t91 = -mrSges(6,1) * t77 + mrSges(6,2) * t75;
t28 = -t55 * t75 + t70;
t90 = -t28 * t77 - t29 * t75;
t88 = m(6) * (t12 * t77 - t75 * t8);
t87 = m(5) * (-t28 * t75 + t29 * t77);
t86 = (Ifges(5,2) * t77 + t132) * t72;
t85 = (Ifges(6,2) * t77 + t130) * t72;
t84 = (mrSges(5,1) * t75 + mrSges(5,2) * t77) * qJD(4);
t2 = -qJD(4) * t89 + t69 * t72 + t120;
t23 = t68 * t72 + qJD(5) - t108;
t30 = Ifges(6,6) * qJD(4) + t85;
t31 = Ifges(5,6) * qJD(4) + t86;
t63 = Ifges(6,4) * t127;
t32 = Ifges(6,1) * t128 + Ifges(6,5) * qJD(4) + t63;
t64 = Ifges(5,4) * t127;
t33 = Ifges(5,1) * t128 + Ifges(5,5) * qJD(4) + t64;
t42 = pkin(4) * t105 + t98;
t56 = -t72 * pkin(3) - t108;
t79 = t2 * t77 * mrSges(6,3) + t6 * t125 + t42 * t91 + t56 * t84 + t92 * t98 + t23 * (mrSges(6,1) * t75 + mrSges(6,2) * t77) * qJD(4) + (t139 * t77 - t130 - t132) * t105 + ((Ifges(5,5) + Ifges(6,5)) * t77 + (-Ifges(5,6) - Ifges(6,6)) * t75) * qJD(4) ^ 2 / 0.2e1 - (t30 + t31 + t85 + t86) * t113 / 0.2e1 + (t32 + t33 + ((t139 - 0.2e1 * Ifges(5,2) - 0.2e1 * Ifges(6,2)) * t75 + 0.3e1 * (Ifges(5,4) + Ifges(6,4)) * t77) * t72) * t112 / 0.2e1;
t71 = t77 * qJ(5);
t67 = -pkin(3) - t134;
t60 = pkin(7) * t77 + t71;
t59 = t123 * t75;
t57 = t68 - t134;
t50 = pkin(2) * t114 + t106;
t49 = t66 * t77 + t71;
t48 = t116 * t75;
t44 = t92 * t72;
t43 = t91 * t72;
t41 = -t75 * qJD(5) + t100 * t77;
t40 = t100 * t75 + t69;
t35 = t72 * t84;
t10 = (-qJD(5) - t107) * t75 + t77 * t99;
t9 = t107 * t77 + t75 * t99 + t69;
t3 = (-qJD(5) * t72 - t97) * t75 - t12 * qJD(4);
t1 = [m(5) * (t6 * t75 + t7 * t77) + m(6) * (t2 * t75 + t3 * t77) + (t87 + t88 + (mrSges(5,3) + mrSges(6,3)) * t72 * (-t75 ^ 2 - t77 ^ 2) + t135) * qJD(4); (-t7 * mrSges(5,3) - t3 * mrSges(6,3)) * t75 + m(5) * t94 * t66 + m(6) * (t10 * t8 + t12 * t9 + t2 * t49 + t23 * t50 + t3 * t48 + t42 * t57) + t79 + (t90 * mrSges(5,3) + (m(5) * t90 - t77 * t52 - t75 * t54) * t66 + ((-t48 * t77 - t49 * t75) * t72 + t93) * mrSges(6,3)) * qJD(4) + t67 * t35 + t50 * t43 + t10 * t51 + t9 * t53 + t57 * t34 + ((m(5) * (qJD(2) * t67 + t56) + t44 + t111 * mrSges(4,1)) * t76 + (mrSges(4,2) * t111 - t75 * t52 + t77 * t54 + t87) * t78) * t117; (-t3 * t75 + ((-t59 * t77 - t60 * t75) * t72 + t93) * qJD(4)) * mrSges(6,3) + (qJD(4) * t90 - t133) * mrSges(5,3) + m(6) * (t23 * t106 + t12 * t40 + t2 * t60 + t3 * t59 + t41 * t8 + t42 * t68) + t79 + (m(5) * (-t112 * t28 - t113 * t29 + t94) - t54 * t113 - t52 * t112) * pkin(7) + t68 * t34 + t41 * t51 + t40 * t53 - pkin(3) * t35 + t43 * t106 + ((mrSges(4,1) * t110 - t43 - t44) * t76 + (mrSges(4,2) * t110 - t135) * t78 - m(6) * (t12 * t124 - t126 * t8 + t23 * t76) + (-pkin(3) * t114 - t124 * t29 + t126 * t28 - t56 * t76) * m(5)) * t118; t7 * mrSges(5,1) - t6 * mrSges(5,2) - t2 * mrSges(6,2) - t11 * t53 - t28 * t54 + t29 * t52 + (m(6) * pkin(4) + mrSges(6,1)) * t3 + (t51 - m(6) * (t11 - t8)) * t12 + ((-t32 / 0.2e1 - t33 / 0.2e1 - t63 / 0.2e1 - t64 / 0.2e1 - t56 * mrSges(5,2) - t23 * mrSges(6,2) + t28 * mrSges(5,3) + t8 * mrSges(6,3) + (Ifges(6,5) / 0.2e1 + Ifges(5,5) / 0.2e1 - pkin(4) * mrSges(6,3)) * qJD(4)) * t77 + (-t56 * mrSges(5,1) - t23 * mrSges(6,1) + t29 * mrSges(5,3) + t12 * mrSges(6,3) + t30 / 0.2e1 + t31 / 0.2e1 + (Ifges(6,4) / 0.2e1 + Ifges(5,4) / 0.2e1) * t128 + (-Ifges(6,6) / 0.2e1 - Ifges(5,6) / 0.2e1) * qJD(4) + (-m(6) * t23 - t43) * pkin(4) + (-Ifges(6,1) / 0.2e1 - Ifges(5,1) / 0.2e1 + Ifges(6,2) / 0.2e1 + Ifges(5,2) / 0.2e1) * t127) * t75) * t72; m(6) * t42 + (t75 * t51 - t77 * t53 - t88) * t72 + t34;];
tauc = t1(:);
