% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PRRRP2
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
% Datum: 2019-12-05 16:42
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRRRP2_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP2_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP2_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP2_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP2_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP2_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRP2_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:41:37
% EndTime: 2019-12-05 16:41:41
% DurationCPUTime: 0.82s
% Computational Cost: add. (841->165), mult. (1642->231), div. (0->0), fcn. (679->4), ass. (0->83)
t131 = mrSges(5,1) + mrSges(6,1);
t130 = mrSges(6,2) + mrSges(5,3);
t65 = sin(qJ(4));
t106 = pkin(2) * qJD(2);
t62 = qJD(2) + qJD(3);
t66 = sin(qJ(3));
t52 = pkin(7) * t62 + t66 * t106;
t67 = cos(qJ(4));
t29 = qJD(1) * t65 + t52 * t67;
t68 = cos(qJ(3));
t105 = pkin(2) * qJD(3);
t90 = qJD(2) * t105;
t86 = t68 * t90;
t9 = t29 * qJD(4) + t65 * t86;
t121 = t9 * t65;
t103 = qJD(4) * t65;
t102 = qJD(4) * t67;
t108 = qJD(1) * t102 + t67 * t86;
t8 = -t52 * t103 + t108;
t129 = t67 * t8 + t121;
t97 = qJD(4) * qJ(5);
t22 = t29 + t97;
t101 = t22 * qJD(4);
t115 = t52 * t65;
t28 = qJD(1) * t67 - t115;
t126 = -t28 + qJD(5);
t16 = -qJD(4) * pkin(4) + t126;
t6 = (qJD(5) - t115) * qJD(4) + t108;
t128 = -t65 * t101 + t16 * t102 + t6 * t67 + t121;
t113 = t62 * t67;
t95 = mrSges(6,2) * t113;
t51 = qJD(4) * mrSges(6,3) + t95;
t109 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t113 + t51;
t114 = t62 * t65;
t110 = t131 * qJD(4) - t130 * t114;
t125 = t109 * t67 - t110 * t65;
t124 = pkin(2) * t68;
t120 = t9 * t67;
t119 = Ifges(5,4) * t65;
t117 = Ifges(6,5) * t65;
t116 = Ifges(6,5) * t67;
t112 = t65 * t68;
t111 = t67 * t68;
t104 = qJD(3) * t66;
t100 = t65 * qJD(5);
t99 = -qJD(2) - t62;
t98 = -qJD(3) + t62;
t94 = t68 * t106;
t91 = t114 / 0.2e1;
t88 = t103 / 0.2e1;
t87 = t102 / 0.2e1;
t85 = t66 * t90;
t84 = t62 * t88;
t82 = -mrSges(5,1) * t67 + mrSges(5,2) * t65;
t81 = -mrSges(6,1) * t67 - mrSges(6,3) * t65;
t80 = pkin(4) * t65 - qJ(5) * t67;
t54 = -t67 * pkin(4) - t65 * qJ(5) - pkin(3);
t79 = (Ifges(6,1) * t65 - t116) * t62;
t78 = (Ifges(5,2) * t67 + t119) * t62;
t77 = (mrSges(5,1) * t65 + mrSges(5,2) * t67) * qJD(4);
t76 = (mrSges(6,1) * t65 - mrSges(6,3) * t67) * qJD(4);
t40 = pkin(4) * t103 - t67 * t97 - t100;
t71 = m(5) * (-t28 * t65 + t29 * t67) + m(6) * (t16 * t65 + t22 * t67) + t125;
t70 = (-t109 * t65 - t110 * t67) * qJD(4) + m(6) * t128 + m(5) * (-t28 * t102 - t29 * t103 + t129);
t11 = t54 * t62 - t94;
t56 = Ifges(6,5) * t114;
t30 = Ifges(6,6) * qJD(4) - Ifges(6,3) * t113 + t56;
t31 = Ifges(5,6) * qJD(4) + t78;
t32 = Ifges(6,4) * qJD(4) + t79;
t57 = Ifges(5,4) * t113;
t33 = Ifges(5,1) * t114 + Ifges(5,5) * qJD(4) + t57;
t53 = -t62 * pkin(3) - t94;
t7 = t85 + (t80 * qJD(4) - t100) * t62;
t69 = t53 * t77 + t11 * t76 + t7 * t81 + t82 * t85 + (-Ifges(6,3) * t67 + t117) * t84 + t30 * t88 + ((Ifges(6,4) + Ifges(5,5)) * t67 + (-Ifges(5,6) + Ifges(6,6)) * t65) * qJD(4) ^ 2 / 0.2e1 - (t78 + t31) * t103 / 0.2e1 + 0.2e1 * (t117 - t119 + (Ifges(5,1) + Ifges(6,1)) * t67) * t84 + ((-t28 * t67 - t29 * t65) * qJD(4) + t129) * mrSges(5,3) + (t79 + t33 + t32) * t87 + t128 * mrSges(6,2) + (-(Ifges(6,3) * t65 + t116) * t102 + (0.3e1 * Ifges(5,4) * t67 + (Ifges(5,1) - 0.2e1 * Ifges(5,2)) * t65) * t87) * t62;
t60 = -pkin(3) - t124;
t47 = t54 - t124;
t43 = t80 * t62;
t42 = t82 * t62;
t41 = t81 * t62;
t35 = t62 * t77;
t34 = t62 * t76;
t23 = pkin(2) * t104 + t40;
t1 = [m(5) * (t8 * t65 - t120) + m(6) * (t6 * t65 - t120) + (t71 + t130 * t62 * (-t65 ^ 2 - t67 ^ 2)) * qJD(4); t60 * t35 + t23 * t41 + t47 * t34 + t69 + ((m(5) * (qJD(2) * t60 + t53) + t42 + t99 * mrSges(4,1)) * t66 + (t99 * mrSges(4,2) + t71) * t68) * t105 + m(6) * (t11 * t23 + t7 * t47) + t70 * (pkin(2) * t66 + pkin(7)); t54 * t34 - pkin(3) * t35 + t40 * t41 + t69 + ((t98 * mrSges(4,1) - t41 - t42) * t66 + (t98 * mrSges(4,2) - t125) * t68 - m(6) * (t11 * t66 + t22 * t111 + t16 * t112) + (-pkin(3) * t104 - t29 * t111 + t28 * t112 - t53 * t66) * m(5)) * t106 + t70 * pkin(7) + m(6) * (t11 * t40 + t7 * t54); -t8 * mrSges(5,2) + t6 * mrSges(6,3) + qJD(5) * t51 - t43 * t41 - t131 * t9 + t110 * t29 - t109 * t28 + ((-t32 / 0.2e1 - t33 / 0.2e1 - t57 / 0.2e1 - t16 * mrSges(6,2) + t28 * mrSges(5,3) - t53 * mrSges(5,2) + t11 * mrSges(6,3) + Ifges(6,5) * t113 / 0.2e1 + (Ifges(6,4) / 0.2e1 + Ifges(5,5) / 0.2e1 - pkin(4) * mrSges(6,2)) * qJD(4)) * t67 + (-t56 / 0.2e1 + t22 * mrSges(6,2) + t29 * mrSges(5,3) - t53 * mrSges(5,1) - t11 * mrSges(6,1) - t30 / 0.2e1 + t31 / 0.2e1 + Ifges(5,4) * t91 + (Ifges(6,6) / 0.2e1 - Ifges(5,6) / 0.2e1 - qJ(5) * mrSges(6,2)) * qJD(4) + (Ifges(6,3) / 0.2e1 - Ifges(5,1) / 0.2e1 - Ifges(6,1) / 0.2e1 + Ifges(5,2) / 0.2e1) * t113) * t65) * t62 + (-t9 * pkin(4) + t6 * qJ(5) - t11 * t43 + t126 * t22 - t16 * t29) * m(6); t41 * t114 + (-t51 + t95) * qJD(4) + 0.2e1 * (t9 / 0.2e1 + t11 * t91 - t101 / 0.2e1) * m(6);];
tauc = t1(:);
