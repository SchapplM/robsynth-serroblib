% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2020-01-03 12:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRPR2_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR2_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR2_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR2_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR2_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR2_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR2_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:07:22
% EndTime: 2020-01-03 12:07:25
% DurationCPUTime: 0.94s
% Computational Cost: add. (2094->171), mult. (4284->258), div. (0->0), fcn. (2304->8), ass. (0->96)
t114 = pkin(2) * qJD(3);
t74 = sin(pkin(9));
t77 = sin(qJ(3));
t123 = t74 * t77;
t113 = qJD(1) * pkin(1);
t78 = sin(qJ(2));
t80 = cos(qJ(3));
t119 = t78 * t80;
t81 = cos(qJ(2));
t93 = -t77 * t81 - t119;
t56 = t93 * t113;
t120 = t77 * t78;
t92 = t80 * t81 - t120;
t57 = t92 * t113;
t75 = cos(pkin(9));
t117 = t56 * t74 + t57 * t75 - (t75 * t80 - t123) * t114;
t76 = sin(qJ(5));
t121 = t76 * mrSges(6,3);
t73 = qJD(1) + qJD(2);
t72 = qJD(3) + t73;
t61 = qJD(5) * mrSges(6,1) - t72 * t121;
t79 = cos(qJ(5));
t124 = t72 * t79;
t62 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t124;
t94 = t76 * t61 - t79 * t62;
t104 = t78 * t113;
t63 = t73 * pkin(2) + t81 * t113;
t38 = -t104 * t77 + t80 * t63;
t31 = pkin(3) * t72 + t38;
t39 = t104 * t80 + t63 * t77;
t34 = t75 * t39;
t16 = t74 * t31 + t34;
t14 = pkin(8) * t72 + t16;
t11 = qJD(4) * t79 - t14 * t76;
t12 = qJD(4) * t76 + t14 * t79;
t97 = -t11 * t76 + t12 * t79;
t139 = m(6) * t97 - t94;
t122 = t75 * t77;
t118 = -t75 * t56 + t57 * t74 - (t74 * t80 + t122) * t114;
t137 = mrSges(3,1) * t78 + mrSges(3,2) * t81;
t70 = pkin(2) * t80 + pkin(3);
t136 = pkin(2) * t123 - t70 * t75;
t109 = qJD(3) * t80;
t110 = qJD(3) * t77;
t86 = (qJD(2) * t92 - t78 * t110) * pkin(1);
t21 = qJD(1) * t86 + t63 * t109;
t85 = (qJD(2) * t93 - t78 * t109) * pkin(1);
t22 = qJD(1) * t85 - t63 * t110;
t7 = t21 * t75 + t22 * t74;
t2 = qJD(5) * t11 + t7 * t79;
t135 = t2 * t79;
t132 = Ifges(6,1) * t76;
t131 = Ifges(6,4) * t76;
t129 = t11 * mrSges(6,3);
t126 = t39 * t74;
t71 = pkin(1) * t81 + pkin(2);
t101 = -pkin(1) * t120 + t80 * t71;
t55 = pkin(3) + t101;
t58 = pkin(1) * t119 + t71 * t77;
t116 = t74 * t55 + t75 * t58;
t115 = pkin(2) * t122 + t74 * t70;
t112 = Ifges(6,5) * qJD(5);
t111 = Ifges(6,6) * qJD(5);
t108 = qJD(5) * t12;
t107 = qJD(5) * t76;
t99 = -mrSges(6,1) * t79 + mrSges(6,2) * t76;
t98 = -t11 * t79 - t12 * t76;
t15 = t31 * t75 - t126;
t96 = t55 * t75 - t58 * t74;
t95 = -t61 * t79 - t62 * t76;
t91 = (Ifges(6,2) * t79 + t131) * t72;
t90 = t98 * mrSges(6,3);
t89 = (mrSges(6,1) * t76 + mrSges(6,2) * t79) * qJD(5);
t3 = -t7 * t76 - t108;
t84 = m(6) * (qJD(5) * t98 - t3 * t76 + t135);
t13 = -pkin(4) * t72 - t15;
t40 = t91 + t111;
t65 = Ifges(6,4) * t124;
t41 = t72 * t132 + t112 + t65;
t6 = t21 * t74 - t75 * t22;
t83 = t22 * mrSges(4,1) - t21 * mrSges(4,2) - t7 * mrSges(5,2) + mrSges(6,3) * t135 + t13 * t89 + qJD(5) ^ 2 * (Ifges(6,5) * t79 - Ifges(6,6) * t76) / 0.2e1 + t72 * (Ifges(6,1) * t79 - t131) * t107 + (-mrSges(5,1) + t99) * t6 - (t40 + t91) * t107 / 0.2e1 + (t41 + (0.3e1 * Ifges(6,4) * t79 - 0.2e1 * Ifges(6,2) * t76 + t132) * t72) * qJD(5) * t79 / 0.2e1;
t82 = -t3 * t121 + t83;
t69 = -pkin(3) * t75 - pkin(4);
t68 = pkin(3) * t74 + pkin(8);
t52 = pkin(8) + t115;
t51 = -pkin(4) + t136;
t50 = t99 * t72;
t42 = t72 * t89;
t29 = -t71 * t110 + t85;
t28 = t71 * t109 + t86;
t24 = pkin(8) + t116;
t23 = -pkin(4) - t96;
t18 = t38 * t75 - t126;
t17 = t38 * t74 + t34;
t9 = t28 * t74 - t75 * t29;
t1 = [m(6) * (t13 * t9 + t23 * t6) + t24 * t84 + (t24 * t95 + t90) * qJD(5) + m(4) * (t101 * t22 + t21 * t58 + t39 * t28 + t38 * t29) + m(5) * (t7 * t116 - t15 * t9 - t6 * t96) + (mrSges(4,1) * t29 - mrSges(5,1) * t9 - mrSges(4,2) * t28) * t72 + t23 * t42 + t9 * t50 + t82 + t137 * pkin(1) * qJD(2) * (-qJD(1) - t73) + (m(5) * t16 - t72 * mrSges(5,2) + t139) * (t28 * t75 + t29 * t74); (-qJD(5) * t52 * t62 + t117 * t61 + (-t3 - t108) * mrSges(6,3)) * t76 + (-t117 * t62 + (-t52 * t61 - t129) * qJD(5)) * t79 + t52 * t84 + t51 * t42 + t83 - t118 * t50 + (t117 * mrSges(5,2) + (-pkin(2) * t109 + t57) * mrSges(4,2) + t118 * mrSges(5,1) + (-pkin(2) * t110 - t56) * mrSges(4,1)) * t72 + t137 * t113 * (-qJD(2) + t73) + (-t117 * t97 - t118 * t13 + t51 * t6) * m(6) + (-t38 * t56 - t39 * t57 + (t39 * t109 - t38 * t110 + t21 * t77 + t22 * t80) * pkin(2)) * m(4) + (t7 * t115 - t117 * t16 + t118 * t15 + t136 * t6) * m(5); t68 * t84 + (t68 * t95 + t90) * qJD(5) + t94 * t18 + (mrSges(4,1) * t39 + mrSges(5,1) * t17 + mrSges(4,2) * t38 + mrSges(5,2) * t18) * t72 + t69 * t42 - t17 * t50 + t82 + (-t13 * t17 - t18 * t97 + t6 * t69) * m(6) + ((-t6 * t75 + t7 * t74) * pkin(3) + t15 * t17 - t16 * t18) * m(5); m(6) * (t2 * t76 + t3 * t79) + ((-t76 ^ 2 - t79 ^ 2) * t72 * mrSges(6,3) + t139) * qJD(5); t3 * mrSges(6,1) - t2 * mrSges(6,2) - t11 * t62 + t12 * t61 + ((t112 / 0.2e1 - t13 * mrSges(6,2) - t41 / 0.2e1 - t65 / 0.2e1 + t129) * t79 + (-t111 / 0.2e1 - t13 * mrSges(6,1) + t40 / 0.2e1 + t12 * mrSges(6,3) + (t131 / 0.2e1 + (-Ifges(6,1) / 0.2e1 + Ifges(6,2) / 0.2e1) * t79) * t72) * t76) * t72;];
tauc = t1(:);
