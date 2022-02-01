% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% m [6x1]
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
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 10:34
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRPRR3_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR3_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR3_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR3_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR3_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR3_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR3_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:34:02
% EndTime: 2022-01-20 10:34:05
% DurationCPUTime: 0.87s
% Computational Cost: add. (2582->135), mult. (5384->187), div. (0->0), fcn. (4320->8), ass. (0->93)
t71 = sin(qJ(5));
t74 = cos(qJ(5));
t138 = -Ifges(6,4) * t71 + (Ifges(6,1) - Ifges(6,2)) * t74;
t67 = t71 ^ 2;
t68 = t74 ^ 2;
t112 = t67 + t68;
t132 = t112 * mrSges(6,3);
t135 = -mrSges(5,2) + t132;
t61 = -mrSges(6,1) * t74 + mrSges(6,2) * t71;
t133 = -mrSges(5,1) + t61;
t70 = cos(pkin(9));
t73 = sin(qJ(2));
t117 = t70 * t73;
t69 = sin(pkin(9));
t76 = cos(qJ(2));
t54 = (-t69 * t76 - t117) * pkin(1);
t118 = t69 * t73;
t55 = (t70 * t76 - t118) * pkin(1);
t72 = sin(qJ(4));
t75 = cos(qJ(4));
t38 = t54 * t72 + t55 * t75;
t131 = t135 * t38 + t54 * mrSges(4,1) - t55 * mrSges(4,2) + (-t73 * mrSges(3,1) - t76 * mrSges(3,2)) * pkin(1);
t130 = m(6) / 0.2e1;
t129 = pkin(4) / 0.2e1;
t64 = pkin(1) * t76 + pkin(2);
t98 = -pkin(1) * t118 + t70 * t64;
t48 = pkin(3) + t98;
t52 = pkin(1) * t117 + t64 * t69;
t33 = t48 * t75 - t52 * t72;
t27 = -pkin(4) - t33;
t128 = -t27 / 0.2e1;
t127 = -t33 / 0.2e1;
t126 = -t38 / 0.2e1;
t122 = pkin(2) * t69;
t63 = pkin(2) * t70 + pkin(3);
t51 = -t72 * t122 + t63 * t75;
t49 = -pkin(4) - t51;
t125 = -t49 / 0.2e1;
t124 = -t51 / 0.2e1;
t37 = -t75 * t54 + t55 * t72;
t120 = t37 * mrSges(5,1);
t119 = t37 * t61;
t116 = t71 * mrSges(6,1);
t114 = t74 * mrSges(6,2);
t103 = t112 * t38;
t34 = t48 * t72 + t52 * t75;
t28 = pkin(8) + t34;
t3 = t133 * t37 + m(6) * (t28 * t103 + t27 * t37) + m(5) * (-t33 * t37 + t34 * t38) + m(4) * (t52 * t55 + t98 * t54) + t131;
t111 = t3 * qJD(1);
t104 = t112 * t33;
t79 = -t33 * mrSges(5,2) + mrSges(6,3) * t104 + t133 * t34;
t4 = m(6) * (t28 * t104 + t27 * t34) + t79;
t110 = t4 * qJD(1);
t90 = t114 + t116;
t85 = t27 * t90;
t66 = Ifges(6,4) * t74;
t99 = t138 * t71 + t74 * t66;
t15 = t85 + t99;
t109 = t15 * qJD(1);
t107 = t127 + t124;
t53 = t75 * t122 + t63 * t72;
t106 = t34 / 0.2e1 + t53 / 0.2e1;
t50 = pkin(8) + t53;
t102 = t112 * t50;
t101 = t112 * t51;
t100 = Ifges(6,5) * t74 - Ifges(6,6) * t71;
t97 = t127 + t128 + t129;
t96 = t124 + t129 + t125;
t93 = t126 + t128 + t125;
t86 = pkin(4) * t90;
t92 = -t86 / 0.2e1 + t99;
t77 = (t28 * t101 + t33 * t102 + t27 * t53 + t34 * t49) * t130;
t80 = m(6) * (-pkin(4) * t37 + pkin(8) * t103);
t2 = t77 - t80 / 0.2e1 + t133 * (-t37 / 0.2e1 + t106) + t135 * (t126 - t107);
t78 = -t51 * mrSges(5,2) + mrSges(6,3) * t101 + t133 * t53;
t6 = m(6) * (t50 * t101 + t49 * t53) + t78;
t89 = t2 * qJD(1) + t6 * qJD(2);
t84 = t49 * t90;
t16 = t84 + t99;
t8 = (t93 * mrSges(6,2) - t66) * t74 + (t93 * mrSges(6,1) - t138) * t71;
t88 = -t8 * qJD(1) + t16 * qJD(2);
t87 = -t114 / 0.2e1 - t116 / 0.2e1;
t10 = (t97 * mrSges(6,2) - t66) * t74 + (t97 * mrSges(6,1) - t138) * t71;
t12 = (t96 * mrSges(6,2) - t66) * t74 + (t96 * mrSges(6,1) - t138) * t71;
t17 = -t86 + t99;
t81 = t10 * qJD(1) + t12 * qJD(2) - t17 * qJD(4);
t39 = t84 / 0.2e1;
t18 = t85 / 0.2e1;
t13 = t87 * t51 + t39 + t92;
t11 = t87 * t33 + t18 + t92;
t9 = t87 * t38 + t18 + t39 + t99;
t1 = t77 + t80 / 0.2e1 + t119 / 0.2e1 - t120 / 0.2e1 + t38 * t132 / 0.2e1 + t133 * t106 + (t33 + t51) * mrSges(6,3) * (t68 / 0.2e1 + t67 / 0.2e1) + (t126 + t107) * mrSges(5,2);
t5 = [qJD(2) * t3 + qJD(4) * t4 + qJD(5) * t15, t1 * qJD(4) + t9 * qJD(5) + t111 + (t119 - t120 + 0.2e1 * (t38 * t102 + t37 * t49) * t130 + m(5) * (-t37 * t51 + t38 * t53) + m(4) * (t54 * t70 + t55 * t69) * pkin(2) + t131) * qJD(2), 0, t110 + t1 * qJD(2) + (m(6) * (-pkin(4) * t34 + pkin(8) * t104) + t79) * qJD(4) + t11 * qJD(5), t109 + t9 * qJD(2) + t11 * qJD(4) + (t61 * t28 + t100) * qJD(5); qJD(4) * t2 - qJD(5) * t8 - t111, qJD(4) * t6 + qJD(5) * t16, 0, (m(6) * (-pkin(4) * t53 + pkin(8) * t101) + t78) * qJD(4) + t13 * qJD(5) + t89, t13 * qJD(4) + (t61 * t50 + t100) * qJD(5) + t88; 0, 0, 0, 0, -t90 * qJD(5); -qJD(2) * t2 - qJD(5) * t10 - t110, -qJD(5) * t12 - t89, 0, t17 * qJD(5), (t61 * pkin(8) + t100) * qJD(5) - t81; qJD(2) * t8 + qJD(4) * t10 - t109, qJD(4) * t12 - t88, 0, t81, 0;];
Cq = t5;
