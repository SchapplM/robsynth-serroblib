% Calculate joint inertia matrix for
% S6RRPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRP1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP1_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP1_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP1_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP1_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP1_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:39:01
% EndTime: 2019-03-09 11:39:04
% DurationCPUTime: 1.49s
% Computational Cost: add. (2136->247), mult. (4006->336), div. (0->0), fcn. (4493->8), ass. (0->96)
t167 = Ifges(6,4) + Ifges(7,4);
t162 = Ifges(6,6) + Ifges(7,6);
t102 = sin(qJ(4));
t143 = cos(qJ(4));
t100 = cos(pkin(10));
t103 = sin(qJ(2));
t105 = cos(qJ(2));
t99 = sin(pkin(10));
t64 = t100 * t105 - t103 * t99;
t65 = t100 * t103 + t105 * t99;
t49 = t102 * t65 - t143 * t64;
t166 = t162 * t49;
t165 = Ifges(6,1) + Ifges(7,1);
t164 = Ifges(6,5) + Ifges(7,5);
t163 = Ifges(6,2) + Ifges(7,2);
t104 = cos(qJ(5));
t161 = t167 * t104;
t101 = sin(qJ(5));
t160 = t167 * t101;
t159 = Ifges(6,3) + Ifges(7,3);
t50 = t102 * t64 + t143 * t65;
t158 = (-t163 * t101 + t161) * t50 + t166;
t157 = (t165 * t104 - t160) * t50 + t164 * t49;
t156 = t163 * t104 + t160;
t155 = t165 * t101 + t161;
t116 = t164 * t101 + t162 * t104;
t133 = t101 ^ 2 + t104 ^ 2;
t154 = 0.2e1 * t133;
t136 = -qJ(3) - pkin(7);
t75 = t136 * t103;
t77 = t136 * t105;
t51 = t100 * t75 + t77 * t99;
t110 = -pkin(8) * t65 + t51;
t52 = -t100 * t77 + t99 * t75;
t39 = pkin(8) * t64 + t52;
t23 = t102 * t39 - t143 * t110;
t153 = t23 ^ 2;
t152 = 0.2e1 * t23;
t88 = -pkin(2) * t105 - pkin(1);
t56 = -pkin(3) * t64 + t88;
t151 = 0.2e1 * t56;
t150 = 0.2e1 * t64;
t90 = t101 * mrSges(7,2);
t73 = -t104 * mrSges(7,1) + t90;
t149 = 0.2e1 * t73;
t148 = m(7) * pkin(5);
t146 = pkin(2) * t99;
t142 = pkin(9) * t104;
t141 = t104 * pkin(5);
t22 = pkin(4) * t49 - pkin(9) * t50 + t56;
t25 = t102 * t110 + t143 * t39;
t6 = t101 * t22 + t104 * t25;
t140 = t104 * t6;
t86 = pkin(2) * t100 + pkin(3);
t60 = -t102 * t146 + t143 * t86;
t139 = t60 * mrSges(5,1);
t61 = t102 * t86 + t143 * t146;
t138 = t61 * mrSges(5,2);
t124 = t104 * t50;
t125 = t101 * t50;
t26 = mrSges(7,1) * t125 + mrSges(7,2) * t124;
t132 = t103 ^ 2 + t105 ^ 2;
t131 = mrSges(6,2) * t101;
t126 = t101 * mrSges(7,3);
t59 = pkin(9) + t61;
t123 = t104 * t59;
t89 = t104 * qJ(6);
t122 = 0.2e1 * mrSges(7,3);
t119 = t133 * t59;
t118 = -t64 * mrSges(4,1) + t65 * mrSges(4,2);
t5 = -t101 * t25 + t104 * t22;
t117 = t164 * t124 + t159 * t49;
t115 = t155 * t101 + t156 * t104 + Ifges(5,3);
t1 = pkin(5) * t49 - t50 * t89 + t5;
t114 = -t5 * mrSges(6,3) - t1 * mrSges(7,3);
t113 = -t101 * t5 + t140;
t112 = mrSges(6,3) * t154;
t111 = mrSges(6,1) * t101 + mrSges(6,2) * t104;
t58 = -pkin(4) - t60;
t3 = -qJ(6) * t125 + t6;
t74 = -mrSges(6,1) * t104 + t131;
t8 = pkin(5) * t125 + t23;
t109 = -t25 * mrSges(5,2) + mrSges(6,3) * t140 + Ifges(5,5) * t50 + t8 * t73 + (-mrSges(5,1) + t74) * t23 + t157 * t101 / 0.2e1 - t156 * t125 / 0.2e1 + t155 * t124 / 0.2e1 + (-Ifges(5,6) + t116 / 0.2e1) * t49 + (t3 * mrSges(7,3) + t158 / 0.2e1) * t104;
t87 = -pkin(4) - t141;
t76 = t89 + t142;
t72 = (-qJ(6) - pkin(9)) * t101;
t55 = t58 - t141;
t54 = t89 + t123;
t53 = (-qJ(6) - t59) * t101;
t44 = t50 * mrSges(5,2);
t31 = mrSges(6,1) * t49 - mrSges(6,3) * t124;
t30 = mrSges(7,1) * t49 - mrSges(7,3) * t124;
t29 = -mrSges(6,2) * t49 - mrSges(6,3) * t125;
t28 = -mrSges(7,2) * t49 - mrSges(7,3) * t125;
t27 = t111 * t50;
t2 = [-0.2e1 * pkin(1) * (-t105 * mrSges(3,1) + t103 * mrSges(3,2)) + t103 * (Ifges(3,1) * t103 + Ifges(3,4) * t105) + t105 * (Ifges(3,4) * t103 + Ifges(3,2) * t105) + t52 * mrSges(4,3) * t150 + 0.2e1 * t88 * t118 + Ifges(4,2) * t64 ^ 2 + t44 * t151 + 0.2e1 * t8 * t26 + t27 * t152 + 0.2e1 * t3 * t28 + 0.2e1 * t6 * t29 + 0.2e1 * t1 * t30 + 0.2e1 * t5 * t31 + Ifges(2,3) + 0.2e1 * t132 * pkin(7) * mrSges(3,3) + (-0.2e1 * t51 * mrSges(4,3) + Ifges(4,1) * t65 + Ifges(4,4) * t150) * t65 + (mrSges(5,1) * t151 - 0.2e1 * t25 * mrSges(5,3) + Ifges(5,2) * t49 + t117) * t49 + (mrSges(5,3) * t152 + Ifges(5,1) * t50 - 0.2e1 * Ifges(5,4) * t49 + t157 * t104 + (-t158 - t166) * t101) * t50 + m(3) * (t132 * pkin(7) ^ 2 + pkin(1) ^ 2) + m(4) * (t51 ^ 2 + t52 ^ 2 + t88 ^ 2) + m(5) * (t25 ^ 2 + t56 ^ 2 + t153) + m(6) * (t5 ^ 2 + t6 ^ 2 + t153) + m(7) * (t1 ^ 2 + t3 ^ 2 + t8 ^ 2); t109 + (-t59 * t31 + t114) * t101 + m(5) * (-t23 * t60 + t25 * t61) + m(7) * (t1 * t53 + t3 * t54 + t55 * t8) + (-t103 * mrSges(3,1) - t105 * mrSges(3,2)) * pkin(7) + (-t61 * t49 - t60 * t50) * mrSges(5,3) + Ifges(3,6) * t105 + Ifges(3,5) * t103 + t29 * t123 + t58 * t27 + Ifges(4,6) * t64 + Ifges(4,5) * t65 + t51 * mrSges(4,1) - t52 * mrSges(4,2) + t53 * t30 + t54 * t28 + t55 * t26 + m(6) * (t113 * t59 + t23 * t58) + ((-t100 * t65 + t64 * t99) * mrSges(4,3) + m(4) * (t100 * t51 + t52 * t99)) * pkin(2); 0.2e1 * t139 - 0.2e1 * t138 + t55 * t149 + 0.2e1 * t58 * t74 + Ifges(3,3) + Ifges(4,3) + (-t53 * t101 + t54 * t104) * t122 + t59 * t112 + m(7) * (t53 ^ 2 + t54 ^ 2 + t55 ^ 2) + m(6) * (t133 * t59 ^ 2 + t58 ^ 2) + m(5) * (t60 ^ 2 + t61 ^ 2) + t115 + (0.2e1 * mrSges(4,1) * t100 - 0.2e1 * mrSges(4,2) * t99 + m(4) * (t100 ^ 2 + t99 ^ 2) * pkin(2)) * pkin(2); t49 * mrSges(5,1) + t44 + (t30 + t31) * t104 + (t28 + t29) * t101 + m(6) * (t101 * t6 + t104 * t5) + m(7) * (t1 * t104 + t101 * t3) + m(5) * t56 + m(4) * t88 + t118; m(7) * (t101 * t54 + t104 * t53); m(4) + m(5) + (m(6) / 0.2e1 + m(7) / 0.2e1) * t154; t109 + (-pkin(9) * t31 + t114) * t101 + m(7) * (t1 * t72 + t3 * t76 + t8 * t87) + t29 * t142 + t87 * t26 + t72 * t30 + t76 * t28 - pkin(4) * t27 + m(6) * (-pkin(4) * t23 + pkin(9) * t113); t139 - t138 + (t58 - pkin(4)) * t74 + (t55 + t87) * t73 + m(7) * (t53 * t72 + t54 * t76 + t55 * t87) + m(6) * (-pkin(4) * t58 + pkin(9) * t119) + ((t54 + t76) * t104 + (-t53 - t72) * t101) * mrSges(7,3) + (pkin(9) * t133 + t119) * mrSges(6,3) + t115; m(7) * (t101 * t76 + t104 * t72); -0.2e1 * pkin(4) * t74 + t87 * t149 + (-t72 * t101 + t76 * t104) * t122 + pkin(9) * t112 + m(7) * (t72 ^ 2 + t76 ^ 2 + t87 ^ 2) + m(6) * (pkin(9) ^ 2 * t133 + pkin(4) ^ 2) + t115; mrSges(6,1) * t5 + mrSges(7,1) * t1 - mrSges(6,2) * t6 - mrSges(7,2) * t3 - t162 * t125 + (m(7) * t1 + t30) * pkin(5) + t117; mrSges(7,1) * t53 - mrSges(7,2) * t54 - t111 * t59 + (m(7) * t53 - t126) * pkin(5) + t116; -t131 - t90 + (mrSges(6,1) + mrSges(7,1) + t148) * t104; mrSges(7,1) * t72 - mrSges(7,2) * t76 - t111 * pkin(9) + (m(7) * t72 - t126) * pkin(5) + t116; (0.2e1 * mrSges(7,1) + t148) * pkin(5) + t159; m(7) * t8 + t26; m(7) * t55 + t73; 0; m(7) * t87 + t73; 0; m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t2(1) t2(2) t2(4) t2(7) t2(11) t2(16); t2(2) t2(3) t2(5) t2(8) t2(12) t2(17); t2(4) t2(5) t2(6) t2(9) t2(13) t2(18); t2(7) t2(8) t2(9) t2(10) t2(14) t2(19); t2(11) t2(12) t2(13) t2(14) t2(15) t2(20); t2(16) t2(17) t2(18) t2(19) t2(20) t2(21);];
Mq  = res;
