% Calculate joint inertia matrix for
% S6RRPRRP2
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
% Datum: 2019-03-09 11:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRP2_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP2_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP2_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP2_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP2_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP2_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:43:07
% EndTime: 2019-03-09 11:43:10
% DurationCPUTime: 1.13s
% Computational Cost: add. (2150->232), mult. (4037->316), div. (0->0), fcn. (4517->8), ass. (0->97)
t167 = Ifges(6,1) + Ifges(7,1);
t166 = Ifges(7,4) + Ifges(6,5);
t165 = Ifges(7,2) + Ifges(6,3);
t104 = cos(qJ(5));
t101 = sin(qJ(5));
t132 = Ifges(7,5) * t101;
t134 = Ifges(6,4) * t101;
t102 = sin(qJ(4));
t151 = cos(qJ(4));
t100 = cos(pkin(10));
t103 = sin(qJ(2));
t105 = cos(qJ(2));
t99 = sin(pkin(10));
t66 = t100 * t105 - t103 * t99;
t67 = t100 * t103 + t105 * t99;
t49 = t102 * t67 - t151 * t66;
t50 = t102 * t66 + t151 * t67;
t164 = (t167 * t104 + t132 - t134) * t50 + t166 * t49;
t131 = Ifges(7,5) * t104;
t133 = Ifges(6,4) * t104;
t163 = t167 * t101 - t131 + t133;
t136 = t101 ^ 2 + t104 ^ 2;
t162 = (Ifges(6,6) - Ifges(7,6)) * t104 + t166 * t101;
t87 = -pkin(2) * t105 - pkin(1);
t58 = -pkin(3) * t66 + t87;
t22 = pkin(4) * t49 - pkin(9) * t50 + t58;
t143 = -qJ(3) - pkin(7);
t76 = t143 * t103;
t77 = t143 * t105;
t52 = t100 * t76 + t77 * t99;
t115 = -pkin(8) * t67 + t52;
t53 = -t100 * t77 + t99 * t76;
t39 = pkin(8) * t66 + t53;
t25 = t102 * t115 + t151 * t39;
t6 = -t101 * t25 + t104 * t22;
t7 = t101 * t22 + t104 * t25;
t121 = -t101 * t6 + t104 * t7;
t161 = (mrSges(6,3) + mrSges(7,2)) * t136;
t23 = t102 * t39 - t151 * t115;
t160 = t23 ^ 2;
t159 = 0.2e1 * t23;
t158 = 0.2e1 * t58;
t157 = 0.2e1 * t66;
t74 = -t104 * mrSges(7,1) - t101 * mrSges(7,3);
t156 = 0.2e1 * t74;
t154 = pkin(2) * t99;
t150 = m(7) * t101;
t149 = Ifges(6,6) * t49;
t3 = qJ(6) * t49 + t7;
t147 = t104 * t3;
t86 = pkin(2) * t100 + pkin(3);
t62 = -t102 * t154 + t151 * t86;
t145 = t62 * mrSges(5,1);
t63 = t102 * t86 + t151 * t154;
t144 = t63 * mrSges(5,2);
t130 = t101 * t50;
t28 = -mrSges(6,2) * t49 - mrSges(6,3) * t130;
t31 = -mrSges(7,2) * t130 + mrSges(7,3) * t49;
t142 = t28 + t31;
t129 = t104 * t50;
t29 = mrSges(6,1) * t49 - mrSges(6,3) * t129;
t30 = -t49 * mrSges(7,1) + mrSges(7,2) * t129;
t141 = t29 - t30;
t61 = pkin(9) + t63;
t140 = t136 * pkin(9) * t61;
t139 = t136 * t61 ^ 2;
t137 = t136 * pkin(9) ^ 2;
t135 = t103 ^ 2 + t105 ^ 2;
t89 = t101 * mrSges(7,2);
t128 = qJ(6) * t104;
t126 = -t66 * mrSges(4,1) + t67 * mrSges(4,2);
t123 = Ifges(7,6) * t130 + t166 * t129 + t165 * t49;
t4 = -pkin(5) * t49 - t6;
t122 = t101 * t4 + t147;
t75 = -t104 * mrSges(6,1) + t101 * mrSges(6,2);
t120 = t101 * mrSges(6,1) + t104 * mrSges(6,2);
t119 = t101 * mrSges(7,1) - t104 * mrSges(7,3);
t118 = t104 * pkin(5) + t101 * qJ(6);
t117 = pkin(5) * t101 - t128;
t73 = -pkin(4) - t118;
t78 = -Ifges(7,3) * t104 + t132;
t79 = Ifges(6,2) * t104 + t134;
t116 = Ifges(5,3) + (-t78 + t79) * t104 + t163 * t101;
t113 = -t141 * t101 + t142 * t104;
t112 = mrSges(7,2) * t128 - pkin(5) * t89 + t162;
t111 = 0.2e1 * t161;
t110 = -m(7) * t117 - t119 - t120;
t14 = Ifges(7,6) * t49 + (Ifges(7,3) * t101 + t131) * t50;
t15 = t149 + (-Ifges(6,2) * t101 + t133) * t50;
t9 = t117 * t50 + t23;
t109 = -t25 * mrSges(5,2) + mrSges(7,2) * t147 + Ifges(5,5) * t50 + t4 * t89 + t9 * t74 + (t75 - mrSges(5,1)) * t23 + t164 * t101 / 0.2e1 + (t78 / 0.2e1 - t79 / 0.2e1) * t130 + t163 * t129 / 0.2e1 + (-t14 / 0.2e1 + t15 / 0.2e1) * t104 + t121 * mrSges(6,3) + (-Ifges(5,6) + t162 / 0.2e1) * t49;
t60 = -pkin(4) - t62;
t51 = t73 - t62;
t44 = t50 * mrSges(5,2);
t27 = t120 * t50;
t26 = t119 * t50;
t1 = [-0.2e1 * pkin(1) * (-t105 * mrSges(3,1) + t103 * mrSges(3,2)) + t103 * (Ifges(3,1) * t103 + Ifges(3,4) * t105) + t105 * (Ifges(3,4) * t103 + Ifges(3,2) * t105) + t53 * mrSges(4,3) * t157 + 0.2e1 * t87 * t126 + Ifges(4,2) * t66 ^ 2 + t44 * t158 + 0.2e1 * t9 * t26 + t27 * t159 + 0.2e1 * t7 * t28 + 0.2e1 * t6 * t29 + 0.2e1 * t4 * t30 + 0.2e1 * t3 * t31 + Ifges(2,3) + 0.2e1 * t135 * pkin(7) * mrSges(3,3) + (-0.2e1 * t52 * mrSges(4,3) + Ifges(4,1) * t67 + Ifges(4,4) * t157) * t67 + (mrSges(5,1) * t158 - 0.2e1 * t25 * mrSges(5,3) + Ifges(5,2) * t49 + t123) * t49 + (mrSges(5,3) * t159 + Ifges(5,1) * t50 - 0.2e1 * Ifges(5,4) * t49 + t164 * t104 + (t14 - t15 - t149) * t101) * t50 + m(3) * (t135 * pkin(7) ^ 2 + pkin(1) ^ 2) + m(4) * (t52 ^ 2 + t53 ^ 2 + t87 ^ 2) + m(5) * (t25 ^ 2 + t58 ^ 2 + t160) + m(6) * (t6 ^ 2 + t7 ^ 2 + t160) + m(7) * (t3 ^ 2 + t4 ^ 2 + t9 ^ 2); t113 * t61 + m(5) * (-t23 * t62 + t25 * t63) + (-t103 * mrSges(3,1) - t105 * mrSges(3,2)) * pkin(7) + (-t63 * t49 - t62 * t50) * mrSges(5,3) + Ifges(3,6) * t105 + Ifges(3,5) * t103 + Ifges(4,6) * t66 + Ifges(4,5) * t67 + t51 * t26 + t52 * mrSges(4,1) - t53 * mrSges(4,2) + t60 * t27 + m(7) * (t122 * t61 + t51 * t9) + m(6) * (t121 * t61 + t23 * t60) + t109 + ((-t100 * t67 + t66 * t99) * mrSges(4,3) + m(4) * (t100 * t52 + t53 * t99)) * pkin(2); 0.2e1 * t145 - 0.2e1 * t144 + t51 * t156 + 0.2e1 * t60 * t75 + Ifges(3,3) + Ifges(4,3) + t111 * t61 + m(7) * (t51 ^ 2 + t139) + m(6) * (t60 ^ 2 + t139) + m(5) * (t62 ^ 2 + t63 ^ 2) + t116 + (0.2e1 * mrSges(4,1) * t100 - 0.2e1 * mrSges(4,2) * t99 + m(4) * (t100 ^ 2 + t99 ^ 2) * pkin(2)) * pkin(2); t49 * mrSges(5,1) + t44 + t141 * t104 + t142 * t101 + m(7) * (t101 * t3 - t104 * t4) + m(6) * (t101 * t7 + t104 * t6) + m(5) * t58 + m(4) * t87 + t126; 0; m(4) + m(5) + 0.2e1 * (m(6) / 0.2e1 + m(7) / 0.2e1) * t136; t113 * pkin(9) + t73 * t26 - pkin(4) * t27 + m(7) * (pkin(9) * t122 + t73 * t9) + m(6) * (-pkin(4) * t23 + pkin(9) * t121) + t109; t145 - t144 + (t60 - pkin(4)) * t75 + (t51 + t73) * t74 + m(7) * (t51 * t73 + t140) + m(6) * (-pkin(4) * t60 + t140) + t116 + (pkin(9) + t61) * t161; 0; -0.2e1 * pkin(4) * t75 + t73 * t156 + m(7) * (t73 ^ 2 + t137) + m(6) * (pkin(4) ^ 2 + t137) + t111 * pkin(9) + t116; -Ifges(6,6) * t130 - pkin(5) * t30 + m(7) * (-pkin(5) * t4 + qJ(6) * t3) + qJ(6) * t31 + t3 * mrSges(7,3) - t7 * mrSges(6,2) + t6 * mrSges(6,1) - t4 * mrSges(7,1) + t123; t110 * t61 + t112; m(7) * t118 - t74 - t75; pkin(9) * t110 + t112; 0.2e1 * pkin(5) * mrSges(7,1) + 0.2e1 * qJ(6) * mrSges(7,3) + m(7) * (pkin(5) ^ 2 + qJ(6) ^ 2) + t165; m(7) * t4 + t30; t61 * t150 + t89; -m(7) * t104; pkin(9) * t150 + t89; -m(7) * pkin(5) - mrSges(7,1); m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
