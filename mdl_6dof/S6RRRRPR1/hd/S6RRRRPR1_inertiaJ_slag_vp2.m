% Calculate joint inertia matrix for
% S6RRRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-03-09 21:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPR1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR1_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR1_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR1_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR1_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR1_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:52:59
% EndTime: 2019-03-09 21:53:01
% DurationCPUTime: 1.08s
% Computational Cost: add. (3593->259), mult. (6674->383), div. (0->0), fcn. (7805->10), ass. (0->111)
t103 = cos(qJ(6));
t99 = sin(qJ(6));
t135 = t99 * mrSges(7,3);
t100 = sin(qJ(4));
t104 = cos(qJ(4));
t101 = sin(qJ(3));
t102 = sin(qJ(2));
t105 = cos(qJ(3));
t106 = cos(qJ(2));
t73 = -t101 * t102 + t105 * t106;
t74 = t101 * t106 + t105 * t102;
t47 = -t100 * t74 + t104 * t73;
t48 = t100 * t73 + t104 * t74;
t97 = sin(pkin(11));
t98 = cos(pkin(11));
t30 = -t98 * t47 + t97 * t48;
t31 = t97 * t47 + t98 * t48;
t15 = -t30 * mrSges(7,2) - t31 * t135;
t127 = t103 * t31;
t16 = t30 * mrSges(7,1) - mrSges(7,3) * t127;
t116 = t103 * t15 - t99 * t16;
t77 = -t103 * mrSges(7,1) + t99 * mrSges(7,2);
t147 = t98 * pkin(4);
t85 = -pkin(5) - t147;
t68 = t85 * t77;
t93 = t99 ^ 2;
t146 = mrSges(7,3) * t93;
t84 = t97 * pkin(4) + pkin(10);
t75 = t84 * t146;
t95 = t103 ^ 2;
t145 = mrSges(7,3) * t95;
t76 = t84 * t145;
t89 = mrSges(6,1) * t147;
t153 = t68 + t75 + t76 + t89;
t126 = t100 * t101;
t87 = t105 * pkin(2) + pkin(3);
t66 = -pkin(2) * t126 + t104 * t87;
t63 = pkin(4) + t66;
t125 = t101 * t104;
t67 = pkin(2) * t125 + t100 * t87;
t42 = t98 * t63 - t97 * t67;
t38 = -pkin(5) - t42;
t33 = t38 * t77;
t43 = t97 * t63 + t98 * t67;
t39 = pkin(10) + t43;
t34 = t39 * t146;
t35 = t39 * t145;
t36 = t42 * mrSges(6,1);
t60 = t66 * mrSges(5,1);
t152 = t33 + t34 + t35 + t36 + t60;
t148 = -pkin(8) - pkin(7);
t22 = -t48 * pkin(9) + ((t100 * t105 + t125) * t106 + (t104 * t105 - t126) * t102) * t148;
t111 = -t48 * qJ(5) + t22;
t122 = t148 * t102;
t123 = t148 * t106;
t51 = t101 * t123 + t105 * t122;
t52 = t101 * t122 - t105 * t123;
t23 = t104 * (t73 * pkin(9) + t52) + t100 * (-t74 * pkin(9) + t51);
t18 = qJ(5) * t47 + t23;
t6 = -t98 * t111 + t97 * t18;
t151 = t6 ^ 2;
t150 = 0.2e1 * t6;
t88 = -t106 * pkin(2) - pkin(1);
t57 = -t73 * pkin(3) + t88;
t32 = -t47 * pkin(4) + t57;
t149 = 0.2e1 * t32;
t144 = Ifges(7,4) * t99;
t143 = t100 * pkin(3);
t13 = t30 * pkin(5) - t31 * pkin(10) + t32;
t8 = t97 * t111 + t98 * t18;
t3 = t103 * t8 + t99 * t13;
t142 = t103 * t3;
t141 = t104 * pkin(3);
t140 = t31 * t99;
t139 = t43 * mrSges(6,2);
t86 = pkin(4) + t141;
t65 = t98 * t143 + t97 * t86;
t138 = t65 * mrSges(6,2);
t137 = t67 * mrSges(5,2);
t136 = t97 * mrSges(6,2);
t133 = Ifges(7,5) * t127 + Ifges(7,3) * t30;
t132 = Ifges(7,5) * t99 + Ifges(7,6) * t103;
t131 = t93 + t95;
t130 = t102 ^ 2 + t106 ^ 2;
t129 = Ifges(7,4) * t103;
t124 = mrSges(5,2) * t143;
t121 = t131 * t84;
t78 = Ifges(7,2) * t103 + t144;
t79 = Ifges(7,1) * t99 + t129;
t120 = t103 * t78 + t99 * t79 + Ifges(5,3) + Ifges(6,3);
t119 = Ifges(4,3) + t120;
t2 = t103 * t13 - t99 * t8;
t118 = -t2 * t99 + t142;
t117 = mrSges(7,1) * t99 + mrSges(7,2) * t103;
t64 = -t97 * t143 + t98 * t86;
t115 = (t105 * mrSges(4,1) - t101 * mrSges(4,2)) * pkin(2);
t61 = -pkin(5) - t64;
t49 = t61 * t77;
t62 = pkin(10) + t65;
t53 = t62 * t146;
t54 = t62 * t145;
t58 = t64 * mrSges(6,1);
t90 = mrSges(5,1) * t141;
t114 = t49 + t53 + t54 + t58 + t90 + t120;
t11 = Ifges(7,6) * t30 + (-Ifges(7,2) * t99 + t129) * t31;
t12 = Ifges(7,5) * t30 + (Ifges(7,1) * t103 - t144) * t31;
t113 = -t23 * mrSges(5,2) - t8 * mrSges(6,2) + mrSges(7,3) * t142 - t2 * t135 + t103 * t11 / 0.2e1 + t22 * mrSges(5,1) - t78 * t140 / 0.2e1 + t79 * t127 / 0.2e1 + Ifges(6,5) * t31 + Ifges(5,6) * t47 + Ifges(5,5) * t48 + t99 * t12 / 0.2e1 + (-mrSges(6,1) + t77) * t6 + (t132 / 0.2e1 - Ifges(6,6)) * t30;
t112 = t51 * mrSges(4,1) - t52 * mrSges(4,2) + Ifges(4,5) * t74 + Ifges(4,6) * t73 + t113;
t26 = t31 * mrSges(6,2);
t14 = t117 * t31;
t1 = [t48 * (Ifges(5,1) * t48 + Ifges(5,4) * t47) + t47 * (Ifges(5,4) * t48 + Ifges(5,2) * t47) + 0.2e1 * t57 * (-mrSges(5,1) * t47 + mrSges(5,2) * t48) + t26 * t149 + t14 * t150 + 0.2e1 * t3 * t15 + 0.2e1 * t2 * t16 + Ifges(2,3) + t106 * (Ifges(3,4) * t102 + Ifges(3,2) * t106) - 0.2e1 * pkin(1) * (-t106 * mrSges(3,1) + t102 * mrSges(3,2)) + t102 * (Ifges(3,1) * t102 + Ifges(3,4) * t106) + 0.2e1 * t88 * (-mrSges(4,1) * t73 + mrSges(4,2) * t74) + t73 * (Ifges(4,4) * t74 + Ifges(4,2) * t73) + t74 * (Ifges(4,1) * t74 + Ifges(4,4) * t73) + (mrSges(6,1) * t149 - 0.2e1 * t8 * mrSges(6,3) + Ifges(6,2) * t30 + t133) * t30 + (mrSges(6,3) * t150 + Ifges(6,1) * t31 + t103 * t12 - t99 * t11 + (-Ifges(7,6) * t99 - (2 * Ifges(6,4))) * t30) * t31 + m(3) * (t130 * pkin(7) ^ 2 + pkin(1) ^ 2) + m(4) * (t51 ^ 2 + t52 ^ 2 + t88 ^ 2) + m(5) * (t22 ^ 2 + t23 ^ 2 + t57 ^ 2) + m(6) * (t32 ^ 2 + t8 ^ 2 + t151) + m(7) * (t2 ^ 2 + t3 ^ 2 + t151) + 0.2e1 * (-t22 * t48 + t23 * t47) * mrSges(5,3) + 0.2e1 * (-t51 * t74 + t52 * t73) * mrSges(4,3) + 0.2e1 * t130 * pkin(7) * mrSges(3,3); m(7) * (t118 * t39 + t38 * t6) + t116 * t39 + m(5) * (t22 * t66 + t23 * t67) + m(6) * (-t42 * t6 + t43 * t8) + t112 + (-t102 * mrSges(3,1) - t106 * mrSges(3,2)) * pkin(7) + (-t43 * t30 - t42 * t31) * mrSges(6,3) + (t67 * t47 - t66 * t48) * mrSges(5,3) + t38 * t14 + Ifges(3,6) * t106 + Ifges(3,5) * t102 + (m(4) * (t101 * t52 + t105 * t51) + (t101 * t73 - t105 * t74) * mrSges(4,3)) * pkin(2); -0.2e1 * t137 - 0.2e1 * t139 + Ifges(3,3) + 0.2e1 * t33 + 0.2e1 * t34 + 0.2e1 * t35 + 0.2e1 * t36 + 0.2e1 * t60 + 0.2e1 * t115 + m(7) * (t131 * t39 ^ 2 + t38 ^ 2) + m(6) * (t42 ^ 2 + t43 ^ 2) + m(5) * (t66 ^ 2 + t67 ^ 2) + m(4) * (t101 ^ 2 + t105 ^ 2) * pkin(2) ^ 2 + t119; m(7) * (t118 * t62 + t61 * t6) + (m(5) * (t100 * t23 + t104 * t22) + (t100 * t47 - t104 * t48) * mrSges(5,3)) * pkin(3) + t112 + t61 * t14 + t116 * t62 + m(6) * (-t6 * t64 + t65 * t8) + (-t65 * t30 - t64 * t31) * mrSges(6,3); m(7) * (t131 * t62 * t39 + t61 * t38) + t114 + m(5) * (t100 * t67 + t104 * t66) * pkin(3) + m(6) * (t42 * t64 + t43 * t65) + t115 + (-t43 - t65) * mrSges(6,2) + (-t67 - t143) * mrSges(5,2) + Ifges(4,3) + t152; -0.2e1 * t124 - 0.2e1 * t138 + 0.2e1 * t49 + 0.2e1 * t53 + 0.2e1 * t54 + 0.2e1 * t58 + 0.2e1 * t90 + m(7) * (t131 * t62 ^ 2 + t61 ^ 2) + m(6) * (t64 ^ 2 + t65 ^ 2) + m(5) * (t100 ^ 2 + t104 ^ 2) * pkin(3) ^ 2 + t119; (m(6) * (-t6 * t98 + t8 * t97) + (-t97 * t30 - t98 * t31) * mrSges(6,3)) * pkin(4) + t113 + (m(7) * t6 + t14) * t85 + (m(7) * t118 + t116) * t84; m(7) * (t39 * t121 + t85 * t38) - t137 - t139 + (m(6) * (t42 * t98 + t43 * t97) - t136) * pkin(4) + t120 + t152 + t153; -t124 + m(7) * (t62 * t121 + t85 * t61) + (-t136 + m(6) * (t64 * t98 + t65 * t97)) * pkin(4) + t114 - t138 + t153; 0.2e1 * t68 + 0.2e1 * t75 + 0.2e1 * t76 + 0.2e1 * t89 + m(7) * (t131 * t84 ^ 2 + t85 ^ 2) + t120 + (-0.2e1 * t136 + m(6) * (t97 ^ 2 + t98 ^ 2) * pkin(4)) * pkin(4); t30 * mrSges(6,1) + t103 * t16 + t99 * t15 + t26 + m(7) * (t103 * t2 + t99 * t3) + m(6) * t32; 0; 0; 0; m(7) * t131 + m(6); t2 * mrSges(7,1) - t3 * mrSges(7,2) - Ifges(7,6) * t140 + t133; -t117 * t39 + t132; -t117 * t62 + t132; -t117 * t84 + t132; -t77; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
