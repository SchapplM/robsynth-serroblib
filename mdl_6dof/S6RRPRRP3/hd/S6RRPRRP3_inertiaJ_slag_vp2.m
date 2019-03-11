% Calculate joint inertia matrix for
% S6RRPRRP3
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
% Datum: 2019-03-09 11:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRP3_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP3_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP3_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP3_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP3_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP3_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:47:22
% EndTime: 2019-03-09 11:47:26
% DurationCPUTime: 1.24s
% Computational Cost: add. (2073->273), mult. (3928->385), div. (0->0), fcn. (4307->8), ass. (0->102)
t149 = Ifges(6,5) + Ifges(7,5);
t148 = Ifges(6,6) + Ifges(7,6);
t134 = Ifges(6,3) + Ifges(7,3);
t107 = cos(qJ(4));
t101 = sin(pkin(10));
t102 = cos(pkin(10));
t105 = sin(qJ(2));
t108 = cos(qJ(2));
t79 = t101 * t108 + t102 * t105;
t122 = t107 * t79;
t78 = t101 * t105 - t102 * t108;
t147 = Ifges(5,5) * t122 + Ifges(5,3) * t78;
t103 = sin(qJ(5));
t104 = sin(qJ(4));
t106 = cos(qJ(5));
t81 = -t103 * t104 + t106 * t107;
t82 = t103 * t107 + t104 * t106;
t146 = t148 * t81 + t149 * t82;
t145 = (t106 * mrSges(6,1) + (-mrSges(6,2) - mrSges(7,2)) * t103) * pkin(4);
t133 = -qJ(3) - pkin(7);
t118 = t133 * t105;
t85 = t133 * t108;
t55 = -t101 * t85 - t102 * t118;
t144 = t55 ^ 2;
t143 = 2 * mrSges(7,1);
t142 = 0.2e1 * t55;
t93 = -pkin(2) * t108 - pkin(1);
t141 = 0.2e1 * t93;
t140 = m(7) * pkin(5);
t138 = pkin(5) * t81;
t90 = pkin(2) * t101 + pkin(8);
t137 = pkin(9) + t90;
t124 = t104 * t79;
t44 = pkin(3) * t78 - pkin(8) * t79 + t93;
t57 = t101 * t118 - t102 * t85;
t20 = t104 * t44 + t107 * t57;
t15 = -pkin(9) * t124 + t20;
t19 = -t104 * t57 + t107 * t44;
t9 = pkin(4) * t78 - pkin(9) * t122 + t19;
t6 = t103 * t9 + t106 * t15;
t136 = m(7) * t103;
t135 = pkin(4) * t106;
t36 = t82 * t79;
t21 = -mrSges(7,2) * t78 - mrSges(7,3) * t36;
t22 = -mrSges(6,2) * t78 - mrSges(6,3) * t36;
t132 = t21 + t22;
t69 = t137 * t104;
t70 = t137 * t107;
t43 = -t103 * t69 + t106 * t70;
t48 = -t81 * mrSges(7,1) + t82 * mrSges(7,2);
t129 = Ifges(5,5) * t104 + Ifges(5,6) * t107;
t128 = t104 ^ 2 + t107 ^ 2;
t127 = Ifges(5,4) * t104;
t126 = Ifges(5,4) * t107;
t125 = t103 * t81;
t121 = t105 ^ 2 + t108 ^ 2;
t37 = t81 * t79;
t5 = -t103 * t15 + t106 * t9;
t2 = pkin(5) * t78 - qJ(6) * t37 + t5;
t23 = mrSges(7,1) * t78 - mrSges(7,3) * t37;
t120 = m(7) * t2 + t23;
t91 = -pkin(2) * t102 - pkin(3);
t17 = t36 * mrSges(7,1) + t37 * mrSges(7,2);
t49 = -t81 * mrSges(6,1) + t82 * mrSges(6,2);
t42 = -t103 * t70 - t106 * t69;
t30 = pkin(4) * t124 + t55;
t28 = -qJ(6) * t82 + t42;
t117 = m(7) * t28 - t82 * mrSges(7,3);
t84 = -t107 * mrSges(5,1) + t104 * mrSges(5,2);
t116 = mrSges(5,1) * t104 + mrSges(5,2) * t107;
t115 = -t49 - t48;
t114 = t134 * t78 - t148 * t36 + t149 * t37;
t83 = -pkin(4) * t107 + t91;
t29 = qJ(6) * t81 + t43;
t113 = t42 * mrSges(6,1) + t28 * mrSges(7,1) - t43 * mrSges(6,2) - t29 * mrSges(7,2) + t146;
t3 = -qJ(6) * t36 + t6;
t112 = t5 * mrSges(6,1) + t2 * mrSges(7,1) - t6 * mrSges(6,2) - t3 * mrSges(7,2) + t114;
t110 = pkin(4) ^ 2;
t96 = t103 ^ 2 * t110;
t92 = pkin(5) + t135;
t87 = Ifges(5,1) * t104 + t126;
t86 = Ifges(5,2) * t107 + t127;
t65 = t79 * mrSges(4,2);
t64 = t103 * pkin(4) * t82;
t58 = t83 - t138;
t53 = Ifges(6,1) * t82 + Ifges(6,4) * t81;
t52 = Ifges(7,1) * t82 + Ifges(7,4) * t81;
t51 = Ifges(6,4) * t82 + Ifges(6,2) * t81;
t50 = Ifges(7,4) * t82 + Ifges(7,2) * t81;
t46 = mrSges(5,1) * t78 - mrSges(5,3) * t122;
t45 = -mrSges(5,2) * t78 - mrSges(5,3) * t124;
t41 = t116 * t79;
t27 = Ifges(5,5) * t78 + (Ifges(5,1) * t107 - t127) * t79;
t26 = Ifges(5,6) * t78 + (-Ifges(5,2) * t104 + t126) * t79;
t24 = mrSges(6,1) * t78 - mrSges(6,3) * t37;
t18 = mrSges(6,1) * t36 + mrSges(6,2) * t37;
t16 = pkin(5) * t36 + t30;
t13 = Ifges(6,1) * t37 - Ifges(6,4) * t36 + Ifges(6,5) * t78;
t12 = Ifges(7,1) * t37 - Ifges(7,4) * t36 + Ifges(7,5) * t78;
t11 = Ifges(6,4) * t37 - Ifges(6,2) * t36 + Ifges(6,6) * t78;
t10 = Ifges(7,4) * t37 - Ifges(7,2) * t36 + Ifges(7,6) * t78;
t1 = [-0.2e1 * pkin(1) * (-t108 * mrSges(3,1) + t105 * mrSges(3,2)) + t105 * (Ifges(3,1) * t105 + Ifges(3,4) * t108) + t108 * (Ifges(3,4) * t105 + Ifges(3,2) * t108) + t65 * t141 + 0.2e1 * t20 * t45 + 0.2e1 * t19 * t46 + t41 * t142 + 0.2e1 * t3 * t21 + 0.2e1 * t6 * t22 + 0.2e1 * t2 * t23 + 0.2e1 * t5 * t24 + 0.2e1 * t30 * t18 + 0.2e1 * t16 * t17 + Ifges(2,3) + (t12 + t13) * t37 - (t10 + t11) * t36 + 0.2e1 * t121 * pkin(7) * mrSges(3,3) + (mrSges(4,3) * t142 + Ifges(4,1) * t79 - t104 * t26 + t107 * t27) * t79 + (mrSges(4,1) * t141 - 0.2e1 * t57 * mrSges(4,3) + Ifges(4,2) * t78 + (-Ifges(5,6) * t104 - (2 * Ifges(4,4))) * t79 + t114 + t147) * t78 + m(3) * (t121 * pkin(7) ^ 2 + pkin(1) ^ 2) + m(4) * (t57 ^ 2 + t93 ^ 2 + t144) + m(5) * (t19 ^ 2 + t20 ^ 2 + t144) + m(6) * (t30 ^ 2 + t5 ^ 2 + t6 ^ 2) + m(7) * (t16 ^ 2 + t2 ^ 2 + t3 ^ 2); (t84 - mrSges(4,1)) * t55 + (t26 / 0.2e1 + t79 * t87 / 0.2e1 + t90 * t45 + t20 * mrSges(5,3)) * t107 + (t27 / 0.2e1 - t79 * t86 / 0.2e1 - t90 * t46 - t19 * mrSges(5,3)) * t104 + ((-t101 * t78 - t102 * t79) * mrSges(4,3) + m(4) * (t101 * t57 - t102 * t55)) * pkin(2) + m(6) * (t30 * t83 + t42 * t5 + t43 * t6) + m(7) * (t16 * t58 + t2 * t28 + t29 * t3) + m(5) * (t55 * t91 + (-t19 * t104 + t20 * t107) * t90) + (t13 / 0.2e1 + t12 / 0.2e1 - t2 * mrSges(7,3) - t5 * mrSges(6,3)) * t82 + (t10 / 0.2e1 + t11 / 0.2e1 + t3 * mrSges(7,3) + t6 * mrSges(6,3)) * t81 + (t52 / 0.2e1 + t53 / 0.2e1) * t37 - (t50 / 0.2e1 + t51 / 0.2e1) * t36 + (t129 + t146) * t78 / 0.2e1 + Ifges(3,6) * t108 + Ifges(3,5) * t105 + t91 * t41 + t83 * t18 - Ifges(4,6) * t78 + Ifges(4,5) * t79 + t43 * t22 + t16 * t48 + t30 * t49 - t57 * mrSges(4,2) + t58 * t17 + t42 * t24 + t28 * t23 + t29 * t21 + (-t105 * mrSges(3,1) - t108 * mrSges(3,2)) * pkin(7); t104 * t87 + t107 * t86 + 0.2e1 * t58 * t48 + 0.2e1 * t83 * t49 + 0.2e1 * t91 * t84 + Ifges(3,3) + Ifges(4,3) + (-0.2e1 * mrSges(6,3) * t42 - 0.2e1 * mrSges(7,3) * t28 + t52 + t53) * t82 + (0.2e1 * mrSges(6,3) * t43 + 0.2e1 * mrSges(7,3) * t29 + t50 + t51) * t81 + m(6) * (t42 ^ 2 + t43 ^ 2 + t83 ^ 2) + m(7) * (t28 ^ 2 + t29 ^ 2 + t58 ^ 2) + m(5) * (t128 * t90 ^ 2 + t91 ^ 2) + m(4) * (t101 ^ 2 + t102 ^ 2) * pkin(2) ^ 2 + 0.2e1 * (mrSges(4,1) * t102 - mrSges(4,2) * t101) * pkin(2) + 0.2e1 * t128 * t90 * mrSges(5,3); t78 * mrSges(4,1) + t104 * t45 + t107 * t46 + t65 + t132 * t82 + (t23 + t24) * t81 + m(7) * (t2 * t81 + t3 * t82) + m(6) * (t5 * t81 + t6 * t82) + m(5) * (t104 * t20 + t107 * t19) + m(4) * t93; m(6) * (t42 * t81 + t43 * t82) + m(7) * (t28 * t81 + t29 * t82); m(4) + m(5) * t128 + 0.2e1 * (m(6) / 0.2e1 + m(7) / 0.2e1) * (t81 ^ 2 + t82 ^ 2); t120 * t92 + (t106 * t24 + t132 * t103 + t3 * t136 + m(6) * (t103 * t6 + t106 * t5)) * pkin(4) + t19 * mrSges(5,1) - t20 * mrSges(5,2) + t112 - Ifges(5,6) * t124 + t147; t117 * t92 - t116 * t90 + (mrSges(7,3) * t125 + (-t106 * t82 + t125) * mrSges(6,3) + m(6) * (t103 * t43 + t106 * t42) + t29 * t136) * pkin(4) + t113 + t129; m(6) * (t81 * t135 + t64) + m(7) * (t81 * t92 + t64) + t115 - t84; t92 * t143 + Ifges(5,3) + m(6) * (t106 ^ 2 * t110 + t96) + m(7) * (t92 ^ 2 + t96) + 0.2e1 * t145 + t134; t120 * pkin(5) + t112; t117 * pkin(5) + t113; m(7) * t138 + t115; t92 * t140 + (pkin(5) + t92) * mrSges(7,1) + t145 + t134; (t143 + t140) * pkin(5) + t134; m(7) * t16 + t17; m(7) * t58 + t48; 0; 0; 0; m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
