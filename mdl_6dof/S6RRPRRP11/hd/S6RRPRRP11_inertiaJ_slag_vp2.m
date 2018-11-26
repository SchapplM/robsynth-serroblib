% Calculate joint inertia matrix for
% S6RRPRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 17:18
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRP11_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP11_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP11_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP11_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP11_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP11_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:18:27
% EndTime: 2018-11-23 17:18:28
% DurationCPUTime: 1.08s
% Computational Cost: add. (1338->277), mult. (2389->364), div. (0->0), fcn. (2231->6), ass. (0->103)
t147 = pkin(3) + pkin(7);
t134 = Ifges(6,3) + Ifges(7,3);
t101 = sin(qJ(2));
t104 = cos(qJ(2));
t146 = t101 ^ 2 + t104 ^ 2;
t100 = sin(qJ(4));
t102 = cos(qJ(5));
t103 = cos(qJ(4));
t99 = sin(qJ(5));
t113 = t100 * t99 - t102 * t103;
t64 = -t102 * t100 - t99 * t103;
t145 = t113 ^ 2 + t64 ^ 2;
t144 = -m(4) * pkin(2) + mrSges(4,2);
t136 = mrSges(6,2) + mrSges(7,2);
t143 = (t102 * mrSges(6,1) - t136 * t99) * pkin(4);
t141 = 2 * mrSges(7,1);
t140 = m(7) * pkin(5);
t105 = -pkin(2) - pkin(8);
t139 = m(7) * t99;
t138 = -t100 / 0.2e1;
t119 = -qJ(3) * t101 - pkin(1);
t52 = t105 * t104 + t119;
t77 = t147 * t101;
t68 = t103 * t77;
t16 = pkin(4) * t101 + t68 + (pkin(9) * t104 - t52) * t100;
t124 = t103 * t104;
t23 = t100 * t77 + t103 * t52;
t20 = -pkin(9) * t124 + t23;
t6 = t102 * t20 + t99 * t16;
t137 = t64 * t99;
t133 = -pkin(9) + t105;
t46 = t113 * t104;
t34 = -mrSges(7,2) * t101 + mrSges(7,3) * t46;
t35 = -mrSges(6,2) * t101 + mrSges(6,3) * t46;
t132 = t34 + t35;
t71 = t133 * t100;
t72 = t133 * t103;
t32 = t102 * t71 + t99 * t72;
t131 = t146 * pkin(7) ^ 2;
t78 = t147 * t104;
t130 = t100 ^ 2 + t103 ^ 2;
t129 = Ifges(5,4) * t100;
t128 = Ifges(5,4) * t103;
t126 = t102 * t113;
t125 = t103 * mrSges(5,1);
t80 = t100 * pkin(4) + qJ(3);
t50 = pkin(4) * t124 + t78;
t47 = t64 * t104;
t5 = t102 * t16 - t20 * t99;
t2 = pkin(5) * t101 - qJ(6) * t47 + t5;
t36 = mrSges(7,1) * t101 - mrSges(7,3) * t47;
t123 = m(7) * t2 + t36;
t122 = m(5) * t130;
t120 = t130 * mrSges(5,3);
t18 = -t46 * mrSges(7,1) + t47 * mrSges(7,2);
t25 = -t64 * mrSges(7,1) - mrSges(7,2) * t113;
t31 = t102 * t72 - t71 * t99;
t14 = qJ(6) * t113 + t31;
t117 = m(7) * t14 + mrSges(7,3) * t113;
t116 = -t100 * mrSges(5,2) + t125;
t115 = -Ifges(5,5) * t100 - Ifges(5,6) * t103;
t22 = -t100 * t52 + t68;
t114 = t23 * t100 + t22 * t103;
t112 = (Ifges(6,5) + Ifges(7,5)) * t47 + (Ifges(6,6) + Ifges(7,6)) * t46 + t134 * t101;
t111 = t136 * t64 + (-mrSges(6,1) - mrSges(7,1)) * t113;
t15 = qJ(6) * t64 + t32;
t56 = Ifges(7,6) * t64;
t57 = Ifges(6,6) * t64;
t58 = Ifges(7,5) * t113;
t59 = Ifges(6,5) * t113;
t110 = t31 * mrSges(6,1) + t14 * mrSges(7,1) - t32 * mrSges(6,2) - t15 * mrSges(7,2) + t56 + t57 - t58 - t59;
t3 = qJ(6) * t46 + t6;
t109 = t5 * mrSges(6,1) + t2 * mrSges(7,1) - t6 * mrSges(6,2) - t3 * mrSges(7,2) + t112;
t108 = pkin(4) ^ 2;
t106 = qJ(3) ^ 2;
t88 = t99 ^ 2 * t108;
t87 = Ifges(5,5) * t103;
t86 = Ifges(5,3) * t101;
t83 = pkin(4) * t102 + pkin(5);
t76 = Ifges(5,1) * t103 - t129;
t75 = -Ifges(5,2) * t100 + t128;
t74 = mrSges(5,1) * t100 + mrSges(5,2) * t103;
t73 = -pkin(2) * t104 + t119;
t70 = -mrSges(5,2) * t101 - mrSges(5,3) * t124;
t69 = mrSges(5,3) * t100 * t104 + mrSges(5,1) * t101;
t51 = t116 * t104;
t49 = pkin(4) * t137;
t45 = Ifges(5,5) * t101 + (-Ifges(5,1) * t100 - t128) * t104;
t44 = Ifges(5,6) * t101 + (-Ifges(5,2) * t103 - t129) * t104;
t38 = -pkin(5) * t64 + t80;
t37 = mrSges(6,1) * t101 - mrSges(6,3) * t47;
t30 = -Ifges(6,1) * t113 + Ifges(6,4) * t64;
t29 = -Ifges(7,1) * t113 + Ifges(7,4) * t64;
t28 = -Ifges(6,4) * t113 + Ifges(6,2) * t64;
t27 = -Ifges(7,4) * t113 + Ifges(7,2) * t64;
t26 = -mrSges(6,1) * t64 - mrSges(6,2) * t113;
t21 = -pkin(5) * t46 + t50;
t19 = -mrSges(6,1) * t46 + mrSges(6,2) * t47;
t13 = Ifges(6,1) * t47 + Ifges(6,4) * t46 + Ifges(6,5) * t101;
t12 = Ifges(7,1) * t47 + Ifges(7,4) * t46 + Ifges(7,5) * t101;
t11 = Ifges(6,4) * t47 + Ifges(6,2) * t46 + Ifges(6,6) * t101;
t10 = Ifges(7,4) * t47 + Ifges(7,2) * t46 + Ifges(7,6) * t101;
t1 = [0.2e1 * t21 * t18 + 0.2e1 * t50 * t19 + 0.2e1 * t2 * t36 + 0.2e1 * t22 * t69 + 0.2e1 * t23 * t70 + 0.2e1 * t3 * t34 + 0.2e1 * t6 * t35 + 0.2e1 * t5 * t37 + 0.2e1 * t78 * t51 + Ifges(2,3) + (t12 + t13) * t47 + (t10 + t11) * t46 + (0.2e1 * pkin(1) * mrSges(3,1) + 0.2e1 * t73 * mrSges(4,2) - t100 * t45 - t103 * t44 + (Ifges(4,3) + Ifges(3,2)) * t104) * t104 + (-0.2e1 * pkin(1) * mrSges(3,2) - 0.2e1 * t73 * mrSges(4,3) + t86 + (Ifges(4,2) + Ifges(3,1)) * t101 + ((2 * Ifges(3,4)) + (2 * Ifges(4,6)) + t115) * t104 + t112) * t101 + m(4) * (t73 ^ 2 + t131) + m(3) * (pkin(1) ^ 2 + t131) + m(5) * (t22 ^ 2 + t23 ^ 2 + t78 ^ 2) + m(6) * (t5 ^ 2 + t50 ^ 2 + t6 ^ 2) + m(7) * (t2 ^ 2 + t21 ^ 2 + t3 ^ 2) + 0.2e1 * (mrSges(4,1) + mrSges(3,3)) * pkin(7) * t146; qJ(3) * t51 + t14 * t36 + t15 * t34 + t38 * t18 + t80 * t19 + t21 * t25 + t50 * t26 + t31 * t37 + t32 * t35 + t78 * t74 + (t29 / 0.2e1 + t30 / 0.2e1) * t47 + (t27 / 0.2e1 + t28 / 0.2e1) * t46 + (t105 * t69 - t22 * mrSges(5,3) + t45 / 0.2e1) * t103 + (t105 * t70 - t23 * mrSges(5,3) - t44 / 0.2e1) * t100 + (-pkin(2) * mrSges(4,1) + Ifges(5,6) * t138 + t87 / 0.2e1 - t58 / 0.2e1 + t56 / 0.2e1 - t59 / 0.2e1 + t57 / 0.2e1 - Ifges(4,4) + Ifges(3,5)) * t101 + m(5) * (qJ(3) * t78 + t105 * t114) + m(6) * (t31 * t5 + t32 * t6 + t50 * t80) + m(7) * (t14 * t2 + t15 * t3 + t38 * t21) - (t12 / 0.2e1 + t13 / 0.2e1 - t2 * mrSges(7,3) - t5 * mrSges(6,3)) * t113 + (t10 / 0.2e1 + t11 / 0.2e1 + t3 * mrSges(7,3) + t6 * mrSges(6,3)) * t64 + (-t103 * t75 / 0.2e1 + t76 * t138 + qJ(3) * mrSges(4,1) - Ifges(4,5) + Ifges(3,6)) * t104 + ((m(4) * qJ(3) - mrSges(3,2) + mrSges(4,3)) * t104 + (-mrSges(3,1) + t144) * t101) * pkin(7); -0.2e1 * pkin(2) * mrSges(4,2) - t100 * t75 + t103 * t76 + 0.2e1 * t38 * t25 + 0.2e1 * t80 * t26 + Ifges(4,1) + Ifges(3,3) - (-0.2e1 * mrSges(6,3) * t31 - 0.2e1 * mrSges(7,3) * t14 + t29 + t30) * t113 + (0.2e1 * mrSges(6,3) * t32 + 0.2e1 * mrSges(7,3) * t15 + t27 + t28) * t64 + m(7) * (t14 ^ 2 + t15 ^ 2 + t38 ^ 2) + m(6) * (t31 ^ 2 + t32 ^ 2 + t80 ^ 2) + m(5) * (t105 ^ 2 * t130 + t106) + m(4) * (pkin(2) ^ 2 + t106) + 0.2e1 * (t74 + mrSges(4,3)) * qJ(3) - 0.2e1 * t105 * t120; t100 * t70 + t103 * t69 - (t36 + t37) * t113 - t132 * t64 + (m(4) * pkin(7) + mrSges(4,1)) * t101 + m(7) * (-t113 * t2 - t3 * t64) + m(6) * (-t113 * t5 - t6 * t64) + m(5) * t114; -t120 + m(7) * (-t113 * t14 - t15 * t64) + m(6) * (-t113 * t31 - t32 * t64) + t105 * t122 + t145 * (-mrSges(6,3) - mrSges(7,3)) + t144; m(4) + t122 + 0.2e1 * (m(6) / 0.2e1 + m(7) / 0.2e1) * t145; t109 + t115 * t104 + t123 * t83 + (t102 * t37 + t132 * t99 + t3 * t139 + m(6) * (t102 * t5 + t6 * t99)) * pkin(4) + t86 + t22 * mrSges(5,1) - t23 * mrSges(5,2); t105 * t125 + t87 + t117 * t83 + (-mrSges(5,2) * t105 - Ifges(5,6)) * t100 + (mrSges(7,3) * t137 + (t126 + t137) * mrSges(6,3) + t15 * t139 + m(6) * (t102 * t31 + t32 * t99)) * pkin(4) + t110; m(6) * (-pkin(4) * t126 - t49) + m(7) * (-t113 * t83 - t49) + t111 + t116; t83 * t141 + Ifges(5,3) + m(7) * (t83 ^ 2 + t88) + m(6) * (t102 ^ 2 * t108 + t88) + 0.2e1 * t143 + t134; pkin(5) * t123 + t109; pkin(5) * t117 + t110; -t113 * t140 + t111; t83 * t140 + (pkin(5) + t83) * mrSges(7,1) + t143 + t134; (t141 + t140) * pkin(5) + t134; m(7) * t21 + t18; m(7) * t38 + t25; 0; 0; 0; m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
