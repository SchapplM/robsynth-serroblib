% Calculate joint inertia matrix for
% S6RRRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
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
% Datum: 2019-03-09 18:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRR3_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR3_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR3_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR3_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR3_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR3_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:10:57
% EndTime: 2019-03-09 18:11:00
% DurationCPUTime: 1.22s
% Computational Cost: add. (1950->257), mult. (3400->342), div. (0->0), fcn. (3538->8), ass. (0->91)
t86 = cos(qJ(6));
t151 = t86 / 0.2e1;
t82 = sin(qJ(6));
t61 = -t86 * mrSges(7,1) + mrSges(7,2) * t82;
t147 = -mrSges(6,1) + t61;
t138 = -pkin(8) - pkin(7);
t85 = sin(qJ(2));
t113 = t138 * t85;
t89 = cos(qJ(2));
t114 = t138 * t89;
t84 = sin(qJ(3));
t88 = cos(qJ(3));
t33 = t84 * t113 - t88 * t114;
t50 = t84 * t85 - t88 * t89;
t22 = pkin(9) * t50 + t33;
t83 = sin(qJ(5));
t87 = cos(qJ(5));
t31 = -t88 * t113 - t84 * t114;
t51 = t84 * t89 + t85 * t88;
t96 = -t51 * pkin(9) + t31;
t11 = t22 * t83 - t87 * t96;
t123 = t86 * mrSges(7,3);
t13 = t87 * t22 + t83 * t96;
t70 = -t89 * pkin(2) - pkin(1);
t24 = t50 * pkin(3) - t51 * qJ(4) + t70;
t17 = -pkin(4) * t50 - t24;
t26 = -t87 * t50 + t51 * t83;
t27 = t50 * t83 + t51 * t87;
t4 = pkin(5) * t26 - pkin(10) * t27 + t17;
t2 = -t13 * t82 + t4 * t86;
t137 = t2 * t82;
t3 = t13 * t86 + t4 * t82;
t62 = Ifges(7,5) * t82 + Ifges(7,6) * t86;
t134 = Ifges(7,4) * t82;
t63 = Ifges(7,2) * t86 + t134;
t133 = Ifges(7,4) * t86;
t64 = Ifges(7,1) * t82 + t133;
t7 = Ifges(7,6) * t26 + (-Ifges(7,2) * t82 + t133) * t27;
t8 = Ifges(7,5) * t26 + (Ifges(7,1) * t86 - t134) * t27;
t150 = (t64 * t151 - t82 * t63 / 0.2e1 + Ifges(6,5)) * t27 + (t62 / 0.2e1 - Ifges(6,6)) * t26 - t13 * mrSges(6,2) - mrSges(7,3) * t137 + t3 * t123 + t7 * t151 + t82 * t8 / 0.2e1 + t147 * t11;
t117 = t82 ^ 2 + t86 ^ 2;
t109 = mrSges(7,3) * t117;
t131 = t27 * t82;
t15 = -mrSges(7,2) * t26 - mrSges(7,3) * t131;
t16 = mrSges(7,1) * t26 - t27 * t123;
t102 = t86 * t15 - t82 * t16;
t105 = -t86 * t63 - t82 * t64 - Ifges(6,3);
t145 = t11 ^ 2;
t144 = 2 * mrSges(5,1);
t143 = 2 * mrSges(5,3);
t142 = 0.2e1 * t11;
t141 = 0.2e1 * t17;
t140 = -0.2e1 * t61;
t139 = 0.2e1 * t70;
t132 = t11 * t87;
t69 = -pkin(2) * t88 - pkin(3);
t66 = -pkin(4) + t69;
t67 = pkin(2) * t84 + qJ(4);
t37 = t66 * t87 - t67 * t83;
t130 = t37 * mrSges(6,1);
t38 = t83 * t66 + t87 * t67;
t129 = t38 * mrSges(6,2);
t90 = -pkin(3) - pkin(4);
t58 = -qJ(4) * t83 + t87 * t90;
t128 = t58 * mrSges(6,1);
t59 = t87 * qJ(4) + t83 * t90;
t127 = t59 * mrSges(6,2);
t119 = mrSges(5,2) + mrSges(4,3);
t118 = Ifges(7,5) * t86 * t27 + Ifges(7,3) * t26;
t116 = t85 ^ 2 + t89 ^ 2;
t115 = t31 ^ 2 + t33 ^ 2;
t110 = t117 * pkin(10);
t36 = -pkin(10) + t38;
t108 = t117 * t36;
t57 = -pkin(10) + t59;
t107 = t117 * t57;
t106 = t117 * t83;
t104 = t3 * t86 - t137;
t103 = mrSges(7,1) * t82 + mrSges(7,2) * t86;
t101 = -0.2e1 * t109;
t99 = (mrSges(4,1) * t88 - mrSges(4,2) * t84) * pkin(2);
t97 = Ifges(5,2) + Ifges(4,3) - t105;
t73 = t83 * mrSges(6,2);
t95 = -t83 * t109 + t147 * t87 - mrSges(5,1) + t73;
t94 = (Ifges(4,5) + Ifges(5,4)) * t51 + (Ifges(5,6) - Ifges(4,6)) * t50 + (-mrSges(4,2) + mrSges(5,3)) * t33 + (-mrSges(5,1) - mrSges(4,1)) * t31 - t150;
t80 = t87 ^ 2;
t77 = t83 ^ 2;
t56 = pkin(5) - t58;
t35 = pkin(5) - t37;
t14 = t103 * t27;
t1 = [t85 * (Ifges(3,1) * t85 + Ifges(3,4) * t89) + t89 * (Ifges(3,4) * t85 + Ifges(3,2) * t89) - 0.2e1 * pkin(1) * (-mrSges(3,1) * t89 + mrSges(3,2) * t85) + t14 * t142 + 0.2e1 * t3 * t15 + 0.2e1 * t2 * t16 + Ifges(2,3) + 0.2e1 * t116 * pkin(7) * mrSges(3,3) + (mrSges(6,1) * t141 - 0.2e1 * t13 * mrSges(6,3) + Ifges(6,2) * t26 + t118) * t26 + (mrSges(6,2) * t141 + mrSges(6,3) * t142 + Ifges(6,1) * t27 - t82 * t7 + t86 * t8 + (-Ifges(7,6) * t82 - (2 * Ifges(6,4))) * t26) * t27 + (mrSges(4,2) * t139 - 0.2e1 * t24 * mrSges(5,3) + (Ifges(5,1) + Ifges(4,1)) * t51 + 0.2e1 * t119 * t31) * t51 + (mrSges(4,1) * t139 + t24 * t144 + (Ifges(5,3) + Ifges(4,2)) * t50 + 0.2e1 * (-Ifges(4,4) + Ifges(5,5)) * t51 - 0.2e1 * t119 * t33) * t50 + m(3) * (t116 * pkin(7) ^ 2 + pkin(1) ^ 2) + m(4) * (t70 ^ 2 + t115) + m(5) * (t24 ^ 2 + t115) + m(6) * (t13 ^ 2 + t17 ^ 2 + t145) + m(7) * (t2 ^ 2 + t3 ^ 2 + t145); m(7) * (t104 * t36 + t11 * t35) + Ifges(3,6) * t89 + Ifges(3,5) * t85 + t35 * t14 + t102 * t36 + m(5) * (t31 * t69 + t33 * t67) + m(6) * (-t11 * t37 + t13 * t38) + t94 + (-t85 * mrSges(3,1) - t89 * mrSges(3,2)) * pkin(7) + (-t38 * t26 - t37 * t27) * mrSges(6,3) + (-t67 * t50 + t69 * t51) * mrSges(5,2) + (m(4) * (-t31 * t88 + t33 * t84) + (-t50 * t84 - t51 * t88) * mrSges(4,3)) * pkin(2); -0.2e1 * t69 * mrSges(5,1) - 0.2e1 * t130 + 0.2e1 * t129 + t67 * t143 + t35 * t140 + Ifges(3,3) + 0.2e1 * t99 + t36 * t101 + m(7) * (t117 * t36 ^ 2 + t35 ^ 2) + m(6) * (t37 ^ 2 + t38 ^ 2) + m(5) * (t67 ^ 2 + t69 ^ 2) + m(4) * (t84 ^ 2 + t88 ^ 2) * pkin(2) ^ 2 + t97; (-t59 * t26 - t58 * t27) * mrSges(6,3) + (-pkin(3) * t51 - qJ(4) * t50) * mrSges(5,2) + m(7) * (t104 * t57 + t11 * t56) + t56 * t14 + t94 + t102 * t57 + m(6) * (-t11 * t58 + t13 * t59) + m(5) * (-pkin(3) * t31 + qJ(4) * t33); (-t35 - t56) * t61 + t99 + (t67 + qJ(4)) * mrSges(5,3) + (t38 + t59) * mrSges(6,2) + (-t37 - t58) * mrSges(6,1) + (-t69 + pkin(3)) * mrSges(5,1) + m(7) * (t36 * t107 + t35 * t56) + m(6) * (t37 * t58 + t38 * t59) + m(5) * (-pkin(3) * t69 + qJ(4) * t67) + t97 + (-t36 - t57) * t109; pkin(3) * t144 - 0.2e1 * t128 + 0.2e1 * t127 + qJ(4) * t143 + t56 * t140 + t57 * t101 + m(7) * (t117 * t57 ^ 2 + t56 ^ 2) + m(6) * (t58 ^ 2 + t59 ^ 2) + m(5) * (pkin(3) ^ 2 + qJ(4) ^ 2) + t97; t51 * mrSges(5,2) + (-t27 * mrSges(6,3) - t14) * t87 + (-t26 * mrSges(6,3) + t102) * t83 + m(7) * (t104 * t83 - t132) + m(6) * (t13 * t83 - t132) + m(5) * t31; m(7) * (t36 * t106 - t35 * t87) + m(6) * (t37 * t87 + t38 * t83) + m(5) * t69 + t95; -m(5) * pkin(3) + m(7) * (t57 * t106 - t56 * t87) + m(6) * (t58 * t87 + t59 * t83) + t95; m(5) + m(6) * (t77 + t80) + m(7) * (t117 * t77 + t80); (-m(7) * t11 - t14) * pkin(5) + (m(7) * t104 + t102) * pkin(10) + t150; m(7) * (-pkin(5) * t35 + pkin(10) * t108) - t129 + t130 + (t35 + pkin(5)) * t61 + (t108 - t110) * mrSges(7,3) + t105; m(7) * (-pkin(5) * t56 + pkin(10) * t107) - t127 + t128 + (pkin(5) + t56) * t61 + (t107 - t110) * mrSges(7,3) + t105; -t73 + (m(7) * pkin(10) + mrSges(7,3)) * t106 + (m(7) * pkin(5) - t147) * t87; pkin(5) * t140 + m(7) * (t117 * pkin(10) ^ 2 + pkin(5) ^ 2) + 0.2e1 * pkin(10) * t109 - t105; mrSges(7,1) * t2 - mrSges(7,2) * t3 - Ifges(7,6) * t131 + t118; -t103 * t36 - t62; -t103 * t57 - t62; -t103 * t83; -t103 * pkin(10) + t62; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
