% Calculate joint inertia matrix for
% S5RRRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRR7_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR7_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR7_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR7_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR7_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR7_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:21:24
% EndTime: 2019-12-31 22:21:25
% DurationCPUTime: 0.65s
% Computational Cost: add. (1378->177), mult. (2619->269), div. (0->0), fcn. (2776->8), ass. (0->89)
t72 = sin(qJ(5));
t104 = t72 * mrSges(6,3);
t74 = sin(qJ(3));
t75 = sin(qJ(2));
t78 = cos(qJ(3));
t79 = cos(qJ(2));
t47 = -t74 * t75 + t78 * t79;
t48 = t74 * t79 + t78 * t75;
t73 = sin(qJ(4));
t77 = cos(qJ(4));
t26 = -t77 * t47 + t73 * t48;
t27 = t73 * t47 + t77 * t48;
t15 = -t26 * mrSges(6,2) - t27 * t104;
t76 = cos(qJ(5));
t106 = t27 * t76;
t16 = t26 * mrSges(6,1) - mrSges(6,3) * t106;
t124 = t76 * t15 - t72 * t16;
t53 = -t76 * mrSges(6,1) + t72 * mrSges(6,2);
t114 = pkin(4) * t53;
t68 = t72 ^ 2;
t111 = mrSges(6,3) * t68;
t63 = pkin(9) * t111;
t70 = t76 ^ 2;
t110 = mrSges(6,3) * t70;
t64 = pkin(9) * t110;
t123 = t63 + t64 - t114;
t116 = -pkin(7) - pkin(6);
t93 = t116 * t75;
t94 = t116 * t79;
t31 = t74 * t93 - t78 * t94;
t21 = t47 * pkin(8) + t31;
t30 = t74 * t94 + t78 * t93;
t86 = -t48 * pkin(8) + t30;
t11 = t73 * t21 - t77 * t86;
t89 = mrSges(6,1) * t72 + mrSges(6,2) * t76;
t14 = t89 * t27;
t122 = m(6) * t11 + t14;
t112 = t77 * pkin(3);
t60 = -pkin(4) - t112;
t42 = t60 * t53;
t59 = t73 * pkin(3) + pkin(9);
t51 = t59 * t111;
t52 = t59 * t110;
t65 = mrSges(5,1) * t112;
t121 = t42 + t51 + t52 + t65;
t62 = -t79 * pkin(2) - pkin(1);
t34 = -t47 * pkin(3) + t62;
t10 = t26 * pkin(4) - t27 * pkin(9) + t34;
t13 = t77 * t21 + t73 * t86;
t3 = t72 * t10 + t76 * t13;
t113 = t3 * t76;
t2 = t76 * t10 - t72 * t13;
t90 = -t2 * t72 + t113;
t120 = m(6) * t90 + t124;
t119 = t11 ^ 2;
t118 = 0.2e1 * t11;
t117 = 0.2e1 * t34;
t115 = pkin(2) * t74;
t109 = Ifges(6,4) * t72;
t108 = Ifges(6,4) * t76;
t107 = t27 * t72;
t61 = t78 * pkin(2) + pkin(3);
t40 = t77 * t115 + t73 * t61;
t105 = t40 * mrSges(5,2);
t102 = t73 * mrSges(5,2);
t100 = Ifges(6,5) * t106 + Ifges(6,3) * t26;
t99 = Ifges(6,5) * t72 + Ifges(6,6) * t76;
t98 = t68 + t70;
t97 = t75 ^ 2 + t79 ^ 2;
t96 = pkin(3) * t102;
t54 = Ifges(6,2) * t76 + t109;
t55 = Ifges(6,1) * t72 + t108;
t95 = t76 * t54 + t72 * t55 + Ifges(5,3);
t92 = t98 * t59;
t91 = Ifges(4,3) + t95;
t39 = -t73 * t115 + t77 * t61;
t88 = (t78 * mrSges(4,1) - t74 * mrSges(4,2)) * pkin(2);
t37 = -pkin(4) - t39;
t28 = t37 * t53;
t38 = pkin(9) + t40;
t32 = t38 * t111;
t33 = t38 * t110;
t35 = t39 * mrSges(5,1);
t87 = t28 + t32 + t33 + t35 + t95 - t105;
t6 = Ifges(6,6) * t26 + (-Ifges(6,2) * t72 + t108) * t27;
t7 = Ifges(6,5) * t26 + (Ifges(6,1) * t76 - t109) * t27;
t85 = -t13 * mrSges(5,2) + mrSges(6,3) * t113 - t2 * t104 + t72 * t7 / 0.2e1 - t54 * t107 / 0.2e1 + t55 * t106 / 0.2e1 + Ifges(5,5) * t27 + t76 * t6 / 0.2e1 + (t99 / 0.2e1 - Ifges(5,6)) * t26 + (t53 - mrSges(5,1)) * t11;
t84 = t30 * mrSges(4,1) - t31 * mrSges(4,2) + Ifges(4,5) * t48 + Ifges(4,6) * t47 + t85;
t1 = [t75 * (Ifges(3,1) * t75 + Ifges(3,4) * t79) + t79 * (Ifges(3,4) * t75 + Ifges(3,2) * t79) - 0.2e1 * pkin(1) * (-t79 * mrSges(3,1) + t75 * mrSges(3,2)) + 0.2e1 * t62 * (-t47 * mrSges(4,1) + t48 * mrSges(4,2)) + t47 * (Ifges(4,4) * t48 + Ifges(4,2) * t47) + t48 * (Ifges(4,1) * t48 + Ifges(4,4) * t47) + t14 * t118 + 0.2e1 * t3 * t15 + 0.2e1 * t2 * t16 + Ifges(2,3) + (mrSges(5,1) * t117 - 0.2e1 * t13 * mrSges(5,3) + Ifges(5,2) * t26 + t100) * t26 + (mrSges(5,2) * t117 + mrSges(5,3) * t118 + Ifges(5,1) * t27 - t72 * t6 + t76 * t7 + (-Ifges(6,6) * t72 - (2 * Ifges(5,4))) * t26) * t27 + m(3) * (t97 * pkin(6) ^ 2 + pkin(1) ^ 2) + m(4) * (t30 ^ 2 + t31 ^ 2 + t62 ^ 2) + m(5) * (t13 ^ 2 + t34 ^ 2 + t119) + m(6) * (t2 ^ 2 + t3 ^ 2 + t119) + 0.2e1 * (-t30 * t48 + t31 * t47) * mrSges(4,3) + 0.2e1 * t97 * pkin(6) * mrSges(3,3); (m(4) * (t30 * t78 + t31 * t74) + (t74 * t47 - t78 * t48) * mrSges(4,3)) * pkin(2) + m(6) * (t37 * t11 + t90 * t38) + t124 * t38 + m(5) * (-t39 * t11 + t40 * t13) + (-t75 * mrSges(3,1) - t79 * mrSges(3,2)) * pkin(6) + (-t40 * t26 - t39 * t27) * mrSges(5,3) + t84 + Ifges(3,5) * t75 + Ifges(3,6) * t79 + t37 * t14; -0.2e1 * t105 + Ifges(3,3) + 0.2e1 * t28 + 0.2e1 * t32 + 0.2e1 * t33 + 0.2e1 * t35 + 0.2e1 * t88 + m(6) * (t98 * t38 ^ 2 + t37 ^ 2) + m(5) * (t39 ^ 2 + t40 ^ 2) + m(4) * (t74 ^ 2 + t78 ^ 2) * pkin(2) ^ 2 + t91; (m(5) * (-t11 * t77 + t13 * t73) + (-t73 * t26 - t77 * t27) * mrSges(5,3)) * pkin(3) + t84 + t122 * t60 + t120 * t59; m(6) * (t60 * t37 + t38 * t92) + (m(5) * (t39 * t77 + t40 * t73) - t102) * pkin(3) + t87 + t88 + Ifges(4,3) + t121; -0.2e1 * t96 + 0.2e1 * t42 + 0.2e1 * t51 + 0.2e1 * t52 + 0.2e1 * t65 + m(6) * (t98 * t59 ^ 2 + t60 ^ 2) + m(5) * (t73 ^ 2 + t77 ^ 2) * pkin(3) ^ 2 + t91; -t122 * pkin(4) + t120 * pkin(9) + t85; m(6) * (t98 * t38 * pkin(9) - pkin(4) * t37) + t87 + t123; -t96 + m(6) * (-pkin(4) * t60 + pkin(9) * t92) + t95 + t121 + t123; 0.2e1 * t64 + 0.2e1 * t63 - 0.2e1 * t114 + m(6) * (t98 * pkin(9) ^ 2 + pkin(4) ^ 2) + t95; t2 * mrSges(6,1) - t3 * mrSges(6,2) - Ifges(6,6) * t107 + t100; -t89 * t38 + t99; -t89 * t59 + t99; -t89 * pkin(9) + t99; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
