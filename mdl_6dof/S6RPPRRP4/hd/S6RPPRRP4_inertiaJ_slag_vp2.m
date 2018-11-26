% Calculate joint inertia matrix for
% S6RPPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
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
% Datum: 2018-11-23 15:45
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RPPRRP4_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP4_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP4_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP4_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP4_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRP4_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:45:17
% EndTime: 2018-11-23 15:45:17
% DurationCPUTime: 0.70s
% Computational Cost: add. (710->214), mult. (1228->279), div. (0->0), fcn. (956->6), ass. (0->92)
t117 = Ifges(7,2) + Ifges(6,3);
t67 = sin(qJ(5));
t69 = cos(qJ(5));
t87 = t67 ^ 2 + t69 ^ 2;
t65 = sin(pkin(9));
t66 = cos(pkin(9));
t70 = cos(qJ(4));
t98 = t67 * t70;
t18 = t65 * t98 + t69 * t66;
t103 = t18 * t67;
t94 = t69 * t70;
t20 = t65 * t94 - t67 * t66;
t116 = t20 * t69 + t103;
t115 = m(6) + m(7);
t114 = mrSges(7,2) + mrSges(6,3);
t113 = -m(7) * pkin(5) - mrSges(7,1);
t112 = t114 * t87;
t111 = pkin(5) * t67;
t68 = sin(qJ(4));
t110 = t68 * pkin(8);
t71 = -pkin(1) - pkin(2);
t31 = -t65 * qJ(2) + t66 * t71;
t25 = pkin(3) - t31;
t58 = t70 * pkin(4);
t8 = t25 + t58 + t110;
t109 = t69 * t8;
t32 = t66 * qJ(2) + t65 * t71;
t26 = -pkin(7) + t32;
t4 = t26 * t94 + t67 * t8;
t108 = Ifges(6,4) * t67;
t107 = Ifges(6,4) * t69;
t106 = Ifges(7,5) * t67;
t105 = Ifges(7,5) * t69;
t104 = Ifges(7,6) * t70;
t101 = t26 * t65;
t100 = t65 * t68;
t99 = t67 * t68;
t97 = t68 * t69;
t96 = t69 * mrSges(6,2);
t95 = t69 * mrSges(7,3);
t93 = -Ifges(7,4) - Ifges(6,5);
t27 = -t70 * mrSges(6,2) + mrSges(6,3) * t99;
t30 = mrSges(7,2) * t99 + t70 * mrSges(7,3);
t92 = t27 + t30;
t28 = t70 * mrSges(6,1) + mrSges(6,3) * t97;
t29 = -t70 * mrSges(7,1) - mrSges(7,2) * t97;
t91 = -t28 + t29;
t35 = -t69 * mrSges(6,1) + t67 * mrSges(6,2);
t90 = t35 - mrSges(5,1);
t89 = t87 * t110;
t88 = t87 * pkin(8) ^ 2;
t62 = t68 ^ 2;
t64 = t70 ^ 2;
t86 = t62 + t64;
t85 = qJ(6) * t69;
t84 = Ifges(6,6) * t99 + t117 * t70;
t22 = -mrSges(7,1) * t99 + t68 * t95;
t41 = t68 * t85;
t5 = t41 + (t26 - t111) * t68;
t83 = m(7) * t5 + t22;
t82 = t86 * mrSges(5,3);
t80 = t116 * pkin(8);
t3 = -t26 * t98 + t109;
t78 = -t3 * t67 + t4 * t69;
t77 = -t67 * mrSges(6,1) - t96;
t76 = t85 - t111;
t1 = t70 * qJ(6) + t4;
t2 = -t109 + (t26 * t67 - pkin(5)) * t70;
t73 = m(7) * (t1 * t69 + t2 * t67) + t92 * t69 + t91 * t67;
t60 = t66 ^ 2;
t59 = t65 ^ 2;
t55 = Ifges(7,4) * t67;
t54 = Ifges(6,5) * t67;
t52 = Ifges(6,6) * t69;
t50 = t70 * mrSges(5,1);
t46 = t62 * t59;
t40 = Ifges(6,1) * t67 + t107;
t39 = Ifges(7,1) * t67 - t105;
t38 = Ifges(6,2) * t69 + t108;
t37 = -Ifges(7,3) * t69 + t106;
t36 = -t68 * mrSges(5,2) + t50;
t34 = -t69 * mrSges(7,1) - t67 * mrSges(7,3);
t33 = -t69 * pkin(5) - t67 * qJ(6) - pkin(4);
t24 = t26 ^ 2;
t23 = t77 * t68;
t21 = t62 * t24;
t14 = t62 * t101;
t13 = Ifges(6,5) * t70 + (-Ifges(6,1) * t69 + t108) * t68;
t12 = Ifges(7,4) * t70 + (-Ifges(7,1) * t69 - t106) * t68;
t11 = Ifges(6,6) * t70 + (Ifges(6,2) * t67 - t107) * t68;
t10 = t104 + (-Ifges(7,3) * t67 - t105) * t68;
t6 = [(2 * pkin(1) * mrSges(3,1)) - 0.2e1 * t31 * mrSges(4,1) + 0.2e1 * t32 * mrSges(4,2) + 0.2e1 * qJ(2) * mrSges(3,3) + 0.2e1 * t1 * t30 + 0.2e1 * t2 * t29 + 0.2e1 * t5 * t22 + 0.2e1 * t25 * t36 + 0.2e1 * t4 * t27 + 0.2e1 * t3 * t28 + Ifges(3,2) + Ifges(2,3) + Ifges(4,3) - 0.2e1 * t26 * t82 + (Ifges(5,2) * t70 + t84) * t70 + m(5) * (t64 * t24 + t25 ^ 2 + t21) + m(3) * ((pkin(1) ^ 2) + qJ(2) ^ 2) + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) + m(6) * (t3 ^ 2 + t4 ^ 2 + t21) + m(4) * (t31 ^ 2 + t32 ^ 2) + (Ifges(5,1) * t68 + 0.2e1 * Ifges(5,4) * t70 + 0.2e1 * t26 * t23 + (-t10 + t11 - t104) * t67 + (t93 * t70 - t12 - t13) * t69) * t68; -m(3) * pkin(1) - mrSges(3,1) + (-mrSges(4,1) - t36) * t66 + t92 * t20 + t91 * t18 + (mrSges(4,2) + (t22 + t23) * t68 - t82) * t65 + m(7) * (t20 * t1 + t5 * t100 + t18 * t2) + m(6) * (-t18 * t3 + t20 * t4 + t14) + m(5) * (t64 * t101 - t66 * t25 + t14) + m(4) * (t66 * t31 + t65 * t32); m(3) + m(5) * (t64 * t59 + t46 + t60) + m(4) * (t59 + t60) + t115 * (t18 ^ 2 + t20 ^ 2 + t46); (-t23 - t83) * t70 + (m(6) * (-t26 * t70 + t78) + t73) * t68; t115 * (t20 * t97 + (-t65 * t70 + t103) * t68); m(4) + m(5) * t86 + 0.2e1 * (m(6) / 0.2e1 + m(7) / 0.2e1) * (t87 * t62 + t64); -pkin(4) * t23 + t5 * t34 + t83 * t33 + (t55 / 0.2e1 + t54 / 0.2e1 + t52 / 0.2e1 - Ifges(5,6) - t26 * mrSges(5,2)) * t70 + (-t104 / 0.2e1 - t10 / 0.2e1 + t11 / 0.2e1 + t1 * mrSges(7,2) + t4 * mrSges(6,3)) * t69 + (t12 / 0.2e1 + t13 / 0.2e1 + t2 * mrSges(7,2) - t3 * mrSges(6,3)) * t67 + (m(6) * t78 + t73) * pkin(8) + (-Ifges(5,5) + (-t39 / 0.2e1 - t40 / 0.2e1) * t69 + (-t37 / 0.2e1 + t38 / 0.2e1) * t67 + (-m(6) * pkin(4) + t90) * t26) * t68; (-t70 * mrSges(5,2) + (t34 + t90) * t68) * t65 + m(6) * (-pkin(4) * t100 + t80) + m(7) * (t33 * t100 + t80) + t114 * t116; t50 + (-t34 - t35) * t70 + m(7) * (-t33 * t70 + t89) + m(6) * (t58 + t89) + (-mrSges(5,2) + t112) * t68; -0.2e1 * pkin(4) * t35 + 0.2e1 * t33 * t34 + Ifges(5,3) + (-t37 + t38) * t69 + (t39 + t40) * t67 + m(7) * (t33 ^ 2 + t88) + m(6) * (pkin(4) ^ 2 + t88) + 0.2e1 * pkin(8) * t112; -pkin(5) * t29 + m(7) * (-pkin(5) * t2 + qJ(6) * t1) + qJ(6) * t30 + t1 * mrSges(7,3) - t4 * mrSges(6,2) + t3 * mrSges(6,1) - t2 * mrSges(7,1) + (-Ifges(7,6) * t67 + t93 * t69) * t68 + t84; (m(7) * qJ(6) - mrSges(6,2) + mrSges(7,3)) * t20 + (-mrSges(6,1) + t113) * t18; -t68 * t96 - mrSges(6,1) * t99 + m(7) * (-pkin(5) * t99 + t41) + t22; -Ifges(7,6) * t69 + t52 + t54 + t55 + t76 * mrSges(7,2) + (m(7) * t76 - t67 * mrSges(7,1) + t77 + t95) * pkin(8); 0.2e1 * pkin(5) * mrSges(7,1) + 0.2e1 * qJ(6) * mrSges(7,3) + m(7) * (pkin(5) ^ 2 + qJ(6) ^ 2) + t117; m(7) * t2 + t29; m(7) * t18; m(7) * t99; (m(7) * pkin(8) + mrSges(7,2)) * t67; t113; m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t6(1) t6(2) t6(4) t6(7) t6(11) t6(16); t6(2) t6(3) t6(5) t6(8) t6(12) t6(17); t6(4) t6(5) t6(6) t6(9) t6(13) t6(18); t6(7) t6(8) t6(9) t6(10) t6(14) t6(19); t6(11) t6(12) t6(13) t6(14) t6(15) t6(20); t6(16) t6(17) t6(18) t6(19) t6(20) t6(21);];
Mq  = res;
