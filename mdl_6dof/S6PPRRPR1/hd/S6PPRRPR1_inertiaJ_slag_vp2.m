% Calculate joint inertia matrix for
% S6PPRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2,theta5]';
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
% Datum: 2019-03-08 18:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PPRRPR1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR1_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRPR1_inertiaJ_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRPR1_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRPR1_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRRPR1_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:45:36
% EndTime: 2019-03-08 18:45:38
% DurationCPUTime: 0.83s
% Computational Cost: add. (1095->242), mult. (2819->363), div. (0->0), fcn. (3164->14), ass. (0->104)
t123 = m(6) + m(7);
t124 = 2 * pkin(9);
t80 = sin(pkin(13));
t84 = cos(pkin(13));
t88 = sin(qJ(6));
t91 = cos(qJ(6));
t53 = t80 * t91 + t84 * t88;
t89 = sin(qJ(4));
t39 = t53 * t89;
t52 = -t80 * t88 + t84 * t91;
t40 = t52 * t89;
t16 = t39 * mrSges(7,1) + t40 * mrSges(7,2);
t109 = t84 * t89;
t112 = t80 * t89;
t46 = mrSges(6,1) * t112 + mrSges(6,2) * t109;
t122 = t16 + t46;
t85 = cos(pkin(12));
t86 = cos(pkin(7));
t108 = t85 * t86;
t82 = sin(pkin(7));
t90 = sin(qJ(3));
t111 = t82 * t90;
t81 = sin(pkin(12));
t83 = sin(pkin(6));
t87 = cos(pkin(6));
t93 = cos(qJ(3));
t20 = t87 * t111 + (t90 * t108 + t81 * t93) * t83;
t41 = -t82 * t83 * t85 + t86 * t87;
t92 = cos(qJ(4));
t11 = t20 * t89 - t41 * t92;
t10 = t11 ^ 2;
t110 = t82 * t93;
t18 = -t87 * t110 + (-t108 * t93 + t81 * t90) * t83;
t121 = t18 ^ 2;
t43 = t89 * t111 - t92 * t86;
t42 = t43 ^ 2;
t120 = m(6) * pkin(4);
t119 = -t80 / 0.2e1;
t118 = t84 / 0.2e1;
t117 = pkin(9) * t92;
t116 = Ifges(6,4) * t80;
t115 = Ifges(6,4) * t84;
t114 = t11 * t89;
t5 = t43 * t11;
t113 = t43 * t89;
t107 = t92 * mrSges(5,3);
t64 = -mrSges(5,1) * t92 + mrSges(5,2) * t89;
t106 = mrSges(4,1) - t64;
t105 = pkin(10) + qJ(5);
t104 = -Ifges(7,5) * t40 + Ifges(7,6) * t39;
t60 = -t84 * mrSges(6,1) + t80 * mrSges(6,2);
t103 = t60 - mrSges(5,1);
t58 = -pkin(4) * t92 - qJ(5) * t89 - pkin(3);
t32 = t84 * t117 + t80 * t58;
t102 = t80 ^ 2 + t84 ^ 2;
t21 = -t52 * mrSges(7,1) + t53 * mrSges(7,2);
t100 = t21 + t103;
t99 = t89 * mrSges(5,3) + t122;
t13 = t20 * t92 + t41 * t89;
t3 = -t13 * t80 + t18 * t84;
t4 = t13 * t84 + t18 * t80;
t98 = -t3 * t80 + t4 * t84;
t45 = t92 * t111 + t86 * t89;
t24 = -t84 * t110 - t45 * t80;
t25 = -t80 * t110 + t45 * t84;
t97 = -t24 * t80 + t25 * t84;
t51 = t84 * t58;
t31 = -t80 * t117 + t51;
t96 = -t31 * t80 + t32 * t84;
t95 = pkin(9) ^ 2;
t79 = t92 ^ 2;
t78 = t89 ^ 2;
t75 = t82 ^ 2;
t73 = t78 * t95;
t71 = -pkin(5) * t84 - pkin(4);
t69 = t75 * t93 ^ 2;
t63 = Ifges(6,1) * t80 + t115;
t62 = Ifges(6,2) * t84 + t116;
t61 = t105 * t84;
t59 = t105 * t80;
t57 = (pkin(5) * t80 + pkin(9)) * t89;
t55 = -mrSges(6,1) * t92 - mrSges(6,3) * t109;
t54 = mrSges(6,2) * t92 - mrSges(6,3) * t112;
t49 = Ifges(7,5) * t53;
t48 = Ifges(7,6) * t52;
t38 = -Ifges(6,5) * t92 + (Ifges(6,1) * t84 - t116) * t89;
t37 = -Ifges(6,6) * t92 + (-Ifges(6,2) * t80 + t115) * t89;
t30 = -mrSges(7,1) * t92 - mrSges(7,3) * t40;
t29 = mrSges(7,2) * t92 - mrSges(7,3) * t39;
t28 = -t59 * t88 + t61 * t91;
t27 = -t59 * t91 - t61 * t88;
t26 = -pkin(10) * t112 + t32;
t23 = Ifges(7,1) * t53 + Ifges(7,4) * t52;
t22 = Ifges(7,4) * t53 + Ifges(7,2) * t52;
t17 = -pkin(10) * t109 + t51 + (-pkin(9) * t80 - pkin(5)) * t92;
t15 = Ifges(7,1) * t40 - Ifges(7,4) * t39 - Ifges(7,5) * t92;
t14 = Ifges(7,4) * t40 - Ifges(7,2) * t39 - Ifges(7,6) * t92;
t9 = t24 * t88 + t25 * t91;
t8 = t24 * t91 - t25 * t88;
t7 = t17 * t88 + t26 * t91;
t6 = t17 * t91 - t26 * t88;
t2 = t3 * t88 + t4 * t91;
t1 = t3 * t91 - t4 * t88;
t12 = [m(2) + m(7) * (t1 ^ 2 + t2 ^ 2 + t10) + m(6) * (t3 ^ 2 + t4 ^ 2 + t10) + m(5) * (t13 ^ 2 + t10 + t121) + m(4) * (t20 ^ 2 + t41 ^ 2 + t121) + m(3) * (t87 ^ 2 + (t81 ^ 2 + t85 ^ 2) * t83 ^ 2); m(3) * t87 + m(7) * (t1 * t8 + t2 * t9 + t5) + m(6) * (t24 * t3 + t25 * t4 + t5) + m(5) * (-t18 * t110 + t13 * t45 + t5) + m(4) * (t86 * t41 + (-t18 * t93 + t20 * t90) * t82); m(3) + m(7) * (t8 ^ 2 + t9 ^ 2 + t42) + m(5) * (t45 ^ 2 + t42 + t69) + m(6) * (t24 ^ 2 + t25 ^ 2 + t42) + m(4) * (t75 * t90 ^ 2 + t86 ^ 2 + t69); t13 * t107 - t20 * mrSges(4,2) + t1 * t30 + t2 * t29 + t3 * t55 + t4 * t54 - t106 * t18 + t99 * t11 + m(7) * (t1 * t6 + t11 * t57 + t2 * t7) + m(6) * (pkin(9) * t114 + t3 * t31 + t32 * t4) + m(5) * (-pkin(3) * t18 + (t13 * t92 + t114) * pkin(9)); t45 * t107 + t24 * t55 + t25 * t54 + t9 * t29 + t8 * t30 + (-t90 * mrSges(4,2) + t106 * t93) * t82 + t99 * t43 + m(7) * (t43 * t57 + t6 * t8 + t7 * t9) + m(5) * (pkin(3) * t110 + (t45 * t92 + t113) * pkin(9)) + m(6) * (pkin(9) * t113 + t24 * t31 + t25 * t32); -0.2e1 * pkin(3) * t64 - t39 * t14 + t40 * t15 + 0.2e1 * t57 * t16 + 0.2e1 * t7 * t29 + 0.2e1 * t6 * t30 + 0.2e1 * t31 * t55 + 0.2e1 * t32 * t54 + Ifges(4,3) + (t78 + t79) * mrSges(5,3) * t124 + m(7) * (t57 ^ 2 + t6 ^ 2 + t7 ^ 2) + m(6) * (t31 ^ 2 + t32 ^ 2 + t73) + m(5) * (pkin(3) ^ 2 + t79 * t95 + t73) + ((Ifges(7,3) + Ifges(6,3) + Ifges(5,2)) * t92 + t104) * t92 + (Ifges(5,1) * t89 + t46 * t124 - t80 * t37 + t84 * t38 + (-Ifges(6,5) * t84 + Ifges(6,6) * t80 + (2 * Ifges(5,4))) * t92) * t89; -t13 * mrSges(5,2) + (-t1 * t53 + t2 * t52) * mrSges(7,3) + t98 * mrSges(6,3) + t100 * t11 + m(7) * (t1 * t27 + t11 * t71 + t2 * t28) + m(6) * (-pkin(4) * t11 + t98 * qJ(5)); -t45 * mrSges(5,2) + (t52 * t9 - t53 * t8) * mrSges(7,3) + t97 * mrSges(6,3) + t100 * t43 + m(7) * (t27 * t8 + t28 * t9 + t43 * t71) + m(6) * (-pkin(4) * t43 + t97 * qJ(5)); m(7) * (t27 * t6 + t28 * t7 + t57 * t71) + t37 * t118 + t80 * t38 / 0.2e1 + t57 * t21 + t71 * t16 - pkin(4) * t46 + t52 * t14 / 0.2e1 + t53 * t15 / 0.2e1 - t39 * t22 / 0.2e1 + t40 * t23 / 0.2e1 + t28 * t29 + t27 * t30 + (-pkin(9) * mrSges(5,2) + Ifges(6,5) * t119 - Ifges(6,6) * t84 / 0.2e1 - t49 / 0.2e1 - t48 / 0.2e1 + Ifges(5,6)) * t92 + (t52 * t7 - t53 * t6) * mrSges(7,3) + t96 * mrSges(6,3) + (m(6) * t96 + t84 * t54 - t80 * t55) * qJ(5) + (t63 * t118 + t62 * t119 + Ifges(5,5) + (t103 - t120) * pkin(9)) * t89; -0.2e1 * pkin(4) * t60 + 0.2e1 * t71 * t21 + t52 * t22 + t53 * t23 + t84 * t62 + t80 * t63 + Ifges(5,3) + m(7) * (t27 ^ 2 + t28 ^ 2 + t71 ^ 2) + m(6) * (t102 * qJ(5) ^ 2 + pkin(4) ^ 2) + 0.2e1 * (-t27 * t53 + t28 * t52) * mrSges(7,3) + 0.2e1 * t102 * qJ(5) * mrSges(6,3); t11 * t123; t43 * t123; m(6) * t89 * pkin(9) + m(7) * t57 + t122; m(7) * t71 - t120 + t21 + t60; t123; mrSges(7,1) * t1 - mrSges(7,2) * t2; mrSges(7,1) * t8 - mrSges(7,2) * t9; mrSges(7,1) * t6 - mrSges(7,2) * t7 - Ifges(7,3) * t92 - t104; mrSges(7,1) * t27 - mrSges(7,2) * t28 + t48 + t49; 0; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t12(1) t12(2) t12(4) t12(7) t12(11) t12(16); t12(2) t12(3) t12(5) t12(8) t12(12) t12(17); t12(4) t12(5) t12(6) t12(9) t12(13) t12(18); t12(7) t12(8) t12(9) t12(10) t12(14) t12(19); t12(11) t12(12) t12(13) t12(14) t12(15) t12(20); t12(16) t12(17) t12(18) t12(19) t12(20) t12(21);];
Mq  = res;
