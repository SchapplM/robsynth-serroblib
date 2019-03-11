% Calculate joint inertia matrix for
% S6RRRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
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
% Datum: 2019-03-09 15:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPPR2_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR2_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR2_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR2_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR2_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR2_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:24:59
% EndTime: 2019-03-09 15:25:01
% DurationCPUTime: 0.83s
% Computational Cost: add. (1790->215), mult. (3274->291), div. (0->0), fcn. (3603->8), ass. (0->85)
t123 = mrSges(6,2) - mrSges(5,1);
t77 = sin(qJ(6));
t80 = cos(qJ(6));
t57 = mrSges(7,1) * t77 + mrSges(7,2) * t80;
t122 = mrSges(6,3) + t57;
t103 = t77 ^ 2 + t80 ^ 2;
t100 = t103 * mrSges(7,3);
t78 = sin(qJ(3));
t115 = pkin(2) * t78;
t81 = cos(qJ(3));
t68 = pkin(2) * t81 + pkin(3);
t75 = sin(pkin(10));
t76 = cos(pkin(10));
t47 = t115 * t76 + t68 * t75;
t43 = qJ(5) + t47;
t121 = t43 ^ 2;
t65 = pkin(3) * t75 + qJ(5);
t120 = t65 ^ 2;
t119 = 2 * mrSges(6,2);
t79 = sin(qJ(2));
t82 = cos(qJ(2));
t52 = -t78 * t79 + t81 * t82;
t53 = t78 * t82 + t79 * t81;
t35 = -t76 * t52 + t53 * t75;
t36 = t52 * t75 + t53 * t76;
t69 = -pkin(2) * t82 - pkin(1);
t40 = -pkin(3) * t52 + t69;
t87 = -qJ(5) * t36 + t40;
t11 = pkin(4) * t35 + t87;
t118 = -0.2e1 * t11;
t117 = 0.2e1 * t40;
t116 = -pkin(8) - pkin(7);
t56 = m(7) * t103;
t114 = m(6) + t56;
t113 = Ifges(7,4) * t77;
t112 = Ifges(7,4) * t80;
t111 = t35 * t77;
t110 = t35 * t80;
t109 = t43 * t65;
t46 = -t75 * t115 + t68 * t76;
t108 = t46 * mrSges(5,1);
t107 = t47 * mrSges(5,2);
t106 = mrSges(6,1) + mrSges(5,3);
t60 = t116 * t79;
t61 = t116 * t82;
t39 = t78 * t60 - t81 * t61;
t104 = t79 ^ 2 + t82 ^ 2;
t102 = Ifges(7,5) * t111 + Ifges(7,6) * t110 + Ifges(7,3) * t36;
t25 = qJ(4) * t52 + t39;
t38 = t60 * t81 + t61 * t78;
t88 = -qJ(4) * t53 + t38;
t14 = t25 * t75 - t76 * t88;
t16 = t76 * t25 + t75 * t88;
t101 = t14 ^ 2 + t16 ^ 2;
t67 = -pkin(3) * t76 - pkin(4);
t64 = -pkin(9) + t67;
t99 = t103 * t64;
t98 = Ifges(7,5) * t80 - Ifges(7,6) * t77;
t97 = 0.2e1 * t122;
t45 = -pkin(4) - t46;
t4 = (pkin(4) + pkin(9)) * t35 + t87;
t5 = pkin(5) * t36 + t14;
t1 = -t4 * t77 + t5 * t80;
t2 = t4 * t80 + t5 * t77;
t96 = t1 * t80 + t2 * t77;
t95 = mrSges(6,2) - t100;
t94 = t76 * mrSges(5,1) - t75 * mrSges(5,2);
t93 = mrSges(7,1) * t80 - t77 * mrSges(7,2);
t19 = mrSges(7,1) * t36 - mrSges(7,3) * t111;
t20 = -mrSges(7,2) * t36 + mrSges(7,3) * t110;
t92 = t80 * t19 + t77 * t20;
t91 = -0.2e1 * t100;
t58 = -Ifges(7,2) * t77 + t112;
t59 = Ifges(7,1) * t80 - t113;
t90 = -t77 * t58 + t80 * t59 + Ifges(6,1) + Ifges(4,3) + Ifges(5,3);
t89 = (mrSges(4,1) * t81 - mrSges(4,2) * t78) * pkin(2);
t10 = Ifges(7,5) * t36 + (Ifges(7,1) * t77 + t112) * t35;
t6 = -t35 * pkin(5) + t16;
t9 = Ifges(7,6) * t36 + (Ifges(7,2) * t80 + t113) * t35;
t86 = t59 * t111 / 0.2e1 + t58 * t110 / 0.2e1 + t6 * t57 + t38 * mrSges(4,1) + Ifges(4,6) * t52 + Ifges(4,5) * t53 - t77 * t9 / 0.2e1 + t80 * t10 / 0.2e1 + (-mrSges(5,2) + mrSges(6,3)) * t16 - t96 * mrSges(7,3) - t39 * mrSges(4,2) + (Ifges(6,5) - Ifges(5,6)) * t35 + t123 * t14 + (t98 / 0.2e1 + Ifges(5,5) - Ifges(6,4)) * t36;
t42 = -pkin(9) + t45;
t29 = t36 * mrSges(6,3);
t28 = t36 * mrSges(5,2);
t18 = t93 * t35;
t3 = [t82 * (Ifges(3,4) * t79 + Ifges(3,2) * t82) + t79 * (Ifges(3,1) * t79 + Ifges(3,4) * t82) - 0.2e1 * pkin(1) * (-t82 * mrSges(3,1) + t79 * mrSges(3,2)) + 0.2e1 * t69 * (-mrSges(4,1) * t52 + mrSges(4,2) * t53) + t53 * (Ifges(4,1) * t53 + Ifges(4,4) * t52) + t52 * (Ifges(4,4) * t53 + Ifges(4,2) * t52) + t28 * t117 + t29 * t118 - 0.2e1 * t6 * t18 + 0.2e1 * t1 * t19 + 0.2e1 * t2 * t20 + Ifges(2,3) + 0.2e1 * (-t38 * t53 + t39 * t52) * mrSges(4,3) + 0.2e1 * t104 * pkin(7) * mrSges(3,3) + (mrSges(5,1) * t117 + mrSges(6,2) * t118 + t77 * t10 + t80 * t9 + (Ifges(6,3) + Ifges(5,2)) * t35 - 0.2e1 * t106 * t16) * t35 + m(3) * (pkin(7) ^ 2 * t104 + pkin(1) ^ 2) + m(4) * (t38 ^ 2 + t39 ^ 2 + t69 ^ 2) + m(5) * (t40 ^ 2 + t101) + m(7) * (t1 ^ 2 + t2 ^ 2 + t6 ^ 2) + m(6) * (t11 ^ 2 + t101) + ((Ifges(5,1) + Ifges(6,2)) * t36 + t102 + 0.2e1 * (-Ifges(5,4) - Ifges(6,6)) * t35 + 0.2e1 * t106 * t14) * t36; t86 + (-t79 * mrSges(3,1) - t82 * mrSges(3,2)) * pkin(7) + (-t47 * t35 - t46 * t36) * mrSges(5,3) + (-t43 * t35 + t45 * t36) * mrSges(6,1) + Ifges(3,5) * t79 + Ifges(3,6) * t82 + m(7) * (t42 * t96 + t43 * t6) - t43 * t18 + t92 * t42 + m(6) * (t14 * t45 + t16 * t43) + m(5) * (-t14 * t46 + t16 * t47) + (m(4) * (t38 * t81 + t39 * t78) + (t52 * t78 - t53 * t81) * mrSges(4,3)) * pkin(2); 0.2e1 * t108 - 0.2e1 * t107 + t45 * t119 + Ifges(3,3) + t43 * t97 + 0.2e1 * t89 + t42 * t91 + m(7) * (t103 * t42 ^ 2 + t121) + m(6) * (t45 ^ 2 + t121) + m(5) * (t46 ^ 2 + t47 ^ 2) + m(4) * (t78 ^ 2 + t81 ^ 2) * pkin(2) ^ 2 + t90; t86 + t92 * t64 + m(6) * (t14 * t67 + t16 * t65) + (-t35 * t65 + t36 * t67) * mrSges(6,1) + (m(5) * (-t14 * t76 + t16 * t75) + (-t35 * t75 - t36 * t76) * mrSges(5,3)) * pkin(3) - t65 * t18 + m(7) * (t6 * t65 + t64 * t96); t108 - t107 + t89 + (t45 + t67) * mrSges(6,2) + m(7) * (t42 * t99 + t109) + m(6) * (t45 * t67 + t109) + (m(5) * (t46 * t76 + t47 * t75) + t94) * pkin(3) + t90 + t122 * (t43 + t65) + (-t42 - t64) * t100; t67 * t119 + t65 * t97 + t64 * t91 + m(7) * (t103 * t64 ^ 2 + t120) + m(6) * (t67 ^ 2 + t120) + t90 + (0.2e1 * t94 + m(5) * (t75 ^ 2 + t76 ^ 2) * pkin(3)) * pkin(3); -t77 * t19 + t80 * t20 + t28 - t29 - t123 * t35 + m(7) * (-t1 * t77 + t2 * t80) + m(6) * t11 + m(5) * t40; 0; 0; m(5) + t114; m(6) * t14 + m(7) * t96 + t36 * mrSges(6,1) + t92; m(6) * t45 + t42 * t56 + t95; m(6) * t67 + m(7) * t99 + t95; 0; t114; mrSges(7,1) * t1 - mrSges(7,2) * t2 + t102; t42 * t93 + t98; t64 * t93 + t98; -t57; t93; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t3(1) t3(2) t3(4) t3(7) t3(11) t3(16); t3(2) t3(3) t3(5) t3(8) t3(12) t3(17); t3(4) t3(5) t3(6) t3(9) t3(13) t3(18); t3(7) t3(8) t3(9) t3(10) t3(14) t3(19); t3(11) t3(12) t3(13) t3(14) t3(15) t3(20); t3(16) t3(17) t3(18) t3(19) t3(20) t3(21);];
Mq  = res;
