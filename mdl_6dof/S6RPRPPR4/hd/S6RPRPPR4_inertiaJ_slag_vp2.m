% Calculate joint inertia matrix for
% S6RPRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
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
% Datum: 2019-03-09 02:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPPR4_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR4_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR4_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR4_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR4_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR4_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:46:48
% EndTime: 2019-03-09 02:46:50
% DurationCPUTime: 0.96s
% Computational Cost: add. (1302->229), mult. (2517->317), div. (0->0), fcn. (2692->8), ass. (0->92)
t123 = Ifges(6,4) + Ifges(5,5);
t81 = sin(pkin(10));
t83 = cos(pkin(10));
t125 = t81 ^ 2 + t83 ^ 2;
t124 = 0.2e1 * t125;
t122 = Ifges(5,6) - Ifges(6,6);
t111 = Ifges(6,5) * t81;
t113 = Ifges(5,4) * t81;
t114 = cos(qJ(3));
t82 = sin(pkin(9));
t84 = cos(pkin(9));
t86 = sin(qJ(3));
t57 = -t114 * t84 + t82 * t86;
t59 = t114 * t82 + t84 * t86;
t121 = (t111 - t113 + (Ifges(5,1) + Ifges(6,1)) * t83) * t59 + t123 * t57;
t106 = pkin(7) + qJ(2);
t65 = t106 * t84;
t97 = t106 * t82;
t39 = t114 * t97 + t65 * t86;
t120 = t39 ^ 2;
t80 = t84 ^ 2;
t119 = 0.2e1 * t39;
t72 = -pkin(2) * t84 - pkin(1);
t118 = 0.2e1 * t72;
t115 = pkin(4) + pkin(5);
t112 = Ifges(5,4) * t83;
t110 = Ifges(6,5) * t83;
t109 = t59 * t81;
t108 = t59 * t83;
t107 = t81 * mrSges(6,3);
t105 = -pkin(8) + qJ(4);
t29 = -mrSges(5,2) * t57 - mrSges(5,3) * t109;
t32 = -mrSges(6,2) * t109 + mrSges(6,3) * t57;
t104 = t29 + t32;
t30 = mrSges(5,1) * t57 - mrSges(5,3) * t108;
t31 = -t57 * mrSges(6,1) + mrSges(6,2) * t108;
t103 = t30 - t31;
t28 = pkin(3) * t57 - qJ(4) * t59 + t72;
t42 = t114 * t65 - t86 * t97;
t13 = t28 * t81 + t42 * t83;
t27 = mrSges(5,1) * t109 + mrSges(5,2) * t108;
t85 = sin(qJ(6));
t87 = cos(qJ(6));
t58 = t81 * t87 - t83 * t85;
t90 = t81 * t85 + t83 * t87;
t102 = Ifges(7,5) * t58 - Ifges(7,6) * t90;
t101 = t125 * qJ(4) ^ 2;
t99 = t82 ^ 2 + t80;
t7 = qJ(5) * t57 + t13;
t23 = t58 * t59;
t24 = t90 * t59;
t98 = -Ifges(7,5) * t24 - Ifges(7,6) * t23 + Ifges(7,3) * t57;
t96 = -t84 * mrSges(3,1) + mrSges(3,2) * t82;
t35 = mrSges(7,1) * t90 + t58 * mrSges(7,2);
t95 = qJ(5) * t81 + pkin(3);
t33 = t81 * t42;
t12 = t28 * t83 - t33;
t26 = mrSges(6,1) * t109 - mrSges(6,3) * t108;
t93 = -qJ(5) * t108 + t39;
t8 = -pkin(4) * t57 - t12;
t92 = t7 * t83 + t8 * t81;
t9 = -t23 * mrSges(7,1) + t24 * mrSges(7,2);
t91 = -t12 * t81 + t13 * t83;
t75 = t81 * mrSges(5,2);
t69 = Ifges(5,1) * t81 + t112;
t68 = Ifges(6,1) * t81 - t110;
t67 = Ifges(5,2) * t83 + t113;
t66 = -Ifges(6,3) * t83 + t111;
t64 = t105 * t83;
t63 = -t83 * mrSges(5,1) + t75;
t62 = -t83 * mrSges(6,1) - t107;
t61 = t105 * t81;
t60 = -pkin(4) * t83 - t95;
t51 = t59 * mrSges(4,2);
t48 = t115 * t83 + t95;
t41 = t61 * t85 + t64 * t87;
t38 = t61 * t87 - t64 * t85;
t37 = Ifges(7,1) * t58 - Ifges(7,4) * t90;
t36 = Ifges(7,4) * t58 - Ifges(7,2) * t90;
t18 = Ifges(5,6) * t57 + (-Ifges(5,2) * t81 + t112) * t59;
t17 = Ifges(6,6) * t57 + (Ifges(6,3) * t81 + t110) * t59;
t16 = -mrSges(7,1) * t57 - mrSges(7,3) * t24;
t15 = mrSges(7,2) * t57 + mrSges(7,3) * t23;
t14 = pkin(4) * t109 + t93;
t10 = t109 * t115 + t93;
t6 = Ifges(7,1) * t24 + Ifges(7,4) * t23 - Ifges(7,5) * t57;
t5 = Ifges(7,4) * t24 + Ifges(7,2) * t23 - Ifges(7,6) * t57;
t4 = pkin(8) * t109 + t7;
t3 = t33 + (-pkin(8) * t59 - t28) * t83 - t115 * t57;
t2 = t3 * t85 + t4 * t87;
t1 = t3 * t87 - t4 * t85;
t11 = [Ifges(3,2) * t80 - 0.2e1 * pkin(1) * t96 + t51 * t118 + 0.2e1 * t8 * t31 + 0.2e1 * t7 * t32 + t27 * t119 + t23 * t5 + t24 * t6 + 0.2e1 * t14 * t26 + 0.2e1 * t13 * t29 + 0.2e1 * t12 * t30 - 0.2e1 * t10 * t9 + 0.2e1 * t2 * t15 + 0.2e1 * t1 * t16 + Ifges(2,3) + (Ifges(3,1) * t82 + 0.2e1 * Ifges(3,4) * t84) * t82 + 0.2e1 * t99 * qJ(2) * mrSges(3,3) + (mrSges(4,1) * t118 - 0.2e1 * t42 * mrSges(4,3) + (Ifges(5,3) + Ifges(6,2) + Ifges(4,2)) * t57 + t98) * t57 + m(3) * (qJ(2) ^ 2 * t99 + pkin(1) ^ 2) + m(4) * (t42 ^ 2 + t72 ^ 2 + t120) + m(5) * (t12 ^ 2 + t13 ^ 2 + t120) + m(6) * (t14 ^ 2 + t7 ^ 2 + t8 ^ 2) + m(7) * (t1 ^ 2 + t10 ^ 2 + t2 ^ 2) + (mrSges(4,3) * t119 + Ifges(4,1) * t59 + t121 * t83 + (t17 - t18) * t81 + (-t122 * t81 + t123 * t83 - (2 * Ifges(4,4))) * t57) * t59; -m(3) * pkin(1) + t57 * mrSges(4,1) + t58 * t15 - t90 * t16 + t51 + t103 * t83 + t104 * t81 + m(7) * (-t1 * t90 + t2 * t58) + m(6) * (t7 * t81 - t8 * t83) + m(5) * (t12 * t83 + t13 * t81) + m(4) * t72 + t96; m(3) + m(4) + m(7) * (t58 ^ 2 + t90 ^ 2) + (m(5) / 0.2e1 + m(6) / 0.2e1) * t124; (-t1 * t58 - t2 * t90) * mrSges(7,3) - t90 * t5 / 0.2e1 + m(7) * (t1 * t38 - t10 * t48 + t2 * t41) + t121 * t81 / 0.2e1 + (t122 * t83 + t123 * t81) * t57 / 0.2e1 + (Ifges(4,5) + (t66 / 0.2e1 - t67 / 0.2e1) * t81) * t59 + (-t102 / 0.2e1 - Ifges(4,6)) * t57 + t91 * mrSges(5,3) + m(5) * (-pkin(3) * t39 + qJ(4) * t91) + t92 * mrSges(6,2) + m(6) * (qJ(4) * t92 + t14 * t60) - t103 * t81 * qJ(4) + t60 * t26 + t14 * t62 + t58 * t6 / 0.2e1 + t48 * t9 - t10 * t35 + t23 * t36 / 0.2e1 + t24 * t37 / 0.2e1 + t38 * t16 + t41 * t15 - t42 * mrSges(4,2) - pkin(3) * t27 + (t104 * qJ(4) - t17 / 0.2e1 + t18 / 0.2e1 + (t68 / 0.2e1 + t69 / 0.2e1) * t59) * t83 + (t63 - mrSges(4,1)) * t39; m(7) * (-t38 * t90 + t41 * t58); -0.2e1 * pkin(3) * t63 + 0.2e1 * t48 * t35 - t90 * t36 + t58 * t37 + 0.2e1 * t60 * t62 + Ifges(4,3) + (-t66 + t67) * t83 + (t68 + t69) * t81 + 0.2e1 * (-t38 * t58 - t41 * t90) * mrSges(7,3) + m(7) * (t38 ^ 2 + t41 ^ 2 + t48 ^ 2) + m(5) * (pkin(3) ^ 2 + t101) + m(6) * (t60 ^ 2 + t101) + (mrSges(6,2) + mrSges(5,3)) * qJ(4) * t124; m(5) * t39 + m(6) * t14 + m(7) * t10 + t26 + t27 - t9; 0; -m(5) * pkin(3) - t107 + t75 + (-mrSges(6,1) - mrSges(5,1)) * t83 + m(6) * t60 - m(7) * t48 - t35; m(5) + m(6) + m(7); t85 * t15 + t87 * t16 + m(7) * (t1 * t87 + t2 * t85) + m(6) * t8 + t31; -m(6) * t83 + m(7) * (t58 * t85 - t87 * t90); m(7) * (t38 * t87 + t41 * t85) + (m(6) * qJ(4) + mrSges(6,2)) * t81 + (-t58 * t87 - t85 * t90) * mrSges(7,3); 0; m(6) + m(7) * (t85 ^ 2 + t87 ^ 2); mrSges(7,1) * t1 - mrSges(7,2) * t2 - t98; -t35; mrSges(7,1) * t38 - mrSges(7,2) * t41 + t102; 0; mrSges(7,1) * t87 - mrSges(7,2) * t85; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t11(1) t11(2) t11(4) t11(7) t11(11) t11(16); t11(2) t11(3) t11(5) t11(8) t11(12) t11(17); t11(4) t11(5) t11(6) t11(9) t11(13) t11(18); t11(7) t11(8) t11(9) t11(10) t11(14) t11(19); t11(11) t11(12) t11(13) t11(14) t11(15) t11(20); t11(16) t11(17) t11(18) t11(19) t11(20) t11(21);];
Mq  = res;
