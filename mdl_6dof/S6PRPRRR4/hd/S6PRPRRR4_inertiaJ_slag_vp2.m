% Calculate joint inertia matrix for
% S6PRPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
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
% Datum: 2019-03-08 20:39
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPRRR4_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR4_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR4_inertiaJ_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR4_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR4_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRR4_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:36:05
% EndTime: 2019-03-08 20:36:06
% DurationCPUTime: 0.77s
% Computational Cost: add. (1490->241), mult. (3143->359), div. (0->0), fcn. (3543->12), ass. (0->95)
t116 = cos(qJ(4));
t80 = sin(pkin(12));
t82 = cos(pkin(12));
t86 = sin(qJ(4));
t57 = -t116 * t82 + t80 * t86;
t58 = t116 * t80 + t82 * t86;
t36 = t57 * mrSges(5,1) + mrSges(5,2) * t58;
t61 = -t82 * mrSges(4,1) + mrSges(4,2) * t80;
t124 = -t36 - t61;
t89 = cos(qJ(5));
t111 = t58 * t89;
t123 = Ifges(6,5) * t111 + Ifges(6,3) * t57;
t81 = sin(pkin(6));
t87 = sin(qJ(2));
t110 = t81 * t87;
t83 = cos(pkin(6));
t48 = -t110 * t80 + t82 * t83;
t49 = t110 * t82 + t80 * t83;
t24 = -t116 * t48 + t49 * t86;
t21 = t24 ^ 2;
t108 = pkin(8) + qJ(3);
t100 = t108 * t80;
t62 = t108 * t82;
t40 = t100 * t116 + t62 * t86;
t122 = t40 ^ 2;
t77 = t82 ^ 2;
t121 = 0.2e1 * t40;
t120 = m(7) * pkin(5);
t118 = -m(4) - m(5);
t117 = -pkin(10) - pkin(9);
t85 = sin(qJ(5));
t115 = Ifges(6,4) * t85;
t114 = Ifges(6,4) * t89;
t113 = t24 * t40;
t112 = t58 * t85;
t90 = cos(qJ(2));
t109 = t81 * t90;
t70 = -pkin(3) * t82 - pkin(2);
t32 = pkin(4) * t57 - pkin(9) * t58 + t70;
t42 = -t100 * t86 + t116 * t62;
t13 = t32 * t85 + t42 * t89;
t84 = sin(qJ(6));
t88 = cos(qJ(6));
t59 = -t84 * t85 + t88 * t89;
t60 = t84 * t89 + t85 * t88;
t107 = Ifges(7,5) * t60 + Ifges(7,6) * t59;
t63 = -t89 * mrSges(6,1) + t85 * mrSges(6,2);
t106 = t63 - mrSges(5,1);
t105 = Ifges(6,5) * t85 + Ifges(6,6) * t89;
t104 = t80 ^ 2 + t77;
t103 = t85 ^ 2 + t89 ^ 2;
t27 = t60 * t58;
t28 = t59 * t58;
t102 = Ifges(7,5) * t28 - Ifges(7,6) * t27 + Ifges(7,3) * t57;
t26 = t116 * t49 + t48 * t86;
t16 = -t109 * t89 - t26 * t85;
t17 = -t109 * t85 + t26 * t89;
t5 = t16 * t88 - t17 * t84;
t6 = t16 * t84 + t17 * t88;
t101 = mrSges(7,1) * t5 - t6 * mrSges(7,2);
t37 = -mrSges(7,1) * t59 + t60 * mrSges(7,2);
t12 = t32 * t89 - t42 * t85;
t99 = mrSges(6,1) * t85 + mrSges(6,2) * t89;
t98 = -t16 * t85 + t17 * t89;
t97 = -t48 * t80 + t49 * t82;
t66 = t117 * t85;
t67 = t117 * t89;
t44 = t66 * t88 + t67 * t84;
t45 = t66 * t84 - t67 * t88;
t96 = mrSges(7,1) * t44 - t45 * mrSges(7,2) + t107;
t10 = -pkin(10) * t112 + t13;
t7 = pkin(5) * t57 - pkin(10) * t111 + t12;
t2 = -t10 * t84 + t7 * t88;
t3 = t10 * t88 + t7 * t84;
t95 = mrSges(7,1) * t2 - t3 * mrSges(7,2) + t102;
t94 = (mrSges(7,1) * t88 - mrSges(7,2) * t84) * pkin(5);
t76 = t81 ^ 2;
t71 = -pkin(5) * t89 - pkin(4);
t68 = t76 * t90 ^ 2;
t65 = Ifges(6,1) * t85 + t114;
t64 = Ifges(6,2) * t89 + t115;
t39 = Ifges(7,1) * t60 + Ifges(7,4) * t59;
t38 = Ifges(7,4) * t60 + Ifges(7,2) * t59;
t34 = mrSges(6,1) * t57 - mrSges(6,3) * t111;
t33 = -mrSges(6,2) * t57 - mrSges(6,3) * t112;
t31 = t99 * t58;
t20 = pkin(5) * t112 + t40;
t19 = Ifges(6,5) * t57 + (Ifges(6,1) * t89 - t115) * t58;
t18 = Ifges(6,6) * t57 + (-Ifges(6,2) * t85 + t114) * t58;
t15 = mrSges(7,1) * t57 - mrSges(7,3) * t28;
t14 = -mrSges(7,2) * t57 - mrSges(7,3) * t27;
t11 = mrSges(7,1) * t27 + mrSges(7,2) * t28;
t9 = Ifges(7,1) * t28 - Ifges(7,4) * t27 + Ifges(7,5) * t57;
t8 = Ifges(7,4) * t28 - Ifges(7,2) * t27 + Ifges(7,6) * t57;
t1 = [m(2) + m(7) * (t5 ^ 2 + t6 ^ 2 + t21) + m(6) * (t16 ^ 2 + t17 ^ 2 + t21) + m(5) * (t26 ^ 2 + t21 + t68) + m(4) * (t48 ^ 2 + t49 ^ 2 + t68) + m(3) * (t76 * t87 ^ 2 + t83 ^ 2 + t68); -t26 * t57 * mrSges(5,3) + t6 * t14 + t5 * t15 + t16 * t34 + t17 * t33 + t97 * mrSges(4,3) + (t58 * mrSges(5,3) + t11 + t31) * t24 + (-t87 * mrSges(3,2) + (mrSges(3,1) + t124) * t90) * t81 + m(7) * (t2 * t5 + t20 * t24 + t3 * t6) + m(6) * (t12 * t16 + t13 * t17 + t113) + m(5) * (-t109 * t70 + t26 * t42 + t113) + m(4) * (pkin(2) * t109 + qJ(3) * t97); Ifges(4,2) * t77 - 0.2e1 * pkin(2) * t61 + 0.2e1 * t20 * t11 + 0.2e1 * t12 * t34 + 0.2e1 * t13 * t33 + 0.2e1 * t3 * t14 + 0.2e1 * t2 * t15 - t27 * t8 + t28 * t9 + t31 * t121 + 0.2e1 * t70 * t36 + Ifges(3,3) + (Ifges(4,1) * t80 + 0.2e1 * Ifges(4,4) * t82) * t80 + 0.2e1 * t104 * qJ(3) * mrSges(4,3) + (mrSges(5,3) * t121 + Ifges(5,1) * t58 - t18 * t85 + t19 * t89) * t58 + (-0.2e1 * t42 * mrSges(5,3) + Ifges(5,2) * t57 + (-Ifges(6,6) * t85 - (2 * Ifges(5,4))) * t58 + t102 + t123) * t57 + m(4) * (qJ(3) ^ 2 * t104 + pkin(2) ^ 2) + m(5) * (t42 ^ 2 + t70 ^ 2 + t122) + m(6) * (t12 ^ 2 + t13 ^ 2 + t122) + m(7) * (t2 ^ 2 + t20 ^ 2 + t3 ^ 2); t118 * t109 + m(7) * (t5 * t59 + t6 * t60) + m(6) * (t16 * t89 + t17 * t85); -m(4) * pkin(2) + t60 * t14 + t59 * t15 + t85 * t33 + t89 * t34 + m(7) * (t2 * t59 + t3 * t60) + m(6) * (t12 * t89 + t13 * t85) + m(5) * t70 - t124; m(6) * t103 + m(7) * (t59 ^ 2 + t60 ^ 2) - t118; -t26 * mrSges(5,2) + (-t5 * t60 + t59 * t6) * mrSges(7,3) + t98 * mrSges(6,3) + (t37 + t106) * t24 + m(7) * (t24 * t71 + t44 * t5 + t45 * t6) + m(6) * (-pkin(4) * t24 + pkin(9) * t98); t71 * t11 + t59 * t8 / 0.2e1 + t60 * t9 / 0.2e1 - Ifges(5,6) * t57 - t27 * t38 / 0.2e1 + t28 * t39 / 0.2e1 - t42 * mrSges(5,2) + t44 * t15 + t45 * t14 - pkin(4) * t31 + t20 * t37 + t106 * t40 + (-t2 * t60 + t3 * t59) * mrSges(7,3) + (t18 / 0.2e1 + pkin(9) * t33 + t13 * mrSges(6,3)) * t89 + (t19 / 0.2e1 - pkin(9) * t34 - t12 * mrSges(6,3)) * t85 + m(6) * (-pkin(4) * t40 + (-t12 * t85 + t13 * t89) * pkin(9)) + m(7) * (t2 * t44 + t20 * t71 + t3 * t45) + (Ifges(5,5) + t89 * t65 / 0.2e1 - t85 * t64 / 0.2e1) * t58 + (t105 + t107) * t57 / 0.2e1; m(7) * (t44 * t59 + t45 * t60); -0.2e1 * pkin(4) * t63 + 0.2e1 * t71 * t37 + t59 * t38 + t60 * t39 + t89 * t64 + t85 * t65 + Ifges(5,3) + m(7) * (t44 ^ 2 + t45 ^ 2 + t71 ^ 2) + m(6) * (pkin(9) ^ 2 * t103 + pkin(4) ^ 2) + 0.2e1 * (-t44 * t60 + t45 * t59) * mrSges(7,3) + 0.2e1 * t103 * pkin(9) * mrSges(6,3); t16 * mrSges(6,1) - t17 * mrSges(6,2) + (t5 * t88 + t6 * t84) * t120 + t101; -Ifges(6,6) * t112 + t12 * mrSges(6,1) - t13 * mrSges(6,2) + (m(7) * (t2 * t88 + t3 * t84) + t84 * t14 + t88 * t15) * pkin(5) + t95 + t123; (t59 * t88 + t60 * t84) * t120 - t63 - t37; -t99 * pkin(9) + (m(7) * (t44 * t88 + t45 * t84) + (t59 * t84 - t60 * t88) * mrSges(7,3)) * pkin(5) + t96 + t105; Ifges(6,3) + Ifges(7,3) + m(7) * (t84 ^ 2 + t88 ^ 2) * pkin(5) ^ 2 + 0.2e1 * t94; t101; t95; -t37; t96; Ifges(7,3) + t94; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
