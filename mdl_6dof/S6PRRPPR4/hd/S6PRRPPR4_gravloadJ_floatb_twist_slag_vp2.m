% Calculate Gravitation load on the joints for
% S6PRRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRPPR4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR4_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPPR4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:12:33
% EndTime: 2019-03-08 21:12:36
% DurationCPUTime: 1.07s
% Computational Cost: add. (547->119), mult. (1423->171), div. (0->0), fcn. (1706->12), ass. (0->74)
t53 = sin(pkin(11));
t55 = cos(pkin(11));
t123 = -pkin(4) * t55 - qJ(5) * t53;
t117 = mrSges(7,3) - mrSges(6,2) - mrSges(5,3);
t122 = m(7) * pkin(9) + t117;
t56 = sin(qJ(6));
t59 = cos(qJ(6));
t121 = t56 * mrSges(7,1) + t59 * mrSges(7,2) - mrSges(5,2) + mrSges(6,3);
t120 = m(7) * pkin(5) + t59 * mrSges(7,1) - t56 * mrSges(7,2) + mrSges(5,1) + mrSges(6,1);
t100 = cos(qJ(3));
t57 = sin(qJ(3));
t118 = -t100 * pkin(3) - qJ(4) * t57;
t115 = m(6) + m(7);
t114 = -mrSges(4,1) * t100 + mrSges(4,2) * t57 - mrSges(3,1);
t113 = mrSges(3,2) - mrSges(4,3);
t58 = sin(qJ(2));
t60 = cos(qJ(2));
t90 = cos(pkin(10));
t91 = cos(pkin(6));
t69 = t91 * t90;
t89 = sin(pkin(10));
t38 = t58 * t69 + t60 * t89;
t54 = sin(pkin(6));
t81 = t54 * t90;
t19 = t100 * t81 + t38 * t57;
t111 = t123 * t19;
t68 = t91 * t89;
t40 = -t58 * t68 + t60 * t90;
t80 = t54 * t89;
t21 = -t100 * t80 + t40 * t57;
t110 = t123 * t21;
t99 = t54 * t58;
t41 = -t100 * t91 + t57 * t99;
t109 = t123 * t41;
t108 = m(5) + t115;
t107 = t120 * t55 + t121 * t53 + mrSges(4,1);
t106 = mrSges(4,2) - m(7) * (-pkin(9) + qJ(4)) + t117;
t104 = -t122 * t57 - t114;
t98 = t54 * t60;
t94 = pkin(2) * t98 + pkin(8) * t99;
t88 = t57 * t98;
t86 = t53 * t100;
t85 = t55 * t100;
t84 = t60 * t100;
t37 = t58 * t89 - t60 * t69;
t83 = -t37 * pkin(2) + pkin(8) * t38;
t39 = t58 * t90 + t60 * t68;
t82 = -t39 * pkin(2) + pkin(8) * t40;
t15 = t19 * pkin(3);
t20 = t100 * t38 - t57 * t81;
t79 = qJ(4) * t20 - t15;
t16 = t21 * pkin(3);
t22 = t100 * t40 + t57 * t80;
t78 = qJ(4) * t22 - t16;
t36 = t41 * pkin(3);
t42 = t100 * t99 + t57 * t91;
t77 = qJ(4) * t42 - t36;
t75 = t54 * t84;
t76 = pkin(3) * t75 + qJ(4) * t88 + t94;
t71 = t118 * t37 + t83;
t70 = t118 * t39 + t82;
t25 = (t53 * t58 + t55 * t84) * t54;
t24 = t53 * t75 - t55 * t99;
t18 = t42 * t55 - t53 * t98;
t17 = t42 * t53 + t55 * t98;
t10 = -t39 * t85 + t40 * t53;
t9 = -t39 * t86 - t40 * t55;
t8 = -t37 * t85 + t38 * t53;
t7 = -t37 * t86 - t38 * t55;
t4 = t22 * t55 + t39 * t53;
t3 = t22 * t53 - t39 * t55;
t2 = t20 * t55 + t37 * t53;
t1 = t20 * t53 - t37 * t55;
t5 = [(-m(2) - m(3) - m(4) - t108) * g(3) (-m(4) * t94 - m(5) * t76 - t115 * (t25 * pkin(4) + t24 * qJ(5) + t76) + (t113 * t58 + t114 * t60) * t54 - t120 * t25 - t121 * t24 + t122 * t88) * g(3) + (-m(4) * t83 - m(5) * t71 - t115 * (t8 * pkin(4) + qJ(5) * t7 + t71) - t120 * t8 - t121 * t7 + t113 * t38 + t104 * t37) * g(2) + (-m(4) * t82 - m(5) * t70 - t115 * (t10 * pkin(4) + qJ(5) * t9 + t70) + t113 * t40 - t121 * t9 - t120 * t10 + t104 * t39) * g(1) (-m(5) * t77 - m(6) * (t77 + t109) - m(7) * (-t36 + t109) + t106 * t42 + t107 * t41) * g(3) + (-m(5) * t79 - m(6) * (t79 + t111) - m(7) * (-t15 + t111) + t106 * t20 + t107 * t19) * g(2) + (-m(5) * t78 - m(6) * (t78 + t110) - m(7) * (-t16 + t110) + t106 * t22 + t107 * t21) * g(1), t108 * (-g(1) * t21 - g(2) * t19 - g(3) * t41) t115 * (-g(1) * t3 - g(2) * t1 - g(3) * t17) -g(1) * ((t3 * t59 - t4 * t56) * mrSges(7,1) + (-t3 * t56 - t4 * t59) * mrSges(7,2)) - g(2) * ((t1 * t59 - t2 * t56) * mrSges(7,1) + (-t1 * t56 - t2 * t59) * mrSges(7,2)) - g(3) * ((t17 * t59 - t18 * t56) * mrSges(7,1) + (-t17 * t56 - t18 * t59) * mrSges(7,2))];
taug  = t5(:);
