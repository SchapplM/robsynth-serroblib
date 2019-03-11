% Calculate Gravitation load on the joints for
% S6RRRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
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
% Datum: 2019-03-09 15:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPPR4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR4_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:31:51
% EndTime: 2019-03-09 15:31:55
% DurationCPUTime: 1.31s
% Computational Cost: add. (551->148), mult. (799->172), div. (0->0), fcn. (805->10), ass. (0->72)
t117 = mrSges(5,1) + mrSges(6,1);
t105 = mrSges(5,2) - mrSges(6,3);
t116 = mrSges(5,3) + mrSges(6,2);
t35 = qJ(3) + pkin(10);
t30 = sin(t35);
t31 = cos(t35);
t115 = pkin(4) * t31 + qJ(5) * t30;
t38 = sin(qJ(3));
t42 = cos(qJ(3));
t37 = sin(qJ(6));
t41 = cos(qJ(6));
t59 = t30 * t37 + t31 * t41;
t60 = t30 * t41 - t31 * t37;
t114 = m(4) * pkin(2) + t42 * mrSges(4,1) + t59 * mrSges(7,1) - t38 * mrSges(4,2) + t60 * mrSges(7,2) - t105 * t30 + t117 * t31;
t44 = cos(qJ(1));
t40 = sin(qJ(1));
t43 = cos(qJ(2));
t83 = t40 * t43;
t10 = t30 * t83 + t31 * t44;
t81 = t44 * t30;
t11 = t31 * t83 - t81;
t108 = t10 * t41 - t11 * t37;
t61 = t10 * t37 + t11 * t41;
t113 = mrSges(7,1) * t108 - t61 * mrSges(7,2);
t109 = m(6) + m(7);
t39 = sin(qJ(2));
t112 = t39 * mrSges(4,3) + mrSges(2,1);
t82 = t43 * t44;
t16 = -t38 * t82 + t40 * t42;
t111 = g(1) * t44 + g(2) * t40;
t101 = m(7) * pkin(5) + t117;
t110 = -m(3) - m(4);
t106 = mrSges(2,2) - mrSges(3,3);
t104 = t115 * t43;
t103 = m(5) + t109;
t65 = t43 * mrSges(3,1) - t39 * mrSges(3,2);
t102 = t116 * t39 + t65;
t98 = pkin(3) * t38;
t96 = pkin(5) * t31;
t93 = g(3) * t39;
t36 = -qJ(4) - pkin(8);
t91 = t36 * t39;
t90 = t38 * t44;
t86 = t39 * t44;
t85 = t40 * t38;
t29 = pkin(3) * t42 + pkin(2);
t22 = t43 * t29;
t80 = t44 * pkin(1) + t40 * pkin(7);
t77 = t36 * t86;
t33 = t44 * pkin(7);
t75 = pkin(3) * t90 + t40 * t91 + t33;
t74 = m(4) * pkin(8) + mrSges(4,3);
t71 = t22 - t91;
t70 = pkin(3) * t85 + t29 * t82 + t80;
t69 = m(7) * (-pkin(9) - t36) - mrSges(7,3);
t12 = -t40 * t31 + t43 * t81;
t13 = t40 * t30 + t31 * t82;
t1 = t12 * t41 - t13 * t37;
t2 = t12 * t37 + t13 * t41;
t68 = mrSges(7,1) * t1 - mrSges(7,2) * t2;
t67 = (mrSges(7,1) * t60 - mrSges(7,2) * t59) * t39;
t66 = pkin(2) * t43 + pkin(8) * t39;
t58 = t16 * pkin(3);
t14 = t38 * t83 + t42 * t44;
t57 = -t29 - t115;
t56 = t69 * t39;
t53 = t13 * pkin(4) + t12 * qJ(5) + t70;
t52 = t14 * pkin(3);
t19 = t39 * t31 * qJ(5);
t17 = t42 * t82 + t85;
t15 = -t42 * t83 + t90;
t3 = [(-t17 * mrSges(4,1) - t16 * mrSges(4,2) - m(5) * (t70 - t77) - m(6) * (t53 - t77) - m(7) * t53 - t2 * mrSges(7,1) - t1 * mrSges(7,2) - t116 * t86 + t110 * t80 + t106 * t40 - t101 * t13 + t105 * t12 + (-m(4) * t66 - t112 - t56 - t65) * t44) * g(2) + (-m(5) * t75 - t15 * mrSges(4,1) + t61 * mrSges(7,1) - t14 * mrSges(4,2) + t108 * mrSges(7,2) - t109 * (-t11 * pkin(4) - qJ(5) * t10 + t75) + t106 * t44 + t110 * t33 + t101 * t11 - t105 * t10 + (m(3) * pkin(1) - m(4) * (-pkin(1) - t66) + (-m(7) * pkin(9) - mrSges(7,3)) * t39 - t103 * (-pkin(1) - t22) + t102 + t112) * t40) * g(1) (-m(5) * t71 - m(6) * (t71 + t104) - m(7) * (t22 + t104) - t56 - t102) * g(3) + ((-m(7) * t96 - t114) * g(3) + t111 * (mrSges(3,2) - t116 - t69 - t74 + (m(5) + m(6)) * t36)) * t43 + (-t74 * g(3) + t111 * (mrSges(3,1) + m(5) * t29 - m(6) * t57 - m(7) * (t57 - t96) + t114)) * t39 (m(5) * t98 + mrSges(4,1) * t38 + mrSges(5,1) * t30 + mrSges(4,2) * t42 + mrSges(5,2) * t31) * t93 + (-m(6) * t19 - (m(6) * (-pkin(4) * t30 - t98) - t30 * mrSges(6,1) + t31 * mrSges(6,3)) * t39 + t67 - (t19 + (-t98 + (-pkin(4) - pkin(5)) * t30) * t39) * m(7)) * g(3) + (m(5) * t52 + t14 * mrSges(4,1) - t15 * mrSges(4,2) - t109 * (-t10 * pkin(4) + t11 * qJ(5) - t52) + t105 * t11 + t101 * t10 + t113) * g(2) + (-m(5) * t58 - t16 * mrSges(4,1) + t17 * mrSges(4,2) + t68 - t109 * (-t12 * pkin(4) + t13 * qJ(5) + t58) + t105 * t13 + t101 * t12) * g(1) (g(3) * t43 - t111 * t39) * t103, t109 * (-g(1) * t12 - g(2) * t10 - t30 * t93) -g(1) * t68 - g(2) * t113 - g(3) * t67];
taug  = t3(:);
