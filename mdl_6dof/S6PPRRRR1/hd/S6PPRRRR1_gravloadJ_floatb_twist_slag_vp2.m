% Calculate Gravitation load on the joints for
% S6PPRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
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
% Datum: 2019-03-08 19:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PPRRRR1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRR1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRR1_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRRR1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:59:44
% EndTime: 2019-03-08 18:59:45
% DurationCPUTime: 0.65s
% Computational Cost: add. (863->101), mult. (2039->154), div. (0->0), fcn. (2566->16), ass. (0->68)
t120 = mrSges(6,2) - mrSges(7,3);
t52 = sin(qJ(6));
t55 = cos(qJ(6));
t115 = -t55 * mrSges(7,1) + t52 * mrSges(7,2) - mrSges(6,1);
t107 = cos(qJ(3));
t90 = sin(pkin(7));
t91 = sin(pkin(6));
t92 = cos(pkin(13));
t94 = cos(pkin(7));
t95 = cos(pkin(6));
t111 = t92 * t94 * t91 + t95 * t90;
t54 = sin(qJ(3));
t88 = sin(pkin(13));
t72 = t91 * t88;
t35 = t107 * t72 + t111 * t54;
t74 = t91 * t90;
t41 = -t92 * t74 + t95 * t94;
t53 = sin(qJ(4));
t56 = cos(qJ(4));
t118 = -t35 * t53 + t41 * t56;
t89 = sin(pkin(12));
t71 = t89 * t92;
t93 = cos(pkin(12));
t75 = t93 * t88;
t62 = t95 * t71 + t75;
t73 = t91 * t89;
t112 = t62 * t94 - t90 * t73;
t70 = t89 * t88;
t77 = t93 * t92;
t43 = -t95 * t70 + t77;
t29 = t43 * t107 - t112 * t54;
t37 = t62 * t90 + t94 * t73;
t117 = -t29 * t53 + t37 * t56;
t61 = -t95 * t77 + t70;
t113 = t61 * t94 + t93 * t74;
t42 = t95 * t75 + t71;
t27 = t42 * t107 - t113 * t54;
t76 = t93 * t91;
t36 = t61 * t90 - t94 * t76;
t116 = -t27 * t53 + t36 * t56;
t114 = -m(6) - m(7);
t51 = qJ(4) + qJ(5);
t49 = sin(t51);
t50 = cos(t51);
t110 = m(5) * pkin(3) + t56 * mrSges(5,1) - t53 * mrSges(5,2) + mrSges(4,1) + (m(7) * pkin(5) - t115) * t50 + (m(7) * pkin(11) - t120) * t49;
t109 = -m(5) * pkin(9) - t52 * mrSges(7,1) - t55 * mrSges(7,2) + mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t108 = m(3) + m(4) + m(5) - t114;
t11 = -t27 * t49 + t36 * t50;
t12 = t27 * t50 + t36 * t49;
t86 = t11 * pkin(5) + t12 * pkin(11);
t13 = -t29 * t49 + t37 * t50;
t14 = t29 * t50 + t37 * t49;
t85 = t13 * pkin(5) + t14 * pkin(11);
t22 = -t35 * t49 + t41 * t50;
t23 = t35 * t50 + t41 * t49;
t84 = t22 * pkin(5) + t23 * pkin(11);
t83 = t116 * pkin(4);
t82 = t117 * pkin(4);
t81 = t118 * pkin(4);
t69 = t115 * t11 + t120 * t12;
t68 = t115 * t13 + t120 * t14;
t66 = t115 * t22 + t120 * t23;
t57 = -pkin(10) - pkin(9);
t48 = t56 * pkin(4) + pkin(3);
t34 = -t111 * t107 + t54 * t72;
t28 = t112 * t107 + t43 * t54;
t26 = t113 * t107 + t42 * t54;
t1 = [(-m(2) - t108) * g(3) (-t73 * g(1) + t76 * g(2) - t95 * g(3)) * t108 (t114 * (-t34 * t48 - t35 * t57) + t109 * t35 + t110 * t34) * g(3) + (t114 * (-t26 * t48 - t27 * t57) + t109 * t27 + t110 * t26) * g(2) + (t114 * (-t28 * t48 - t29 * t57) + t109 * t29 + t110 * t28) * g(1) (-t118 * mrSges(5,1) - (-t35 * t56 - t41 * t53) * mrSges(5,2) - m(6) * t81 - m(7) * (t81 + t84) + t66) * g(3) + (-t116 * mrSges(5,1) - (-t27 * t56 - t36 * t53) * mrSges(5,2) - m(6) * t83 - m(7) * (t83 + t86) + t69) * g(2) + (-t117 * mrSges(5,1) - (-t29 * t56 - t37 * t53) * mrSges(5,2) - m(6) * t82 - m(7) * (t82 + t85) + t68) * g(1) (-m(7) * t84 + t66) * g(3) + (-m(7) * t86 + t69) * g(2) + (-m(7) * t85 + t68) * g(1), -g(1) * ((-t14 * t52 + t28 * t55) * mrSges(7,1) + (-t14 * t55 - t28 * t52) * mrSges(7,2)) - g(2) * ((-t12 * t52 + t26 * t55) * mrSges(7,1) + (-t12 * t55 - t26 * t52) * mrSges(7,2)) - g(3) * ((-t23 * t52 + t34 * t55) * mrSges(7,1) + (-t23 * t55 - t34 * t52) * mrSges(7,2))];
taug  = t1(:);
