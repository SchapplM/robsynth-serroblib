% Calculate Gravitation load on the joints for
% S5PRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d2,d3,d4,d5,theta1]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:27
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRRR10_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(11,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR10_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR10_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR10_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR10_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR10_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:23:16
% EndTime: 2019-12-05 17:23:19
% DurationCPUTime: 0.92s
% Computational Cost: add. (663->107), mult. (1835->178), div. (0->0), fcn. (2307->14), ass. (0->62)
t105 = m(5) + m(6);
t52 = sin(qJ(5));
t56 = cos(qJ(5));
t107 = m(6) * pkin(4) + t56 * mrSges(6,1) - t52 * mrSges(6,2) + mrSges(5,1);
t99 = -m(6) * pkin(10) + mrSges(5,2) - mrSges(6,3);
t96 = -t52 * mrSges(6,1) - t56 * mrSges(6,2) + mrSges(4,2) - mrSges(5,3);
t53 = sin(qJ(4));
t57 = cos(qJ(4));
t106 = pkin(3) * t105 + t107 * t57 - t99 * t53 + mrSges(4,1);
t50 = sin(pkin(6));
t55 = sin(qJ(2));
t84 = cos(pkin(11));
t86 = cos(pkin(5));
t68 = t86 * t84;
t83 = sin(pkin(11));
t94 = cos(qJ(2));
t61 = t83 * t55 - t68 * t94;
t51 = sin(pkin(5));
t78 = t51 * t84;
t85 = cos(pkin(6));
t103 = t50 * t78 + t61 * t85;
t67 = t86 * t83;
t62 = t84 * t55 + t67 * t94;
t77 = t51 * t83;
t102 = -t50 * t77 + t62 * t85;
t97 = -t105 * pkin(9) + t96;
t93 = cos(qJ(3));
t40 = t55 * t68 + t83 * t94;
t92 = t40 * t50;
t41 = -t55 * t67 + t84 * t94;
t91 = t41 * t50;
t90 = t50 * t53;
t89 = t50 * t57;
t88 = t51 * t55;
t80 = t51 * t94;
t82 = t50 * t88;
t87 = pkin(2) * t80 + pkin(8) * t82;
t79 = t50 * t86;
t54 = sin(qJ(3));
t76 = t54 * t85;
t75 = -t61 * pkin(2) + pkin(8) * t92;
t74 = -t62 * pkin(2) + pkin(8) * t91;
t71 = t85 * t93;
t39 = -t50 * t80 + t85 * t86;
t36 = (-t55 * t76 + t93 * t94) * t51;
t35 = (t54 * t94 + t55 * t71) * t51;
t29 = t50 * t62 + t77 * t85;
t28 = t50 * t61 - t78 * t85;
t27 = t54 * t79 + (t93 * t55 + t76 * t94) * t51;
t26 = t54 * t88 - t71 * t80 - t79 * t93;
t22 = -t41 * t76 - t62 * t93;
t21 = t41 * t71 - t54 * t62;
t20 = -t40 * t76 - t61 * t93;
t19 = t40 * t71 - t54 * t61;
t16 = t27 * t57 + t39 * t53;
t14 = -t102 * t54 + t41 * t93;
t13 = t102 * t93 + t41 * t54;
t12 = -t103 * t54 + t40 * t93;
t11 = t103 * t93 + t40 * t54;
t4 = t14 * t57 + t29 * t53;
t2 = t12 * t57 + t28 * t53;
t1 = [(-m(2) - m(3) - m(4) - t105) * g(3), (-(mrSges(3,1) * t94 - mrSges(3,2) * t55) * t51 - m(4) * t87 - t36 * mrSges(4,1) - mrSges(4,3) * t82 - t105 * (t36 * pkin(3) + pkin(9) * t35 + t87) - t107 * (t36 * t57 + t53 * t82) + t96 * t35 + t99 * (t36 * t53 - t57 * t82)) * g(3) + (-m(4) * t75 + t61 * mrSges(3,1) - t20 * mrSges(4,1) + t40 * mrSges(3,2) - mrSges(4,3) * t92 - t105 * (t20 * pkin(3) + pkin(9) * t19 + t75) + t99 * (t20 * t53 - t40 * t89) - t107 * (t20 * t57 + t40 * t90) + t96 * t19) * g(2) + (-m(4) * t74 + t62 * mrSges(3,1) - t22 * mrSges(4,1) + t41 * mrSges(3,2) - mrSges(4,3) * t91 - t105 * (t22 * pkin(3) + pkin(9) * t21 + t74) + t99 * (t22 * t53 - t41 * t89) - t107 * (t22 * t57 + t41 * t90) + t96 * t21) * g(1), (t106 * t26 + t97 * t27) * g(3) + (t106 * t11 + t97 * t12) * g(2) + (t106 * t13 + t97 * t14) * g(1), (t99 * t16 - t107 * (-t27 * t53 + t39 * t57)) * g(3) + (t99 * t2 - t107 * (-t12 * t53 + t28 * t57)) * g(2) + (t99 * t4 - t107 * (-t14 * t53 + t29 * t57)) * g(1), -g(1) * ((t13 * t56 - t4 * t52) * mrSges(6,1) + (-t13 * t52 - t4 * t56) * mrSges(6,2)) - g(2) * ((t11 * t56 - t2 * t52) * mrSges(6,1) + (-t11 * t52 - t2 * t56) * mrSges(6,2)) - g(3) * ((-t16 * t52 + t26 * t56) * mrSges(6,1) + (-t16 * t56 - t26 * t52) * mrSges(6,2))];
taug = t1(:);
