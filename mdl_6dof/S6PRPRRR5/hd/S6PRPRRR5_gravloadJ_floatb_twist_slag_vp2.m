% Calculate Gravitation load on the joints for
% S6PRPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 15:06
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6PRPRRR5_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR5_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR5_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR5_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR5_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR5_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:06:08
% EndTime: 2018-11-23 15:06:09
% DurationCPUTime: 0.72s
% Computational Cost: add. (967->108), mult. (1023->144), div. (0->0), fcn. (961->16), ass. (0->69)
t109 = -m(6) - m(7);
t110 = m(4) + m(5);
t117 = t109 - t110;
t116 = mrSges(6,2) - mrSges(7,3);
t57 = sin(qJ(6));
t60 = cos(qJ(6));
t111 = -mrSges(7,1) * t60 + mrSges(7,2) * t57 - mrSges(6,1);
t53 = sin(pkin(11));
t55 = cos(pkin(11));
t59 = sin(qJ(2));
t92 = pkin(6) - qJ(2);
t77 = cos(t92) / 0.2e1;
t91 = pkin(6) + qJ(2);
t81 = cos(t91);
t72 = t77 + t81 / 0.2e1;
t31 = t53 * t59 - t55 * t72;
t61 = cos(qJ(4));
t54 = sin(pkin(6));
t58 = sin(qJ(4));
t97 = t54 * t58;
t115 = t31 * t61 + t55 * t97;
t34 = t53 * t72 + t55 * t59;
t114 = t34 * t61 - t53 * t97;
t79 = sin(t91);
t75 = t79 / 0.2e1;
t80 = sin(t92);
t76 = t80 / 0.2e1;
t43 = t75 + t76;
t56 = cos(pkin(6));
t113 = -t43 * t61 - t56 * t58;
t108 = m(5) * pkin(8) + mrSges(7,1) * t57 + mrSges(7,2) * t60 + mrSges(3,1) - mrSges(4,2) + mrSges(5,3) + mrSges(6,3);
t52 = qJ(4) + qJ(5);
t50 = sin(t52);
t51 = cos(t52);
t107 = -t58 * mrSges(5,1) - t61 * mrSges(5,2) + mrSges(3,2) - mrSges(4,3) + (-m(7) * pkin(5) + t111) * t50 + (m(7) * pkin(10) - t116) * t51 + t117 * qJ(3);
t106 = pkin(4) * t58;
t100 = t53 * t54;
t62 = cos(qJ(2));
t99 = t53 * t62;
t98 = t54 * t55;
t96 = t54 * t61;
t95 = t55 * t62;
t93 = t115 * pkin(4);
t11 = -t50 * t100 + t34 * t51;
t12 = t51 * t100 + t34 * t50;
t84 = t11 * pkin(5) + pkin(10) * t12;
t13 = t31 * t51 + t50 * t98;
t14 = -t31 * t50 + t51 * t98;
t83 = t13 * pkin(5) - pkin(10) * t14;
t25 = -t43 * t51 - t50 * t56;
t26 = -t43 * t50 + t51 * t56;
t82 = t25 * pkin(5) + pkin(10) * t26;
t78 = t113 * pkin(4);
t73 = t114 * pkin(4);
t71 = t76 - t79 / 0.2e1;
t70 = t75 - t80 / 0.2e1;
t69 = t111 * t11 + t116 * t12;
t68 = t111 * t13 - t116 * t14;
t67 = t111 * t25 + t116 * t26;
t63 = -pkin(9) - pkin(8);
t44 = t77 - t81 / 0.2e1;
t42 = t43 * pkin(2);
t36 = t53 * t71 + t95;
t35 = -t53 * t70 + t95;
t33 = -t55 * t71 + t99;
t32 = t55 * t70 + t99;
t30 = t34 * pkin(2);
t29 = t31 * pkin(2);
t1 = [(-m(2) - m(3) + t117) * g(3) (t109 * (t44 * t106 - t43 * t63 + t42) - t110 * t42 + t107 * t44 - t108 * t43) * g(3) + (t109 * (t33 * t106 + t31 * t63 - t29) + t110 * t29 + t107 * t33 + t108 * t31) * g(2) + (t109 * (t36 * t106 + t34 * t63 - t30) + t110 * t30 + t107 * t36 + t108 * t34) * g(1) -(-g(1) * t34 - g(2) * t31 + g(3) * t43) * t117 (-t113 * mrSges(5,1) - (t43 * t58 - t56 * t61) * mrSges(5,2) - m(6) * t78 - m(7) * (t78 + t82) + t67) * g(3) + (-t115 * mrSges(5,1) - (-t31 * t58 + t55 * t96) * mrSges(5,2) - m(6) * t93 - m(7) * (t83 + t93) + t68) * g(2) + (-t114 * mrSges(5,1) - (-t34 * t58 - t53 * t96) * mrSges(5,2) - m(6) * t73 - m(7) * (t73 + t84) + t69) * g(1) (-m(7) * t82 + t67) * g(3) + (-m(7) * t83 + t68) * g(2) + (-m(7) * t84 + t69) * g(1), -g(1) * ((-t12 * t57 + t35 * t60) * mrSges(7,1) + (-t12 * t60 - t35 * t57) * mrSges(7,2)) - g(2) * ((t14 * t57 + t32 * t60) * mrSges(7,1) + (t14 * t60 - t32 * t57) * mrSges(7,2)) - g(3) * ((-t26 * t57 + t44 * t60) * mrSges(7,1) + (-t26 * t60 - t44 * t57) * mrSges(7,2))];
taug  = t1(:);
