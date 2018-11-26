% Calculate Gravitation load on the joints for
% S6RRPRPR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6]';
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
% Datum: 2018-11-23 17:10
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRPRPR14_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR14_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR14_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR14_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR14_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR14_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:09:48
% EndTime: 2018-11-23 17:09:49
% DurationCPUTime: 1.06s
% Computational Cost: add. (1176->119), mult. (1419->148), div. (0->0), fcn. (1368->14), ass. (0->71)
t49 = sin(qJ(6));
t53 = cos(qJ(6));
t111 = m(6) + m(7);
t69 = -qJ(5) * t111 + mrSges(5,2) - mrSges(6,3);
t114 = -t49 * mrSges(7,1) - t53 * mrSges(7,2) + t69;
t92 = m(5) + t111;
t88 = m(4) + t92;
t107 = t88 * qJ(3) - mrSges(3,2) + mrSges(4,3);
t104 = m(7) * pkin(10) + mrSges(5,1) - mrSges(6,2) + mrSges(7,3);
t112 = pkin(4) * t111 + t104;
t106 = mrSges(4,2) - mrSges(3,1) - mrSges(6,1) - mrSges(5,3);
t73 = t53 * mrSges(7,1) - t49 * mrSges(7,2);
t105 = m(7) * (pkin(5) + pkin(9)) + t73 - t106;
t108 = -m(7) * pkin(5) - t92 * pkin(9) + t106;
t50 = sin(qJ(4));
t54 = cos(qJ(4));
t102 = -t104 * t50 - t114 * t54 - t107;
t99 = pkin(4) * t50;
t47 = sin(pkin(6));
t52 = sin(qJ(1));
t98 = t47 * t52;
t56 = cos(qJ(1));
t97 = t47 * t56;
t55 = cos(qJ(2));
t96 = t52 * t55;
t95 = t56 * t55;
t93 = t56 * pkin(1) + pkin(8) * t98;
t91 = pkin(6) - qJ(2);
t90 = pkin(6) + qJ(2);
t80 = sin(t90);
t75 = t80 / 0.2e1;
t81 = sin(t91);
t67 = t75 - t81 / 0.2e1;
t29 = -t52 * t67 + t95;
t89 = t29 * pkin(2) + t93;
t87 = -pkin(1) * t52 + pkin(8) * t97;
t51 = sin(qJ(2));
t77 = cos(t91) / 0.2e1;
t82 = cos(t90);
t59 = t77 + t82 / 0.2e1;
t25 = t51 * t52 - t56 * t59;
t19 = t25 * pkin(2);
t86 = -pkin(9) * t25 - t19;
t28 = t56 * t51 + t52 * t59;
t21 = t28 * pkin(2);
t85 = -pkin(9) * t28 - t21;
t76 = t81 / 0.2e1;
t37 = t75 + t76;
t36 = t37 * pkin(2);
t84 = pkin(9) * t37 + t36;
t83 = pkin(3) * t98 + t89;
t26 = t56 * t67 + t96;
t79 = -t26 * pkin(2) + t87;
t74 = mrSges(2,2) + (-mrSges(4,1) - mrSges(3,3)) * t47;
t70 = pkin(3) * t97 + t79;
t68 = t76 - t80 / 0.2e1;
t66 = -t25 * t50 + t54 * t97;
t9 = t25 * t54 + t50 * t97;
t48 = cos(pkin(6));
t38 = t77 - t82 / 0.2e1;
t31 = t38 * t99;
t30 = t52 * t68 + t95;
t27 = -t56 * t68 + t96;
t23 = t37 * t54 + t48 * t50;
t14 = t30 * t99;
t13 = t27 * t99;
t8 = t28 * t50 + t54 * t98;
t7 = -t28 * t54 + t50 * t98;
t2 = t29 * t53 + t49 * t7;
t1 = -t29 * t49 + t53 * t7;
t3 = [(-t56 * mrSges(2,1) - m(3) * t93 - m(4) * t89 - m(5) * t83 - t2 * mrSges(7,1) - t1 * mrSges(7,2) + t69 * t7 - t107 * t28 - t104 * t8 + t74 * t52 + t108 * t29 - t111 * (t8 * pkin(4) + t83)) * g(2) + (t52 * mrSges(2,1) - m(3) * t87 - m(4) * t79 - m(5) * t70 + t107 * t25 + t74 * t56 - t104 * t66 + t114 * t9 + (-t108 + t73) * t26 + t111 * (-pkin(4) * t66 - t70)) * g(1) (-m(4) * t36 - m(5) * t84 - m(6) * (t31 + t84) - m(7) * (t31 + t36) + t102 * t38 - t105 * t37) * g(3) + (m(4) * t19 - m(5) * t86 - m(6) * (t13 + t86) - m(7) * (t13 - t19) + t102 * t27 + t105 * t25) * g(2) + (m(4) * t21 - m(5) * t85 - m(6) * (t14 + t85) - m(7) * (t14 - t21) + t102 * t30 + t105 * t28) * g(1) (-g(1) * t28 - g(2) * t25 + g(3) * t37) * t88 (t114 * (-t37 * t50 + t48 * t54) + t112 * t23) * g(3) + (-t112 * t9 - t114 * t66) * g(2) + (t112 * t7 + t114 * t8) * g(1), t111 * (-g(1) * t7 + g(2) * t9 - g(3) * t23) -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * ((-t26 * t49 - t53 * t9) * mrSges(7,1) + (-t26 * t53 + t49 * t9) * mrSges(7,2)) - g(3) * ((t23 * t53 - t38 * t49) * mrSges(7,1) + (-t23 * t49 - t38 * t53) * mrSges(7,2))];
taug  = t3(:);
