% Calculate Gravitation load on the joints for
% S6PRRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
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
% Datum: 2018-11-23 15:12
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6PRRPRP2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP2_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRP2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:11:43
% EndTime: 2018-11-23 15:11:44
% DurationCPUTime: 0.98s
% Computational Cost: add. (1244->98), mult. (1316->134), div. (0->0), fcn. (1279->16), ass. (0->55)
t65 = sin(qJ(5));
t68 = cos(qJ(5));
t86 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t89 = m(7) * pkin(5) + mrSges(6,1) + mrSges(7,1);
t134 = -t65 * t86 + t68 * t89 + mrSges(5,1);
t121 = -m(6) - m(7);
t125 = mrSges(5,2) - mrSges(7,2) - mrSges(6,3);
t60 = qJ(3) + pkin(11);
t58 = sin(t60);
t59 = cos(t60);
t66 = sin(qJ(3));
t69 = cos(qJ(3));
t131 = t121 * (pkin(4) * t59 + pkin(9) * t58) - m(4) * pkin(2) - t69 * mrSges(4,1) + t66 * mrSges(4,2) + t125 * t58 - mrSges(3,1) - t134 * t59;
t120 = m(5) - t121;
t61 = sin(pkin(10));
t62 = sin(pkin(6));
t108 = t61 * t62;
t101 = cos(pkin(10));
t99 = pkin(6) + qJ(2);
t83 = sin(t99) / 0.2e1;
t100 = pkin(6) - qJ(2);
t90 = sin(t100);
t49 = t83 - t90 / 0.2e1;
t70 = cos(qJ(2));
t79 = t101 * t70 - t61 * t49;
t129 = t69 * t108 - t79 * t66;
t84 = cos(t99) / 0.2e1;
t91 = cos(t100);
t50 = t84 - t91 / 0.2e1;
t63 = cos(pkin(6));
t128 = t50 * t66 + t63 * t69;
t124 = m(4) * pkin(8) + t65 * t89 + t68 * t86 - mrSges(3,2) + mrSges(4,3) + mrSges(5,3);
t80 = t101 * t49 + t61 * t70;
t95 = t62 * t101;
t74 = -t80 * t66 - t69 * t95;
t73 = t74 * pkin(3);
t88 = t129 * pkin(3);
t87 = t128 * pkin(3);
t72 = t91 / 0.2e1 + t84;
t67 = sin(qJ(2));
t64 = -qJ(4) - pkin(8);
t57 = pkin(3) * t69 + pkin(2);
t48 = t83 + t90 / 0.2e1;
t37 = t101 * t67 + t61 * t72;
t34 = -t101 * t72 + t61 * t67;
t27 = -t50 * t59 + t58 * t63;
t26 = t50 * t58 + t59 * t63;
t16 = t58 * t108 + t59 * t79;
t15 = t59 * t108 - t58 * t79;
t14 = -t58 * t95 + t59 * t80;
t13 = -t58 * t80 - t59 * t95;
t9 = t27 * t65 + t48 * t68;
t3 = t16 * t65 - t37 * t68;
t1 = t14 * t65 - t34 * t68;
t2 = [(-m(2) - m(3) - m(4) - t120) * g(3) (-t120 * (t48 * t57 + t50 * t64) + t124 * t50 + t131 * t48) * g(3) + (-t120 * (-t34 * t57 - t64 * t80) - t124 * t80 - t131 * t34) * g(2) + (-t120 * (-t37 * t57 - t64 * t79) - t124 * t79 - t131 * t37) * g(1) (-t128 * mrSges(4,1) - (t50 * t69 - t63 * t66) * mrSges(4,2) - m(5) * t87 + t121 * (t26 * pkin(4) + pkin(9) * t27 + t87) + t125 * t27 - t134 * t26) * g(3) + (-t74 * mrSges(4,1) - (t66 * t95 - t69 * t80) * mrSges(4,2) - m(5) * t73 + t125 * t14 - t134 * t13 + t121 * (t13 * pkin(4) + t14 * pkin(9) + t73)) * g(2) + (-t129 * mrSges(4,1) - (-t66 * t108 - t69 * t79) * mrSges(4,2) - m(5) * t88 + t121 * (t15 * pkin(4) + pkin(9) * t16 + t88) + t125 * t16 - t134 * t15) * g(1), t120 * (-g(1) * t37 - g(2) * t34 + g(3) * t48) (t89 * t9 + t86 * (t27 * t68 - t48 * t65)) * g(3) + (t86 * (t14 * t68 + t34 * t65) + t89 * t1) * g(2) + (t86 * (t16 * t68 + t37 * t65) + t89 * t3) * g(1) (-g(1) * t3 - g(2) * t1 - g(3) * t9) * m(7)];
taug  = t2(:);
