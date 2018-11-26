% Calculate Gravitation load on the joints for
% S6RPRRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
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
% Datum: 2018-11-23 16:28
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RPRRRP8_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP8_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP8_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP8_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP8_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP8_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:28:24
% EndTime: 2018-11-23 16:28:25
% DurationCPUTime: 0.71s
% Computational Cost: add. (367->113), mult. (501->135), div. (0->0), fcn. (447->8), ass. (0->61)
t92 = mrSges(6,2) - mrSges(7,3);
t88 = mrSges(6,3) + mrSges(7,2);
t91 = -m(3) - m(4);
t90 = -m(6) - m(7);
t33 = qJ(3) + qJ(4);
t29 = cos(t33);
t34 = sin(qJ(5));
t39 = cos(qJ(1));
t69 = t34 * t39;
t28 = sin(t33);
t73 = t28 * t39;
t89 = -t29 * mrSges(6,2) * t69 - mrSges(5,2) * t73;
t36 = sin(qJ(1));
t77 = g(2) * t39;
t87 = -g(1) * t36 + t77;
t37 = cos(qJ(5));
t67 = t37 * mrSges(6,1);
t59 = t29 * t67;
t68 = t36 * t37;
t60 = t29 * t68;
t70 = t34 * t36;
t61 = t29 * t70;
t72 = t29 * t36;
t74 = t28 * t36;
t86 = -mrSges(5,1) * t72 - mrSges(7,1) * t60 - t36 * t59 + t92 * t61 - t88 * t74;
t43 = m(7) * (-pkin(5) * t37 - qJ(6) * t34 - pkin(4)) - t37 * mrSges(7,1) - t34 * mrSges(7,3);
t51 = t28 * mrSges(5,1) + t29 * mrSges(5,2);
t85 = t51 - t88 * t29 + (-mrSges(6,2) * t34 - t43 + t67) * t28;
t84 = mrSges(2,1) + mrSges(4,3) + mrSges(5,3) - mrSges(3,2);
t35 = sin(qJ(3));
t38 = cos(qJ(3));
t52 = t35 * mrSges(4,1) + t38 * mrSges(4,2);
t83 = mrSges(2,2) - mrSges(3,3) - t52 - t51;
t82 = m(7) * pkin(5) + mrSges(6,1) + mrSges(7,1);
t81 = m(7) * qJ(6) - t92;
t80 = pkin(3) * t38;
t79 = pkin(9) * t28;
t76 = g(3) * t29;
t75 = t35 * pkin(3);
t71 = t29 * t39;
t66 = t39 * t37;
t65 = pkin(4) * t72 + pkin(9) * t74;
t64 = t39 * pkin(1) + t36 * qJ(2);
t63 = m(5) * t80;
t31 = t39 * qJ(2);
t58 = -pkin(1) * t36 + t31;
t25 = t29 * pkin(9);
t57 = -t28 * pkin(4) + t25;
t56 = pkin(5) * t60 + qJ(6) * t61 + t65;
t54 = -pkin(4) * t29 - t79;
t40 = -pkin(8) - pkin(7);
t50 = t36 * t40 + t39 * t75 + t58;
t48 = t36 * t75 - t39 * t40 + t64;
t47 = -t28 * mrSges(6,3) - t59;
t42 = t43 * t29;
t21 = t36 * t80;
t4 = t28 * t66 - t70;
t3 = t28 * t69 + t68;
t2 = t28 * t68 + t69;
t1 = t28 * t70 - t66;
t5 = [(-m(5) * t48 + t88 * t72 + t91 * t64 + t90 * (pkin(4) * t74 - pkin(9) * t72 + t48) - t82 * t2 - t81 * t1 + (-m(4) * pkin(7) - t84) * t39 + t83 * t36) * g(2) + (-m(3) * t58 - m(4) * t31 - m(5) * t50 + t88 * t71 + t90 * (pkin(4) * t73 - pkin(9) * t71 + t50) - t82 * t4 - t81 * t3 + t83 * t39 + (-m(4) * (-pkin(1) - pkin(7)) + t84) * t36) * g(1), t87 * (m(5) - t90 - t91) -(m(7) * (-t79 - t80) - t28 * mrSges(7,2) + t42) * t77 + t87 * (mrSges(4,1) * t38 - mrSges(4,2) * t35) + ((mrSges(5,1) * t29 + t63 - m(6) * (t54 - t80) - t47) * t39 + t89) * g(2) + (-(-mrSges(5,2) * t28 + t63) * t36 - m(6) * (t21 + t65) - m(7) * (t21 + t56) + t86) * g(1) + (t52 + m(5) * t75 - m(6) * (t57 - t75) - m(7) * (t25 - t75) + t85) * g(3) -((-m(7) * pkin(9) - mrSges(7,2)) * t28 + t42) * t77 + (mrSges(5,1) * t71 - (m(6) * t54 + t47) * t39 + t89) * g(2) + (-m(6) * t57 - m(7) * t25 + t85) * g(3) + (-m(6) * t65 - m(7) * t56 + mrSges(5,2) * t74 + t86) * g(1) (t82 * t34 - t81 * t37) * t76 + (-t82 * t3 + t81 * t4) * g(2) + (t82 * t1 - t81 * t2) * g(1) (-g(1) * t1 + g(2) * t3 - t34 * t76) * m(7)];
taug  = t5(:);
