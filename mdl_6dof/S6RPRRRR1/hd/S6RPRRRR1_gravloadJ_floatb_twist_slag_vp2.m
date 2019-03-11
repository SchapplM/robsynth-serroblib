% Calculate Gravitation load on the joints for
% S6RPRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 06:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRRR1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR1_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:54:10
% EndTime: 2019-03-09 06:54:11
% DurationCPUTime: 0.54s
% Computational Cost: add. (519->102), mult. (381->108), div. (0->0), fcn. (314->12), ass. (0->60)
t33 = qJ(3) + qJ(4);
t28 = qJ(5) + t33;
t21 = sin(t28);
t22 = cos(t28);
t34 = sin(qJ(6));
t68 = t34 * mrSges(7,2);
t99 = t21 * t68 + t22 * (m(7) * pkin(10) + mrSges(7,3));
t91 = t22 * pkin(5) + t21 * pkin(10);
t98 = m(7) * t91;
t97 = -t22 * mrSges(6,1) + (mrSges(6,2) - mrSges(7,3)) * t21;
t26 = sin(t33);
t27 = cos(t33);
t77 = mrSges(6,2) * t22;
t96 = mrSges(5,1) * t26 + mrSges(6,1) * t21 + mrSges(5,2) * t27 + t77;
t31 = qJ(1) + pkin(11);
t24 = sin(t31);
t25 = cos(t31);
t95 = g(1) * t25 + g(2) * t24;
t92 = -m(6) - m(7);
t61 = t27 * mrSges(5,1) - t26 * mrSges(5,2);
t37 = cos(qJ(6));
t67 = t37 * mrSges(7,1);
t90 = -(t67 - t68) * t22 + t97;
t88 = m(3) + m(4) + m(5);
t86 = -t61 + t90;
t38 = cos(qJ(3));
t29 = t38 * pkin(3);
t35 = sin(qJ(3));
t54 = t38 * mrSges(4,1) - t35 * mrSges(4,2);
t85 = mrSges(3,1) + m(5) * (t29 + pkin(2)) + t61 + m(4) * pkin(2) + t54 - t97;
t40 = -pkin(8) - pkin(7);
t84 = -m(4) * pkin(7) + m(5) * t40 + mrSges(3,2) - mrSges(4,3) - mrSges(5,3) - mrSges(6,3);
t83 = pkin(4) * t26;
t20 = pkin(4) * t27;
t82 = pkin(5) * t21;
t79 = t35 * pkin(3);
t36 = sin(qJ(1));
t78 = t36 * pkin(1);
t39 = cos(qJ(1));
t30 = t39 * pkin(1);
t73 = t24 * t34;
t72 = t24 * t37;
t71 = t25 * t34;
t70 = t25 * t37;
t65 = t20 + t29;
t63 = t21 * t67;
t62 = t20 + t91;
t59 = t99 * t24;
t58 = t99 * t25;
t13 = -t79 - t83;
t43 = m(7) * (t13 - t82) - t63;
t42 = m(7) * (-t82 - t83) - t63;
t41 = t77 + (m(7) * pkin(5) + mrSges(6,1) + t67) * t21;
t32 = -pkin(9) + t40;
t12 = pkin(2) + t65;
t4 = t22 * t70 + t73;
t3 = -t22 * t71 + t72;
t2 = -t22 * t72 + t71;
t1 = t22 * t73 + t70;
t5 = [(-t39 * mrSges(2,1) - t4 * mrSges(7,1) + t36 * mrSges(2,2) - t3 * mrSges(7,2) + t92 * (t25 * t12 - t24 * t32 + t30) - t88 * t30 + t84 * t24 + (-t85 - t98) * t25) * g(2) + (t36 * mrSges(2,1) - t2 * mrSges(7,1) + t39 * mrSges(2,2) - t1 * mrSges(7,2) + t92 * (-t25 * t32 - t78) + t88 * t78 + t84 * t25 + (m(6) * t12 - m(7) * (-t12 - t91) + t85) * t24) * g(1) (-t88 + t92) * g(3), -g(1) * (t43 * t25 + t58) - g(2) * (t43 * t24 + t59) + (-t54 - m(5) * t29 - m(6) * t65 - m(7) * (t29 + t62) + t86) * g(3) + t95 * (m(5) * t79 - m(6) * t13 + mrSges(4,1) * t35 + mrSges(4,2) * t38 + t96) -g(1) * (t42 * t25 + t58) - g(2) * (t42 * t24 + t59) + (-m(6) * t20 - m(7) * t62 + t86) * g(3) + (m(6) * t83 + t96) * t95 (t90 - t98) * g(3) + (t41 * t24 - t59) * g(2) + (t41 * t25 - t58) * g(1), -g(1) * (t3 * mrSges(7,1) - t4 * mrSges(7,2)) - g(2) * (-t1 * mrSges(7,1) + t2 * mrSges(7,2)) - g(3) * (-mrSges(7,1) * t34 - mrSges(7,2) * t37) * t21];
taug  = t5(:);
