% Calculate Gravitation load on the joints for
% S6RPRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
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
% Datum: 2019-03-09 07:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRRR9_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR9_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR9_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR9_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR9_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR9_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:23:04
% EndTime: 2019-03-09 07:23:06
% DurationCPUTime: 0.80s
% Computational Cost: add. (394->104), mult. (509->123), div. (0->0), fcn. (465->10), ass. (0->58)
t44 = -pkin(9) - pkin(8);
t96 = mrSges(4,2) - m(5) * pkin(8) + m(6) * t44 + m(7) * (-pkin(10) + t44) - mrSges(5,3) - mrSges(6,3) - mrSges(7,3);
t39 = sin(qJ(3));
t42 = cos(qJ(3));
t95 = t39 * mrSges(4,1) + t96 * t42;
t37 = qJ(4) + qJ(5);
t29 = cos(t37);
t41 = cos(qJ(4));
t33 = t41 * pkin(4);
t23 = pkin(5) * t29 + t33;
t94 = -m(5) * pkin(3) - m(6) * (t33 + pkin(3)) - m(7) * (pkin(3) + t23);
t40 = sin(qJ(1));
t43 = cos(qJ(1));
t92 = g(1) * t40 - g(2) * t43;
t91 = m(6) * pkin(4);
t90 = -mrSges(2,1) + mrSges(3,2) - mrSges(4,3);
t89 = -m(4) - m(5) - m(6) - m(7);
t88 = m(3) - t89;
t28 = sin(t37);
t72 = pkin(5) * t28;
t38 = sin(qJ(4));
t73 = pkin(4) * t38;
t22 = t72 + t73;
t87 = m(6) * t73 + m(7) * t22;
t86 = t94 * t39 + mrSges(2,2) - mrSges(3,3) - t95;
t84 = -mrSges(5,1) - t91;
t32 = qJ(6) + t37;
t25 = sin(t32);
t26 = cos(t32);
t51 = -mrSges(7,1) * t25 - mrSges(7,2) * t26;
t83 = mrSges(6,1) * t28 + mrSges(6,2) * t29 - t51;
t65 = t39 * t43;
t15 = t28 * t65 + t29 * t40;
t16 = -t28 * t40 + t29 * t65;
t7 = t25 * t65 + t26 * t40;
t8 = -t25 * t40 + t26 * t65;
t74 = t7 * mrSges(7,1) + t8 * mrSges(7,2);
t82 = -t15 * mrSges(6,1) - t16 * mrSges(6,2) - t74;
t66 = t39 * t40;
t13 = -t28 * t66 + t29 * t43;
t14 = t28 * t43 + t29 * t66;
t5 = -t25 * t66 + t26 * t43;
t6 = t25 * t43 + t26 * t66;
t75 = t5 * mrSges(7,1) - t6 * mrSges(7,2);
t81 = -t13 * mrSges(6,1) + t14 * mrSges(6,2) - t75;
t80 = mrSges(5,1) * t41 + mrSges(6,1) * t29 + mrSges(7,1) * t26 - mrSges(5,2) * t38 - mrSges(6,2) * t28 - mrSges(7,2) * t25 - t94;
t77 = m(7) * pkin(5);
t69 = g(3) * t42;
t67 = t38 * t43;
t64 = t40 * t41;
t63 = t41 * t43;
t62 = t43 * t22;
t60 = t43 * pkin(1) + t40 * qJ(2);
t19 = t38 * t65 + t64;
t17 = -t38 * t66 + t63;
t20 = -t38 * t40 + t39 * t63;
t18 = t39 * t64 + t67;
t1 = [(-t67 * t91 - m(3) * t60 - m(7) * t62 - t18 * mrSges(5,1) - t14 * mrSges(6,1) - t6 * mrSges(7,1) - t17 * mrSges(5,2) - t13 * mrSges(6,2) - t5 * mrSges(7,2) + t89 * (t43 * pkin(7) + t60) + t90 * t43 + t86 * t40) * g(2) + (-t20 * mrSges(5,1) - t16 * mrSges(6,1) - t8 * mrSges(7,1) + t19 * mrSges(5,2) + t15 * mrSges(6,2) + t7 * mrSges(7,2) + (m(3) * pkin(1) + t89 * (-pkin(1) - pkin(7)) + t87 - t90) * t40 + (-t88 * qJ(2) + t86) * t43) * g(1), -t92 * t88 (t39 * t80 + t95) * g(3) + t92 * ((-mrSges(4,1) - t80) * t42 + t96 * t39) (mrSges(5,1) * t38 + mrSges(5,2) * t41 + t83 + t87) * t69 + (-t20 * mrSges(5,2) - m(7) * (t23 * t40 + t39 * t62) + t84 * t19 + t82) * g(2) + (t18 * mrSges(5,2) - m(7) * (-t22 * t66 + t23 * t43) + t84 * t17 + t81) * g(1) (m(7) * t72 + t83) * t69 + (-t15 * t77 + t82) * g(2) + (-t13 * t77 + t81) * g(1), -g(1) * t75 - g(2) * t74 - t51 * t69];
taug  = t1(:);
