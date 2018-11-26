% Calculate Gravitation load on the joints for
% S6RPRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
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
% Datum: 2018-11-23 16:17
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RPRRPR3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR3_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:16:53
% EndTime: 2018-11-23 16:16:54
% DurationCPUTime: 0.91s
% Computational Cost: add. (505->117), mult. (636->139), div. (0->0), fcn. (650->10), ass. (0->57)
t87 = mrSges(5,1) + mrSges(6,1);
t85 = mrSges(5,2) - mrSges(6,3);
t93 = -mrSges(4,2) + mrSges(6,2) + mrSges(5,3);
t33 = sin(qJ(4));
t37 = cos(qJ(4));
t32 = sin(qJ(6));
t36 = cos(qJ(6));
t48 = t32 * t33 + t36 * t37;
t49 = t32 * t37 - t33 * t36;
t92 = t48 * mrSges(7,1) - t49 * mrSges(7,2) - t85 * t33 + t87 * t37;
t31 = qJ(1) + pkin(10);
t26 = sin(t31);
t27 = cos(t31);
t38 = cos(qJ(3));
t68 = t38 * t33;
t8 = t26 * t68 + t27 * t37;
t67 = t38 * t37;
t9 = t26 * t67 - t27 * t33;
t55 = t32 * t8 + t36 * t9;
t89 = -t32 * t9 + t8 * t36;
t91 = t89 * mrSges(7,1) - t55 * mrSges(7,2);
t90 = m(6) + m(7);
t86 = mrSges(3,2) - mrSges(4,3);
t84 = -m(5) - t90;
t34 = sin(qJ(3));
t83 = t38 * mrSges(4,1) + t93 * t34;
t82 = pkin(8) * t84;
t65 = qJ(5) * t33;
t46 = -pkin(4) * t37 - pkin(3) - t65;
t76 = pkin(5) * t37;
t81 = (m(7) * pkin(9) + mrSges(7,3) - t93) * t38 + (-m(7) * (t46 - t76) - m(6) * t46 + m(5) * pkin(3) + mrSges(4,1) + t92) * t34;
t80 = t34 * mrSges(7,3) - t83;
t79 = m(7) * pkin(5) + t87;
t35 = sin(qJ(1));
t77 = pkin(1) * t35;
t28 = t34 * pkin(8);
t29 = t38 * pkin(3);
t39 = cos(qJ(1));
t30 = t39 * pkin(1);
t74 = t27 * t34;
t73 = t27 * t38;
t72 = t33 * t34;
t66 = t29 + t28;
t64 = t27 * pkin(2) + t26 * pkin(7) + t30;
t63 = -pkin(2) - t29;
t61 = t27 * pkin(7) - t77;
t58 = pkin(4) * t67 + t38 * t65 + t66;
t57 = pkin(3) * t73 + pkin(8) * t74 + t64;
t10 = -t26 * t37 + t27 * t68;
t11 = t26 * t33 + t27 * t67;
t1 = t10 * t36 - t11 * t32;
t2 = t10 * t32 + t11 * t36;
t56 = mrSges(7,1) * t1 - mrSges(7,2) * t2;
t50 = (-mrSges(7,1) * t49 - mrSges(7,2) * t48) * t34;
t43 = t11 * pkin(4) + qJ(5) * t10 + t57;
t20 = t37 * t34 * qJ(5);
t3 = [(-t39 * mrSges(2,1) + t35 * mrSges(2,2) - m(3) * t30 - m(4) * t64 - m(5) * t57 - m(6) * t43 - m(7) * (-pkin(9) * t74 + t43) - t2 * mrSges(7,1) - t1 * mrSges(7,2) + t86 * t26 - t79 * t11 + t85 * t10 + (-mrSges(3,1) + t80) * t27) * g(2) + (m(3) * t77 + t35 * mrSges(2,1) + t55 * mrSges(7,1) + t39 * mrSges(2,2) + t89 * mrSges(7,2) - t90 * (-t9 * pkin(4) - qJ(5) * t8 + t61) + t79 * t9 - t85 * t8 + (-m(4) - m(5)) * t61 + t86 * t27 + (mrSges(3,1) + m(4) * pkin(2) - m(7) * t63 - (m(7) * (-pkin(8) + pkin(9)) + mrSges(7,3)) * t34 + (-m(5) - m(6)) * (t63 - t28) + t83) * t26) * g(1) (-m(3) - m(4) + t84) * g(3) (-m(5) * t66 - m(6) * t58 - m(7) * (-pkin(9) * t34 + t58) + (-m(7) * t76 - t92) * t38 + t80) * g(3) + (t81 * t27 + t73 * t82) * g(1) + (t38 * t82 + t81) * g(2) * t26 (-m(6) * t20 - m(7) * (t20 + (-pkin(4) - pkin(5)) * t72) + t50 + (t85 * t37 + (m(6) * pkin(4) + t87) * t33) * t34) * g(3) + (t85 * t9 - t90 * (-t8 * pkin(4) + qJ(5) * t9) + t79 * t8 + t91) * g(2) + (t56 - t90 * (-t10 * pkin(4) + qJ(5) * t11) + t85 * t11 + t79 * t10) * g(1), t90 * (-g(1) * t10 - g(2) * t8 - g(3) * t72) -g(1) * t56 - g(2) * t91 - g(3) * t50];
taug  = t3(:);
