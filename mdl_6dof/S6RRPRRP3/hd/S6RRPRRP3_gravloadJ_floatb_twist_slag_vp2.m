% Calculate Gravitation load on the joints for
% S6RRPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-03-09 11:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRRP3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP3_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:47:21
% EndTime: 2019-03-09 11:47:23
% DurationCPUTime: 0.86s
% Computational Cost: add. (486->111), mult. (532->125), div. (0->0), fcn. (480->10), ass. (0->65)
t103 = mrSges(6,1) + mrSges(7,1);
t102 = mrSges(6,2) + mrSges(7,2);
t31 = qJ(2) + pkin(10);
t24 = sin(t31);
t36 = sin(qJ(1));
t39 = cos(qJ(1));
t98 = g(1) * t39 + g(2) * t36;
t104 = t98 * t24;
t32 = qJ(4) + qJ(5);
t27 = cos(t32);
t37 = cos(qJ(4));
t28 = t37 * pkin(4);
t19 = pkin(5) * t27 + t28;
t17 = pkin(3) + t19;
t22 = t28 + pkin(3);
t26 = sin(t32);
t34 = sin(qJ(4));
t101 = -m(5) * pkin(3) - m(6) * t22 - m(7) * t17 - t37 * mrSges(5,1) + t34 * mrSges(5,2) + t102 * t26 - t103 * t27;
t90 = m(4) + m(5) + m(6) + m(7);
t100 = -mrSges(4,2) + mrSges(5,3) + mrSges(6,3) + mrSges(7,3);
t25 = cos(t31);
t35 = sin(qJ(2));
t38 = cos(qJ(2));
t99 = -t38 * mrSges(3,1) - t25 * mrSges(4,1) + t35 * mrSges(3,2) - t100 * t24;
t85 = m(6) * pkin(4);
t97 = mrSges(5,1) + t85;
t94 = mrSges(6,1) * t26 + t102 * t27;
t64 = t27 * t36;
t65 = t26 * t39;
t11 = -t25 * t65 + t64;
t63 = t27 * t39;
t66 = t26 * t36;
t12 = t25 * t63 + t66;
t93 = t102 * t12 - t103 * t11;
t10 = -t25 * t64 + t65;
t9 = t25 * t66 + t63;
t92 = -t10 * t102 + t103 * t9;
t91 = -m(3) * pkin(7) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t78 = pkin(4) * t34;
t18 = pkin(5) * t26 + t78;
t87 = -m(6) * t78 - m(7) * t18;
t86 = m(3) * pkin(1) + mrSges(2,1) - t99;
t84 = m(7) * pkin(5);
t40 = -pkin(9) - pkin(8);
t77 = pkin(8) * t24;
t74 = g(3) * t24;
t29 = t38 * pkin(2);
t72 = t18 * t36;
t30 = -qJ(6) + t40;
t68 = t24 * t30;
t67 = t24 * t40;
t62 = t34 * t36;
t61 = t34 * t39;
t60 = t36 * t37;
t59 = t37 * t39;
t55 = pkin(3) * t25 + t77;
t50 = t17 * t25 - t68;
t49 = t25 * t22 - t67;
t15 = -t25 * t61 + t60;
t13 = t25 * t62 + t59;
t33 = -qJ(3) - pkin(7);
t23 = t29 + pkin(1);
t16 = t25 * t59 + t62;
t14 = -t25 * t60 + t61;
t1 = [(-t62 * t85 - m(7) * t72 - t16 * mrSges(5,1) - t15 * mrSges(5,2) - t90 * (t39 * t23 - t36 * t33) - t103 * t12 - t102 * t11 + t91 * t36 + (-m(5) * t55 - m(6) * t49 - m(7) * t50 - t86) * t39) * g(2) + (-t14 * mrSges(5,1) - t13 * mrSges(5,2) - t102 * t9 - t103 * t10 + (t90 * t33 + t87 + t91) * t39 + (m(4) * t23 - m(5) * (-t23 - t55) - m(6) * (-t23 - t49) - m(7) * (-t23 - t50) + t86) * t36) * g(1) (mrSges(4,1) - t101) * t104 + (-m(4) * t29 - m(5) * (t29 + t77) - m(6) * (t29 - t67) - m(7) * (t29 - t68) + t99 + t101 * t25) * g(3) + (mrSges(3,2) * t38 + (-m(5) * pkin(8) + m(6) * t40 + m(7) * t30 - t100) * t25 + (t90 * pkin(2) + mrSges(3,1)) * t35) * t98 (-g(1) * t36 + g(2) * t39) * t90 (mrSges(5,1) * t34 + mrSges(7,1) * t26 + mrSges(5,2) * t37 - t87 + t94) * t74 + (-t14 * mrSges(5,2) - m(7) * (-t19 * t39 - t25 * t72) + t97 * t13 + t92) * g(2) + (t16 * mrSges(5,2) - m(7) * (-t18 * t25 * t39 + t19 * t36) - t97 * t15 + t93) * g(1) (-(-mrSges(7,1) - t84) * t26 + t94) * t74 + (t84 * t9 + t92) * g(2) + (-t11 * t84 + t93) * g(1) (g(3) * t25 - t104) * m(7)];
taug  = t1(:);
