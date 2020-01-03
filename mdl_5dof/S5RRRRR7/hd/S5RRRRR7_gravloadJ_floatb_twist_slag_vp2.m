% Calculate Gravitation load on the joints for
% S5RRRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2019-12-31 22:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRRR7_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR7_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR7_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR7_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR7_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR7_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:21:22
% EndTime: 2019-12-31 22:21:24
% DurationCPUTime: 0.51s
% Computational Cost: add. (377->91), mult. (358->98), div. (0->0), fcn. (300->10), ass. (0->54)
t29 = qJ(2) + qJ(3);
t26 = qJ(4) + t29;
t21 = sin(t26);
t22 = cos(t26);
t30 = sin(qJ(5));
t68 = t30 * mrSges(6,2);
t94 = t21 * t68 + t22 * (m(6) * pkin(9) + mrSges(6,3));
t86 = t22 * pkin(4) + t21 * pkin(9);
t93 = m(6) * t86;
t92 = -t22 * mrSges(5,1) + (mrSges(5,2) - mrSges(6,3)) * t21;
t24 = sin(t29);
t25 = cos(t29);
t72 = mrSges(5,2) * t22;
t91 = mrSges(4,1) * t24 + mrSges(5,1) * t21 + mrSges(4,2) * t25 + t72;
t32 = sin(qJ(1));
t35 = cos(qJ(1));
t90 = g(1) * t35 + g(2) * t32;
t87 = m(5) + m(6);
t55 = t25 * mrSges(4,1) - mrSges(4,2) * t24;
t33 = cos(qJ(5));
t74 = mrSges(6,1) * t33;
t85 = -(-t68 + t74) * t22 + t92;
t82 = -t55 + t85;
t36 = -pkin(7) - pkin(6);
t81 = -m(3) * pkin(6) + m(4) * t36 + mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3);
t34 = cos(qJ(2));
t27 = t34 * pkin(2);
t31 = sin(qJ(2));
t50 = t34 * mrSges(3,1) - t31 * mrSges(3,2);
t80 = mrSges(2,1) + m(4) * (t27 + pkin(1)) + t55 + m(3) * pkin(1) + t50 - t92;
t79 = pkin(2) * t31;
t78 = pkin(3) * t24;
t20 = pkin(3) * t25;
t77 = pkin(4) * t21;
t67 = t30 * t32;
t66 = t30 * t35;
t65 = t32 * t33;
t64 = t33 * t35;
t62 = t20 + t27;
t60 = t21 * t74;
t59 = t20 + t86;
t53 = t94 * t32;
t52 = t94 * t35;
t9 = -t78 - t79;
t39 = m(6) * (t9 - t77) - t60;
t38 = m(6) * (-t77 - t78) - t60;
t37 = t72 + (m(6) * pkin(4) + mrSges(5,1) + t74) * t21;
t28 = -pkin(8) + t36;
t8 = pkin(1) + t62;
t4 = t22 * t64 + t67;
t3 = -t22 * t66 + t65;
t2 = -t22 * t65 + t66;
t1 = t22 * t67 + t64;
t5 = [(-t4 * mrSges(6,1) - t3 * mrSges(6,2) - t87 * (-t28 * t32 + t35 * t8) + t81 * t32 + (-t80 - t93) * t35) * g(2) + (-t2 * mrSges(6,1) - t1 * mrSges(6,2) + (t87 * t28 + t81) * t35 + (m(5) * t8 - m(6) * (-t86 - t8) + t80) * t32) * g(1), -g(1) * (t35 * t39 + t52) - g(2) * (t32 * t39 + t53) + (-t50 - m(4) * t27 - m(5) * t62 - m(6) * (t27 + t59) + t82) * g(3) + t90 * (m(4) * t79 - m(5) * t9 + mrSges(3,1) * t31 + mrSges(3,2) * t34 + t91), -g(1) * (t35 * t38 + t52) - g(2) * (t32 * t38 + t53) + (-m(5) * t20 - m(6) * t59 + t82) * g(3) + (m(5) * t78 + t91) * t90, (t85 - t93) * g(3) + (t32 * t37 - t53) * g(2) + (t35 * t37 - t52) * g(1), -g(1) * (mrSges(6,1) * t3 - mrSges(6,2) * t4) - g(2) * (-mrSges(6,1) * t1 + mrSges(6,2) * t2) - g(3) * (-mrSges(6,1) * t30 - mrSges(6,2) * t33) * t21];
taug = t5(:);
