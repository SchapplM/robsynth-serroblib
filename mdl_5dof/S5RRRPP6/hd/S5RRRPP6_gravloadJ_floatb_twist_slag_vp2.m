% Calculate Gravitation load on the joints for
% S5RRRPP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
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
% Datum: 2019-12-31 21:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRPP6_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP6_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP6_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP6_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP6_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP6_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:59:59
% EndTime: 2019-12-31 21:00:02
% DurationCPUTime: 0.72s
% Computational Cost: add. (290->89), mult. (454->107), div. (0->0), fcn. (425->8), ass. (0->48)
t80 = mrSges(5,1) + mrSges(6,1);
t79 = -mrSges(5,2) + mrSges(6,3);
t71 = m(5) + m(6);
t73 = pkin(3) * t71 + mrSges(4,1);
t78 = -mrSges(5,3) - mrSges(6,2);
t23 = qJ(3) + pkin(8);
t18 = sin(t23);
t19 = cos(t23);
t25 = sin(qJ(3));
t28 = cos(qJ(3));
t77 = m(4) * pkin(2) + t28 * mrSges(4,1) - t25 * mrSges(4,2) + t79 * t18 + t80 * t19;
t26 = sin(qJ(2));
t75 = t26 * mrSges(4,3) + mrSges(2,1);
t27 = sin(qJ(1));
t29 = cos(qJ(2));
t30 = cos(qJ(1));
t52 = t29 * t30;
t7 = -t25 * t52 + t27 * t28;
t74 = g(1) * t30 + g(2) * t27;
t72 = -m(3) - m(4);
t70 = mrSges(2,2) - mrSges(3,3);
t42 = t29 * mrSges(3,1) - t26 * mrSges(3,2);
t69 = t78 * t26 - t42;
t67 = m(6) * pkin(4) + t80;
t66 = m(6) * qJ(5) + t79;
t62 = g(3) * t26;
t24 = -qJ(4) - pkin(7);
t61 = t24 * t26;
t60 = t25 * t30;
t56 = t26 * t30;
t55 = t27 * t25;
t53 = t27 * t29;
t17 = pkin(3) * t28 + pkin(2);
t10 = t29 * t17;
t51 = t30 * t18;
t50 = t30 * pkin(1) + t27 * pkin(6);
t47 = m(4) * pkin(7) + mrSges(4,3);
t44 = pkin(2) * t29 + pkin(7) * t26;
t38 = pkin(4) * t19 + qJ(5) * t18;
t5 = t25 * t53 + t28 * t30;
t21 = t30 * pkin(6);
t8 = t28 * t52 + t55;
t6 = -t28 * t53 + t60;
t4 = t27 * t18 + t19 * t52;
t3 = -t27 * t19 + t29 * t51;
t2 = t19 * t53 - t51;
t1 = t18 * t53 + t19 * t30;
t9 = [(-t8 * mrSges(4,1) - t7 * mrSges(4,2) + t78 * t56 + t72 * t50 - t71 * (pkin(3) * t55 + t17 * t52 - t24 * t56 + t50) - t67 * t4 - t66 * t3 + t70 * t27 + (-m(4) * t44 - t42 - t75) * t30) * g(2) + (-t6 * mrSges(4,1) - t5 * mrSges(4,2) - t71 * (pkin(3) * t60 + t27 * t61 + t21) + t70 * t30 + t72 * t21 + t67 * t2 + t66 * t1 + (m(3) * pkin(1) - m(4) * (-pkin(1) - t44) - t71 * (-pkin(1) - t10) - t69 + t75) * t27) * g(1), (-t71 * (t10 - t61) + t69) * g(3) + ((-m(6) * t38 - t77) * g(3) + t74 * (t71 * t24 + mrSges(3,2) - t47 + t78)) * t29 + (-t47 * g(3) + t74 * (mrSges(3,1) + m(5) * t17 - m(6) * (-t17 - t38) + t77)) * t26, (mrSges(4,2) * t28 + t67 * t18 - t66 * t19 + t73 * t25) * t62 + (-t6 * mrSges(4,2) + t67 * t1 - t66 * t2 + t73 * t5) * g(2) + (t8 * mrSges(4,2) + t67 * t3 - t66 * t4 - t73 * t7) * g(1), (t29 * g(3) - t26 * t74) * t71, (-g(1) * t3 - g(2) * t1 - t18 * t62) * m(6)];
taug = t9(:);
