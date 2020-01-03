% Calculate Gravitation load on the joints for
% S5RRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-31 21:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRPR5_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR5_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR5_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR5_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR5_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR5_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:13:01
% EndTime: 2019-12-31 21:13:03
% DurationCPUTime: 0.50s
% Computational Cost: add. (318->83), mult. (316->94), div. (0->0), fcn. (263->10), ass. (0->49)
t29 = qJ(2) + qJ(3);
t24 = pkin(9) + t29;
t20 = sin(t24);
t21 = cos(t24);
t30 = sin(qJ(5));
t68 = t30 * mrSges(6,2);
t90 = t20 * t68 + t21 * (m(6) * pkin(8) + mrSges(6,3));
t25 = sin(t29);
t26 = cos(t29);
t89 = mrSges(4,1) * t25 + mrSges(5,1) * t20 + mrSges(4,2) * t26 + mrSges(5,2) * t21;
t32 = sin(qJ(1));
t35 = cos(qJ(1));
t88 = g(1) * t35 + g(2) * t32;
t87 = -t26 * mrSges(4,1) - t21 * mrSges(5,1) + mrSges(4,2) * t25 + (mrSges(5,2) - mrSges(6,3)) * t20;
t84 = m(5) + m(6);
t51 = t21 * pkin(4) + t20 * pkin(8);
t33 = cos(qJ(5));
t64 = t33 * mrSges(6,1);
t81 = -(t64 - t68) * t21 + t87;
t36 = -pkin(7) - pkin(6);
t79 = -m(3) * pkin(6) + m(4) * t36 + mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3);
t34 = cos(qJ(2));
t27 = t34 * pkin(2);
t31 = sin(qJ(2));
t49 = t34 * mrSges(3,1) - t31 * mrSges(3,2);
t78 = mrSges(2,1) + m(4) * (t27 + pkin(1)) + m(3) * pkin(1) + t49 - t87;
t77 = pkin(2) * t31;
t76 = pkin(3) * t25;
t22 = pkin(3) * t26;
t75 = pkin(4) * t20;
t67 = t30 * t32;
t66 = t30 * t35;
t65 = t32 * t33;
t63 = t33 * t35;
t62 = t22 + t27;
t60 = t20 * t64;
t59 = t22 + t51;
t53 = t90 * t32;
t52 = t90 * t35;
t9 = -t76 - t77;
t38 = m(6) * (t9 - t75) - t60;
t37 = m(6) * (-t75 - t76) - t60;
t28 = -qJ(4) + t36;
t8 = pkin(1) + t62;
t4 = t21 * t63 + t67;
t3 = -t21 * t66 + t65;
t2 = -t21 * t65 + t66;
t1 = t21 * t67 + t63;
t5 = [(-t4 * mrSges(6,1) - t3 * mrSges(6,2) - t84 * (-t28 * t32 + t35 * t8) + t79 * t32 + (-m(6) * t51 - t78) * t35) * g(2) + (-t2 * mrSges(6,1) - t1 * mrSges(6,2) + (t84 * t28 + t79) * t35 + (m(5) * t8 - m(6) * (-t51 - t8) + t78) * t32) * g(1), -g(1) * (t38 * t35 + t52) - g(2) * (t38 * t32 + t53) + (-t49 - m(4) * t27 - m(5) * t62 - m(6) * (t27 + t59) + t81) * g(3) + t88 * (m(4) * t77 - m(5) * t9 + mrSges(3,1) * t31 + mrSges(3,2) * t34 + t89), -g(1) * (t37 * t35 + t52) - g(2) * (t37 * t32 + t53) + (-m(5) * t22 - m(6) * t59 + t81) * g(3) + (m(5) * t76 + t89) * t88, t84 * (-g(1) * t32 + g(2) * t35), -g(1) * (mrSges(6,1) * t3 - mrSges(6,2) * t4) - g(2) * (-mrSges(6,1) * t1 + mrSges(6,2) * t2) - g(3) * (-mrSges(6,1) * t30 - mrSges(6,2) * t33) * t20];
taug = t5(:);
