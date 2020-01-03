% Calculate Gravitation load on the joints for
% S5RRRPP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
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
% Datum: 2019-12-31 21:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRPP7_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP7_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP7_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP7_gravloadJ_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP7_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP7_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:03:57
% EndTime: 2019-12-31 21:03:59
% DurationCPUTime: 0.79s
% Computational Cost: add. (222->91), mult. (504->108), div. (0->0), fcn. (489->6), ass. (0->43)
t75 = mrSges(4,1) + mrSges(5,1);
t67 = mrSges(4,2) - mrSges(6,2) - mrSges(5,3);
t74 = mrSges(5,2) + mrSges(4,3);
t73 = mrSges(6,3) - t74;
t23 = sin(qJ(3));
t26 = cos(qJ(3));
t72 = -t67 * t23 + t75 * t26;
t24 = sin(qJ(2));
t59 = g(3) * t24;
t71 = m(5) + m(6);
t70 = mrSges(2,2) - mrSges(3,3);
t25 = sin(qJ(1));
t28 = cos(qJ(1));
t68 = g(1) * t28 + g(2) * t25;
t27 = cos(qJ(2));
t39 = m(6) * (-pkin(3) - pkin(4)) - mrSges(6,1);
t48 = qJ(4) * t23;
t43 = -pkin(2) - t48;
t66 = (m(6) * qJ(5) + t73) * t27 + (-m(6) * t43 - t39 * t26 - m(5) * (-pkin(3) * t26 + t43) + m(4) * pkin(2) + t72) * t24;
t38 = t27 * mrSges(3,1) - t24 * mrSges(3,2);
t65 = t74 * t24 + t38;
t64 = m(6) * pkin(4) + mrSges(6,1);
t63 = pkin(7) * (-m(4) - t71);
t62 = t64 + t75;
t18 = t24 * pkin(7);
t20 = t27 * pkin(2);
t55 = t24 * t28;
t54 = t25 * t27;
t53 = t26 * t27;
t52 = t27 * t28;
t51 = t28 * t23;
t50 = t20 + t18;
t49 = t28 * pkin(1) + t25 * pkin(6);
t47 = qJ(5) * t24;
t46 = -pkin(1) - t20;
t42 = pkin(3) * t53 + t27 * t48 + t50;
t41 = pkin(2) * t52 + pkin(7) * t55 + t49;
t21 = t28 * pkin(6);
t8 = t25 * t23 + t26 * t52;
t7 = -t25 * t26 + t27 * t51;
t6 = t25 * t53 - t51;
t5 = t23 * t54 + t26 * t28;
t1 = [(-m(3) * t49 - m(4) * t41 - t71 * (t8 * pkin(3) + t7 * qJ(4) + t41) + (m(6) * t47 - mrSges(2,1) - t38) * t28 + t70 * t25 - t62 * t8 + t67 * t7 + t73 * t55) * g(2) + (-t71 * (-t6 * pkin(3) - qJ(4) * t5 + t21) + t70 * t28 + (-m(3) - m(4)) * t21 + t62 * t6 - t67 * t5 + (mrSges(2,1) + m(3) * pkin(1) - m(6) * t46 - (m(6) * (-pkin(7) + qJ(5)) + mrSges(6,3)) * t24 + (-m(4) - m(5)) * (t46 - t18) + t65) * t25) * g(1), t68 * (mrSges(3,1) * t24 + mrSges(3,2) * t27) + (t66 * t25 + t54 * t63) * g(2) + (t66 * t28 + t52 * t63) * g(1) + (-m(4) * t50 - m(5) * t42 - m(6) * (t42 - t47) + t24 * mrSges(6,3) + (-t64 * t26 - t72) * t27 - t65) * g(3), (-t71 * (-t5 * pkin(3) + qJ(4) * t6) + t67 * t6 + t62 * t5) * g(2) + (-t71 * (-t7 * pkin(3) + qJ(4) * t8) + t67 * t8 + t62 * t7) * g(1) + ((-t71 * qJ(4) + t67) * t26 + (m(5) * pkin(3) - t39 + t75) * t23) * t59, t71 * (-g(1) * t7 - g(2) * t5 - t23 * t59), (-g(3) * t27 + t68 * t24) * m(6)];
taug = t1(:);
