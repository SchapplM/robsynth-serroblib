% Calculate Gravitation load on the joints for
% S6RRPPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3,theta4]';
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
% Datum: 2018-11-23 16:45
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRPPRP1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRP1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP1_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRP1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:45:10
% EndTime: 2018-11-23 16:45:11
% DurationCPUTime: 0.73s
% Computational Cost: add. (449->96), mult. (490->102), div. (0->0), fcn. (442->10), ass. (0->52)
t83 = mrSges(6,1) + mrSges(7,1);
t82 = -mrSges(6,2) + mrSges(7,3);
t73 = -m(6) - m(7);
t21 = qJ(2) + pkin(9);
t16 = sin(t21);
t27 = sin(qJ(1));
t29 = cos(qJ(1));
t75 = g(1) * t29 + g(2) * t27;
t81 = t75 * t16;
t80 = -mrSges(6,3) - mrSges(7,2);
t18 = cos(t21);
t26 = sin(qJ(2));
t28 = cos(qJ(2));
t22 = sin(pkin(10));
t23 = cos(pkin(10));
t37 = m(5) * pkin(3) + t23 * mrSges(5,1) - t22 * mrSges(5,2);
t79 = -t28 * mrSges(3,1) + t26 * mrSges(3,2) - t37 * t18;
t20 = pkin(10) + qJ(5);
t15 = sin(t20);
t17 = cos(t20);
t78 = t82 * t15 + t83 * t17;
t77 = m(5) * qJ(4) + mrSges(5,3);
t13 = pkin(4) * t23 + pkin(3);
t25 = -pkin(8) - qJ(4);
t76 = t18 * t13 - t16 * t25;
t74 = -m(3) * pkin(1) - mrSges(2,1) + t79;
t24 = -qJ(3) - pkin(7);
t72 = -m(3) * pkin(7) + m(5) * t24 - t22 * mrSges(5,1) - t23 * mrSges(5,2) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t42 = t18 * mrSges(4,1) - t16 * mrSges(4,2);
t71 = t42 + (mrSges(5,3) - t80) * t16;
t70 = m(7) * pkin(5) + t83;
t69 = m(7) * qJ(6) + t82;
t67 = pkin(2) * t26;
t66 = pkin(4) * t22;
t63 = g(3) * t16;
t19 = t28 * pkin(2);
t58 = t16 * t29;
t57 = t18 * t29;
t56 = t27 * t15;
t55 = t27 * t17;
t54 = t29 * t15;
t53 = qJ(4) * t16;
t52 = m(5) - t73;
t14 = t19 + pkin(1);
t8 = t29 * t14;
t49 = -t27 * t24 + t8;
t39 = pkin(5) * t17 + qJ(6) * t15;
t4 = t17 * t57 + t56;
t3 = t18 * t54 - t55;
t2 = t18 * t55 - t54;
t1 = t17 * t29 + t18 * t56;
t5 = [(-m(4) * t49 - m(5) * t8 + t80 * t58 + t73 * (t13 * t57 - t25 * t58 + t27 * t66 + t49) - t70 * t4 - t69 * t3 + (-t77 * t16 - t42 + t74) * t29 + t72 * t27) * g(2) + (t69 * t1 + t70 * t2 + (m(4) * t14 - m(5) * (-t14 - t53) + t71 + (-t14 - t76) * t73 - t74) * t27 + (t73 * (-t24 + t66) + m(4) * t24 + t72) * t29) * g(1) (-m(4) * t19 - m(5) * (t19 + t53) + t73 * (t19 + t76) - t71 + t79) * g(3) + t75 * (mrSges(3,1) * t26 + mrSges(3,2) * t28 + (m(4) + m(5)) * t67 + t73 * (-t18 * t25 - t67)) + ((-m(7) * t39 - t78) * g(3) + t75 * (mrSges(4,2) - t77 + t80)) * t18 + (mrSges(4,1) + t37 + m(6) * t13 - m(7) * (-t13 - t39) + t78) * t81 (-g(1) * t27 + g(2) * t29) * (m(4) + t52) (t18 * g(3) - t81) * t52 (t70 * t15 - t69 * t17) * t63 + (t70 * t1 - t69 * t2) * g(2) + (t70 * t3 - t69 * t4) * g(1) (-g(1) * t3 - g(2) * t1 - t15 * t63) * m(7)];
taug  = t5(:);
