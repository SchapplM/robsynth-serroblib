% Calculate Gravitation load on the joints for
% S6RRPPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta5]';
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
% Datum: 2018-11-23 16:43
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRPPPR2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR2_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPPR2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:42:48
% EndTime: 2018-11-23 16:42:49
% DurationCPUTime: 0.74s
% Computational Cost: add. (349->97), mult. (414->99), div. (0->0), fcn. (351->10), ass. (0->50)
t74 = m(6) + m(7);
t48 = m(5) + t74;
t84 = m(4) + t48;
t83 = mrSges(4,1) - mrSges(5,2);
t82 = -mrSges(4,2) + mrSges(5,3);
t19 = pkin(10) + qJ(6);
t14 = sin(t19);
t16 = cos(t19);
t21 = sin(pkin(10));
t22 = cos(pkin(10));
t36 = t21 * mrSges(6,1) + t22 * mrSges(6,2);
t62 = pkin(5) * t21;
t81 = -m(7) * t62 - t14 * mrSges(7,1) - t16 * mrSges(7,2) - t36;
t24 = -pkin(8) - qJ(5);
t80 = -m(6) * (-pkin(3) - qJ(5)) + mrSges(6,3) - m(7) * (-pkin(3) + t24) + mrSges(7,3);
t26 = sin(qJ(1));
t28 = cos(qJ(1));
t70 = g(1) * t28 + g(2) * t26;
t20 = qJ(2) + pkin(9);
t15 = sin(t20);
t17 = cos(t20);
t25 = sin(qJ(2));
t27 = cos(qJ(2));
t77 = -t27 * mrSges(3,1) + t25 * mrSges(3,2) - t82 * t15 - t83 * t17;
t75 = m(5) + m(7);
t10 = t15 * qJ(4);
t53 = t17 * t28;
t73 = pkin(3) * t53 + t28 * t10;
t72 = -t17 * pkin(3) - t10;
t67 = -m(3) * pkin(1) - mrSges(2,1) + t77;
t23 = -qJ(3) - pkin(7);
t65 = -m(6) * (pkin(4) - t23) - m(7) * (pkin(5) * t22 + pkin(4)) - t22 * mrSges(6,1) + t21 * mrSges(6,2) - mrSges(5,1) + mrSges(2,2) - mrSges(4,3) - m(3) * pkin(7) - mrSges(3,3);
t59 = g(3) * t17;
t18 = t27 * pkin(2);
t57 = t15 * t28;
t56 = t16 * t28;
t55 = t17 * mrSges(7,3);
t54 = t17 * t24;
t52 = t26 * t14;
t51 = t26 * t16;
t49 = qJ(5) * t17;
t47 = t18 - t72;
t13 = t18 + pkin(1);
t9 = t28 * t13;
t44 = -t26 * t23 + t9;
t4 = -t15 * t52 + t56;
t3 = t14 * t28 + t15 * t51;
t2 = t14 * t57 + t51;
t1 = t15 * t56 - t52;
t5 = [(-m(4) * t44 - m(6) * (t9 + t73) - mrSges(6,3) * t53 - t2 * mrSges(7,1) - t1 * mrSges(7,2) - t36 * t57 - t75 * (t44 + t73) + t65 * t26 + (-m(6) * t49 - m(7) * (t15 * t62 - t54) - t55 + t67) * t28) * g(2) + (-t4 * mrSges(7,1) + t3 * mrSges(7,2) + ((m(4) + t75) * t23 + t65) * t28 + (-m(5) * t72 + t80 * t17 + (m(6) * qJ(4) + t36 - m(7) * (-qJ(4) - t62)) * t15 + t84 * t13 - t67) * t26) * g(1) (-m(4) * t18 - m(5) * t47 - m(6) * (t47 + t49) - t17 * mrSges(6,3) - m(7) * (t47 - t54) - t55 + t81 * t15 + t77) * g(3) + (mrSges(3,2) * t27 + (m(5) * pkin(3) + t80 + t83) * t15 + (-qJ(4) * t48 + t81 - t82) * t17 + (t84 * pkin(2) + mrSges(3,1)) * t25) * t70 (-g(1) * t26 + g(2) * t28) * t84 (-t70 * t15 + t59) * t48 (-t15 * g(3) - t70 * t17) * t74, -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * (mrSges(7,1) * t3 + mrSges(7,2) * t4) - (-mrSges(7,1) * t16 + mrSges(7,2) * t14) * t59];
taug  = t5(:);
