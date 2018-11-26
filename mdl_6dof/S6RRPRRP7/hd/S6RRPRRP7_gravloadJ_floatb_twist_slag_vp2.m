% Calculate Gravitation load on the joints for
% S6RRPRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
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
% Datum: 2018-11-23 17:15
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRPRRP7_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP7_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP7_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP7_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP7_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP7_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:15:13
% EndTime: 2018-11-23 17:15:14
% DurationCPUTime: 1.06s
% Computational Cost: add. (409->104), mult. (969->130), div. (0->0), fcn. (1040->8), ass. (0->56)
t102 = m(7) * pkin(5) + mrSges(6,1) + mrSges(7,1);
t114 = m(7) * qJ(6) - mrSges(6,2) + mrSges(7,3);
t39 = sin(qJ(5));
t42 = cos(qJ(5));
t123 = t102 * t42 + t114 * t39 + mrSges(5,1);
t116 = mrSges(5,2) - mrSges(6,3) - mrSges(7,2);
t40 = sin(qJ(2));
t43 = cos(qJ(2));
t44 = cos(qJ(1));
t89 = sin(qJ(4));
t71 = t44 * t89;
t90 = cos(qJ(4));
t72 = t44 * t90;
t15 = -t40 * t71 - t43 * t72;
t16 = -t40 * t72 + t43 * t71;
t122 = t116 * t15 - t123 * t16;
t20 = t40 * t90 - t43 * t89;
t41 = sin(qJ(1));
t13 = t20 * t41;
t19 = t40 * t89 + t43 * t90;
t14 = t19 * t41;
t121 = -t116 * t14 + t123 * t13;
t120 = -t116 * t20 - t123 * t19;
t113 = -m(6) - m(7);
t100 = -pkin(2) - pkin(3);
t37 = t44 * pkin(7);
t33 = t40 * qJ(3);
t67 = -pkin(1) - t33;
t119 = (t100 * t43 + t67) * t41 - pkin(8) * t44 + t37;
t118 = (mrSges(3,1) + mrSges(4,1)) * t43 + (-mrSges(3,2) + mrSges(4,3)) * t40;
t117 = -mrSges(2,1) - t118;
t115 = mrSges(2,2) - mrSges(3,3) - mrSges(4,2) + mrSges(5,3);
t109 = -t19 * pkin(4) + pkin(9) * t20;
t108 = -t16 * pkin(4) - t15 * pkin(9);
t107 = t13 * pkin(4) + pkin(9) * t14;
t106 = g(1) * t44 + g(2) * t41;
t95 = g(3) * t20;
t36 = t43 * pkin(2);
t82 = t43 * t44;
t81 = t36 + t33;
t80 = t44 * pkin(1) + t41 * pkin(7);
t79 = qJ(3) * t43;
t78 = t43 * pkin(3) + t81;
t77 = t40 * t100;
t1 = t14 * t39 - t44 * t42;
t66 = pkin(2) * t82 + t44 * t33 + t80;
t2 = t14 * t42 + t39 * t44;
t25 = t41 * t79;
t56 = t41 * t77 + t25;
t26 = t44 * t79;
t55 = t44 * t77 + t26;
t52 = pkin(3) * t82 - pkin(8) * t41 + t66;
t51 = t43 * mrSges(4,3) + (-m(4) * pkin(2) - mrSges(4,1)) * t40;
t6 = -t15 * t42 - t41 * t39;
t5 = -t15 * t39 + t41 * t42;
t3 = [(-m(3) * t80 - m(4) * t66 - m(5) * t52 + t15 * mrSges(5,1) + t113 * (-t15 * pkin(4) + pkin(9) * t16 + t52) - t102 * t6 - t114 * t5 + t117 * t44 + t116 * t16 + t115 * t41) * g(2) + (t14 * mrSges(5,1) - t119 * m(5) + t113 * (-t14 * pkin(4) + t13 * pkin(9) + t119) + t102 * t2 + (-m(3) - m(4)) * t37 + t114 * t1 + (m(3) * pkin(1) - m(4) * (t67 - t36) - t117) * t41 + t116 * t13 + t115 * t44) * g(1), t106 * (mrSges(3,1) * t40 + mrSges(3,2) * t43) + (-m(4) * t25 - m(5) * t56 - t51 * t41 + t113 * (t56 - t107) + t121) * g(2) + (-m(4) * t26 - m(5) * t55 - t51 * t44 + t113 * (t55 - t108) + t122) * g(1) + (-m(4) * t81 - m(5) * t78 + t113 * (t78 - t109) - t118 + t120) * g(3) (t43 * g(3) - t106 * t40) * (m(4) + m(5) - t113) (t113 * t109 - t120) * g(3) + (t113 * t107 - t121) * g(2) + (t113 * t108 - t122) * g(1) (t102 * t39 - t114 * t42) * t95 + (t102 * t1 - t114 * t2) * g(2) + (t102 * t5 - t114 * t6) * g(1) (-g(1) * t5 - g(2) * t1 - t39 * t95) * m(7)];
taug  = t3(:);
