% Calculate Gravitation load on the joints for
% S6RRRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4,theta5]';
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
% Datum: 2018-11-23 17:33
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRRPPR1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR1_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:32:52
% EndTime: 2018-11-23 17:32:53
% DurationCPUTime: 0.74s
% Computational Cost: add. (538->124), mult. (475->131), div. (0->0), fcn. (402->12), ass. (0->68)
t42 = qJ(2) + qJ(3);
t36 = pkin(10) + t42;
t29 = sin(t36);
t30 = cos(t36);
t43 = sin(pkin(11));
t87 = t43 * mrSges(6,2);
t41 = pkin(11) + qJ(6);
t34 = sin(t41);
t94 = t34 * mrSges(7,2);
t126 = t30 * (-mrSges(6,3) - mrSges(7,3)) + (-t87 - t94) * t29;
t37 = sin(t42);
t38 = cos(t42);
t123 = mrSges(4,1) * t37 + mrSges(5,1) * t29 + mrSges(4,2) * t38 + mrSges(5,2) * t30;
t35 = cos(t41);
t91 = t35 * mrSges(7,1);
t44 = cos(pkin(11));
t99 = mrSges(6,1) * t44;
t122 = (t91 + t99) * t29;
t121 = -t87 + t99;
t47 = sin(qJ(1));
t49 = cos(qJ(1));
t120 = g(1) * t49 + g(2) * t47;
t119 = -t38 * mrSges(4,1) - t30 * mrSges(5,1) + t37 * mrSges(4,2) + (mrSges(5,2) - mrSges(7,3)) * t29;
t108 = m(6) + m(7);
t31 = pkin(5) * t44 + pkin(4);
t45 = -pkin(9) - qJ(5);
t64 = -t29 * t45 + t30 * t31;
t105 = pkin(4) * t29;
t106 = pkin(3) * t37;
t63 = -t29 * t31 - t30 * t45;
t117 = -m(7) * (t63 - t106) - m(6) * (-t105 - t106) + t122;
t115 = t126 * t49;
t114 = t126 * t47;
t25 = t29 * mrSges(6,3);
t112 = -t25 + t119 + (-t91 + t94 - t121) * t30;
t111 = m(6) * t105 - m(7) * t63 + t122;
t50 = -pkin(8) - pkin(7);
t40 = -qJ(4) + t50;
t110 = -m(3) * pkin(7) + m(4) * t50 + m(6) * t40 - t43 * mrSges(6,1) - t44 * mrSges(6,2) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3);
t48 = cos(qJ(2));
t39 = t48 * pkin(2);
t46 = sin(qJ(2));
t69 = t48 * mrSges(3,1) - t46 * mrSges(3,2);
t109 = -(m(6) * pkin(4) + t121) * t30 - mrSges(2,1) - m(4) * (t39 + pkin(1)) - m(3) * pkin(1) - t69 + t119;
t107 = pkin(2) * t46;
t32 = pkin(3) * t38;
t104 = pkin(5) * t43;
t93 = t34 * t47;
t92 = t34 * t49;
t90 = t35 * t47;
t89 = t35 * t49;
t84 = t32 + t39;
t83 = qJ(5) * t30;
t23 = t29 * qJ(5);
t77 = t30 * pkin(4) + t23 + t32;
t71 = t32 + t64;
t18 = t49 * t83;
t17 = t47 * t83;
t15 = -t106 - t107;
t14 = pkin(1) + t84;
t9 = t49 * t15;
t8 = t47 * t15;
t7 = t49 * t14;
t4 = t30 * t89 + t93;
t3 = -t30 * t92 + t90;
t2 = -t30 * t90 + t92;
t1 = t30 * t93 + t89;
t5 = [(-m(6) * t7 - t4 * mrSges(7,1) - t3 * mrSges(7,2) + (-m(5) - m(7)) * (-t40 * t47 + t7) + (-m(7) * t104 + t110) * t47 + (-(m(6) * qJ(5) + mrSges(6,3)) * t29 - m(7) * t64 + t109) * t49) * g(2) + (-t2 * mrSges(7,1) - t1 * mrSges(7,2) + (m(5) * t40 - m(7) * (-t40 + t104) + t110) * t49 + (m(5) * t14 - m(6) * (-t14 - t23) + t25 - m(7) * (-t14 - t64) - t109) * t47) * g(1) (-m(6) * (t17 + t8) - m(7) * t8 + t111 * t47 + t114) * g(2) + (-m(6) * (t18 + t9) - m(7) * t9 + t111 * t49 + t115) * g(1) + (-t69 - m(4) * t39 - m(5) * t84 - m(6) * (t39 + t77) - m(7) * (t39 + t71) + t112) * g(3) + t120 * (m(4) * t107 - m(5) * t15 + mrSges(3,1) * t46 + mrSges(3,2) * t48 + t123) (-m(6) * t17 + t117 * t47 + t114) * g(2) + (-m(6) * t18 + t117 * t49 + t115) * g(1) + (-m(5) * t32 - m(6) * t77 - m(7) * t71 + t112) * g(3) + t120 * (m(5) * t106 + t123) (-g(1) * t47 + g(2) * t49) * (m(5) + t108) (t30 * g(3) - t29 * t120) * t108, -g(1) * (mrSges(7,1) * t3 - mrSges(7,2) * t4) - g(2) * (-mrSges(7,1) * t1 + mrSges(7,2) * t2) - g(3) * (-mrSges(7,1) * t34 - mrSges(7,2) * t35) * t29];
taug  = t5(:);
