% Calculate Gravitation load on the joints for
% S6RRRRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2018-11-23 18:21
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRRRPR11_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR11_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR11_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR11_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR11_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR11_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:20:31
% EndTime: 2018-11-23 18:20:32
% DurationCPUTime: 1.50s
% Computational Cost: add. (1767->154), mult. (1981->193), div. (0->0), fcn. (1987->18), ass. (0->74)
t57 = qJ(4) + pkin(12);
t52 = cos(t57);
t63 = cos(qJ(4));
t54 = t63 * pkin(4);
t38 = pkin(5) * t52 + t54;
t34 = pkin(3) + t38;
t53 = qJ(6) + t57;
t48 = sin(t53);
t49 = cos(t53);
t50 = t54 + pkin(3);
t51 = sin(t57);
t59 = sin(qJ(4));
t121 = m(5) * pkin(3) + m(6) * t50 + m(7) * t34 + t63 * mrSges(5,1) + t52 * mrSges(6,1) + t49 * mrSges(7,1) - t59 * mrSges(5,2) - t51 * mrSges(6,2) - t48 * mrSges(7,2) + mrSges(4,1);
t58 = -qJ(5) - pkin(10);
t69 = mrSges(4,2) - m(5) * pkin(10) - mrSges(5,3) + m(7) * (-pkin(11) + t58) - mrSges(7,3) + m(6) * t58 - mrSges(6,3);
t128 = m(6) + m(7);
t134 = -m(4) - m(5);
t101 = t128 - t134;
t60 = sin(qJ(3));
t64 = cos(qJ(3));
t129 = pkin(2) * t101 + t121 * t64 - t69 * t60 + mrSges(3,1);
t130 = -t51 * mrSges(6,1) - t48 * mrSges(7,1) - t63 * mrSges(5,2) - t52 * mrSges(6,2) - t49 * mrSges(7,2);
t113 = pkin(4) * t59;
t126 = mrSges(3,2) - mrSges(4,3);
t37 = pkin(5) * t51 + t113;
t117 = -t59 * mrSges(5,1) + t126 + t134 * pkin(9) - m(6) * (pkin(9) + t113) - m(7) * (pkin(9) + t37) + t130;
t125 = m(6) * pkin(4) + mrSges(5,1);
t122 = m(7) * t37 + t101 * pkin(9) - t126;
t62 = sin(qJ(1));
t65 = cos(qJ(2));
t110 = t62 * t65;
t111 = cos(qJ(1));
t103 = pkin(6) + qJ(2);
t91 = sin(t103);
t85 = t91 / 0.2e1;
t104 = pkin(6) - qJ(2);
t92 = sin(t104);
t78 = t85 - t92 / 0.2e1;
t26 = t111 * t78 + t110;
t105 = sin(pkin(6));
t88 = t111 * t105;
t14 = t26 * t64 - t60 * t88;
t61 = sin(qJ(2));
t87 = cos(t103) / 0.2e1;
t93 = cos(t104);
t71 = t93 / 0.2e1 + t87;
t25 = -t111 * t71 + t61 * t62;
t115 = (-t14 * t48 + t25 * t49) * mrSges(7,1) + (-t14 * t49 - t25 * t48) * mrSges(7,2);
t98 = t111 * t65;
t29 = -t62 * t78 + t98;
t94 = t62 * t105;
t18 = t29 * t64 + t60 * t94;
t28 = t111 * t61 + t62 * t71;
t5 = -t18 * t48 + t28 * t49;
t6 = t18 * t49 + t28 * t48;
t114 = t5 * mrSges(7,1) - t6 * mrSges(7,2);
t106 = cos(pkin(6));
t36 = t87 - t93 / 0.2e1;
t24 = t106 * t60 - t36 * t64;
t86 = t92 / 0.2e1;
t35 = t85 + t86;
t108 = (-t24 * t48 - t35 * t49) * mrSges(7,1) + (-t24 * t49 + t35 * t48) * mrSges(7,2);
t107 = t111 * pkin(1) + pkin(8) * t94;
t102 = t29 * pkin(2) + t107;
t95 = -pkin(1) * t62 + pkin(8) * t88;
t9 = -t18 * t59 + t28 * t63;
t79 = t86 - t91 / 0.2e1;
t13 = t26 * t60 + t64 * t88;
t23 = -t106 * t64 - t36 * t60;
t17 = t29 * t60 - t64 * t94;
t10 = t18 * t63 + t28 * t59;
t8 = t18 * t52 + t28 * t51;
t7 = -t18 * t51 + t28 * t52;
t1 = [(-t111 * mrSges(2,1) - m(3) * t107 - t29 * mrSges(3,1) - m(4) * t102 - t18 * mrSges(4,1) - m(5) * (pkin(3) * t18 + t102) - t10 * mrSges(5,1) - t9 * mrSges(5,2) - m(6) * (t18 * t50 + t102) - t8 * mrSges(6,1) - t7 * mrSges(6,2) - m(7) * (t18 * t34 + t102) - t6 * mrSges(7,1) - t5 * mrSges(7,2) + (-t105 * mrSges(3,3) + mrSges(2,2)) * t62 + (-m(6) * t113 - t122) * t28 + t69 * t17) * g(2) + (t62 * mrSges(2,1) + t111 * mrSges(2,2) - m(3) * t95 + t26 * mrSges(3,1) - mrSges(3,3) * t88 + t121 * t14 + (t125 * t59 + t122 - t130) * t25 - t69 * t13 + t101 * (t26 * pkin(2) - t95)) * g(1) (-t117 * t36 - t129 * t35) * g(3) + (t117 * (-t111 * t79 + t110) + t129 * t25) * g(2) + (t117 * (t62 * t79 + t98) + t129 * t28) * g(1) (t121 * t23 + t69 * t24) * g(3) + (t121 * t13 + t69 * t14) * g(2) + (t121 * t17 + t69 * t18) * g(1) (-(-t24 * t63 + t35 * t59) * mrSges(5,2) - (-t24 * t51 - t35 * t52) * mrSges(6,1) - (-t24 * t52 + t35 * t51) * mrSges(6,2) - m(7) * (-t24 * t37 - t35 * t38) - t108 - t125 * (-t24 * t59 - t35 * t63)) * g(3) + (-(-t14 * t63 - t25 * t59) * mrSges(5,2) - (-t14 * t51 + t25 * t52) * mrSges(6,1) - (-t14 * t52 - t25 * t51) * mrSges(6,2) - m(7) * (-t14 * t37 + t25 * t38) - t115 - t125 * (-t14 * t59 + t25 * t63)) * g(2) + (t10 * mrSges(5,2) - t7 * mrSges(6,1) + t8 * mrSges(6,2) - m(7) * (-t18 * t37 + t28 * t38) - t114 - t125 * t9) * g(1), t128 * (-g(1) * t17 - g(2) * t13 - g(3) * t23) -g(1) * t114 - g(2) * t115 - g(3) * t108];
taug  = t1(:);
