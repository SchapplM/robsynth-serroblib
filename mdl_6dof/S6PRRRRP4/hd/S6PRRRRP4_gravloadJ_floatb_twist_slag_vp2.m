% Calculate Gravitation load on the joints for
% S6PRRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
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
% Datum: 2018-11-23 15:30
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6PRRRRP4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP4_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRP4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:29:41
% EndTime: 2018-11-23 15:29:42
% DurationCPUTime: 1.08s
% Computational Cost: add. (1566->118), mult. (1758->155), div. (0->0), fcn. (1761->16), ass. (0->76)
t153 = mrSges(6,2) - mrSges(7,3);
t140 = -m(7) * qJ(6) + t153;
t154 = mrSges(6,1) + mrSges(7,1);
t141 = -m(7) * pkin(5) - t154;
t74 = qJ(4) + qJ(5);
t72 = sin(t74);
t73 = cos(t74);
t78 = sin(qJ(4));
t81 = cos(qJ(4));
t143 = -m(5) * pkin(3) - t81 * mrSges(5,1) + t78 * mrSges(5,2) + t140 * t72 + t141 * t73 - mrSges(4,1);
t142 = -m(5) * pkin(9) + mrSges(4,2) - mrSges(7,2) - mrSges(5,3) - mrSges(6,3);
t145 = -m(6) - m(7);
t71 = pkin(4) * t81 + pkin(3);
t79 = sin(qJ(3));
t82 = cos(qJ(3));
t84 = -pkin(10) - pkin(9);
t146 = t145 * (-t71 * t82 + t79 * t84) - t142 * t79 + mrSges(3,1) - t143 * t82;
t155 = -m(4) + t145;
t119 = pkin(6) - qJ(2);
t105 = cos(t119);
t118 = pkin(6) + qJ(2);
t99 = cos(t118) / 0.2e1;
t67 = t99 - t105 / 0.2e1;
t77 = cos(pkin(6));
t57 = -t67 * t82 + t77 * t79;
t103 = sin(t118);
t97 = t103 / 0.2e1;
t104 = sin(t119);
t98 = t104 / 0.2e1;
t66 = t97 + t98;
t151 = -t57 * t78 - t66 * t81;
t75 = sin(pkin(11));
t76 = sin(pkin(6));
t128 = t75 * t76;
t120 = cos(pkin(11));
t83 = cos(qJ(2));
t109 = t120 * t83;
t94 = t97 - t104 / 0.2e1;
t54 = -t75 * t94 + t109;
t35 = t79 * t128 + t54 * t82;
t80 = sin(qJ(2));
t88 = t105 / 0.2e1 + t99;
t53 = t120 * t80 + t75 * t88;
t150 = -t35 * t78 + t53 * t81;
t110 = t76 * t120;
t127 = t75 * t83;
t51 = t120 * t94 + t127;
t33 = -t79 * t110 + t51 * t82;
t50 = -t120 * t88 + t75 * t80;
t149 = -t33 * t78 + t50 * t81;
t147 = -m(5) * pkin(8) - t81 * mrSges(5,2) - t140 * t73 + t141 * t72 + mrSges(3,2) - mrSges(4,3) + (pkin(4) * t145 - mrSges(5,1)) * t78;
t11 = t33 * t72 - t50 * t73;
t12 = t33 * t73 + t50 * t72;
t117 = t154 * t11 + t153 * t12;
t13 = t35 * t72 - t53 * t73;
t14 = t35 * t73 + t53 * t72;
t116 = t154 * t13 + t153 * t14;
t111 = -t11 * pkin(5) + qJ(6) * t12;
t108 = -t13 * pkin(5) + qJ(6) * t14;
t24 = t57 * t72 + t66 * t73;
t25 = t57 * t73 - t66 * t72;
t107 = -t24 * pkin(5) + qJ(6) * t25;
t106 = t153 * t25 + t154 * t24;
t102 = t149 * pkin(4);
t101 = t150 * pkin(4);
t100 = t151 * pkin(4);
t87 = t98 - t103 / 0.2e1;
t65 = t66 * pkin(2);
t56 = t67 * t79 + t77 * t82;
t55 = t75 * t87 + t109;
t52 = -t120 * t87 + t127;
t49 = t53 * pkin(2);
t48 = t50 * pkin(2);
t34 = t82 * t128 - t54 * t79;
t32 = -t82 * t110 - t51 * t79;
t1 = [(-m(2) - m(3) - m(5) + t155) * g(3) (-m(5) * t65 + t155 * (-t67 * pkin(8) + t65) - t147 * t67 - t146 * t66) * g(3) + (m(5) * t48 + t155 * (t52 * pkin(8) - t48) + t147 * t52 + t146 * t50) * g(2) + (m(5) * t49 + t155 * (t55 * pkin(8) - t49) + t147 * t55 + t146 * t53) * g(1) (t145 * (t56 * t71 - t57 * t84) + t142 * t57 + t143 * t56) * g(3) + (t145 * (t32 * t71 - t33 * t84) + t142 * t33 + t143 * t32) * g(2) + (t145 * (t34 * t71 - t35 * t84) + t142 * t35 + t143 * t34) * g(1) (-t151 * mrSges(5,1) - (-t57 * t81 + t66 * t78) * mrSges(5,2) - m(6) * t100 - m(7) * (t100 + t107) + t106) * g(3) + (-t149 * mrSges(5,1) - (-t33 * t81 - t50 * t78) * mrSges(5,2) - m(6) * t102 - m(7) * (t102 + t111) + t117) * g(2) + (-t150 * mrSges(5,1) - (-t35 * t81 - t53 * t78) * mrSges(5,2) - m(6) * t101 - m(7) * (t101 + t108) + t116) * g(1) (-m(7) * t107 + t106) * g(3) + (-m(7) * t111 + t117) * g(2) + (-m(7) * t108 + t116) * g(1) (-g(1) * t13 - g(2) * t11 - g(3) * t24) * m(7)];
taug  = t1(:);
