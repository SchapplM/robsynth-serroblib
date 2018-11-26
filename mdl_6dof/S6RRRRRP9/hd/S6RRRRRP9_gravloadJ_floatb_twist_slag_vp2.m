% Calculate Gravitation load on the joints for
% S6RRRRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
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
% Datum: 2018-11-23 18:33
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRRRRP9_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP9_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP9_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP9_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP9_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP9_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:32:50
% EndTime: 2018-11-23 18:32:51
% DurationCPUTime: 1.39s
% Computational Cost: add. (1811->156), mult. (2077->202), div. (0->0), fcn. (2087->16), ass. (0->83)
t70 = qJ(4) + qJ(5);
t66 = cos(t70);
t75 = cos(qJ(4));
t67 = t75 * pkin(4);
t54 = pkin(5) * t66 + t67;
t50 = pkin(3) + t54;
t64 = t67 + pkin(3);
t71 = sin(qJ(4));
t156 = -m(5) * pkin(3) - m(6) * t64 - m(7) * t50 - t75 * mrSges(5,1) + t71 * mrSges(5,2) - mrSges(4,1);
t78 = -pkin(11) - pkin(10);
t80 = m(6) * t78 + m(7) * (-qJ(6) + t78) + mrSges(4,2) - mrSges(6,3) - mrSges(7,3) - m(5) * pkin(10) - mrSges(5,3);
t158 = mrSges(6,1) + mrSges(7,1);
t157 = mrSges(6,2) + mrSges(7,2);
t74 = sin(qJ(1));
t77 = cos(qJ(2));
t125 = t74 * t77;
t132 = cos(qJ(1));
t116 = pkin(6) - qJ(2);
t102 = sin(t116);
t115 = pkin(6) + qJ(2);
t101 = sin(t115);
t95 = t101 / 0.2e1;
t88 = t95 - t102 / 0.2e1;
t42 = t132 * t88 + t125;
t72 = sin(qJ(3));
t76 = cos(qJ(3));
t117 = sin(pkin(6));
t98 = t132 * t117;
t30 = t42 * t76 - t72 * t98;
t73 = sin(qJ(2));
t103 = cos(t116);
t97 = cos(t115) / 0.2e1;
t81 = t103 / 0.2e1 + t97;
t41 = -t132 * t81 + t73 * t74;
t65 = sin(t70);
t155 = t30 * t65 - t41 * t66;
t10 = -t30 * t66 - t41 * t65;
t154 = t156 * t76 + t80 * t72 - mrSges(3,1);
t133 = pkin(4) * t71;
t153 = m(6) * t133;
t150 = m(6) * pkin(4) + mrSges(5,1);
t118 = cos(pkin(6));
t52 = t97 - t103 / 0.2e1;
t40 = t118 * t72 - t52 * t76;
t96 = t102 / 0.2e1;
t51 = t95 + t96;
t25 = -t40 * t65 - t51 * t66;
t149 = -t157 * (-t40 * t66 + t51 * t65) - t158 * t25;
t104 = t74 * t117;
t111 = t132 * t77;
t45 = -t74 * t88 + t111;
t34 = t104 * t72 + t45 * t76;
t44 = t132 * t73 + t74 * t81;
t13 = -t34 * t65 + t44 * t66;
t14 = t34 * t66 + t44 * t65;
t148 = -t158 * t13 + t157 * t14;
t147 = -t157 * t10 + t155 * t158;
t146 = -m(4) - m(6) - m(7);
t113 = m(5) - t146;
t53 = pkin(5) * t65 + t133;
t100 = -m(7) * t53 + mrSges(3,2) - mrSges(4,3);
t145 = t113 * pkin(9) - t100;
t144 = -t157 * t65 + t158 * t66 - t156;
t124 = t75 * mrSges(5,2);
t142 = m(5) * pkin(9) + t71 * mrSges(5,1) - t100 + t124 + t153;
t140 = m(7) * pkin(5);
t129 = t65 * t76;
t128 = t66 * t76;
t119 = t132 * pkin(1) + pkin(8) * t104;
t114 = t45 * pkin(2) + t119;
t108 = -t74 * pkin(1) + pkin(8) * t98;
t15 = -t34 * t71 + t44 * t75;
t89 = t96 - t101 / 0.2e1;
t29 = t42 * t72 + t76 * t98;
t49 = t51 * pkin(2);
t46 = t74 * t89 + t111;
t43 = -t132 * t89 + t125;
t39 = -t118 * t76 - t52 * t72;
t37 = t44 * pkin(2);
t35 = t41 * pkin(2);
t33 = -t104 * t76 + t45 * t72;
t16 = t34 * t75 + t44 * t71;
t1 = [(-t132 * mrSges(2,1) - m(3) * t119 - t45 * mrSges(3,1) - m(4) * t114 - t34 * mrSges(4,1) - m(5) * (pkin(3) * t34 + t114) - t16 * mrSges(5,1) - t15 * mrSges(5,2) - m(6) * (t34 * t64 + t114) - m(7) * (t34 * t50 + t114) + (-mrSges(3,3) * t117 + mrSges(2,2)) * t74 - t158 * t14 - t157 * t13 + (-t145 - t153) * t44 + t80 * t33) * g(2) + (t74 * mrSges(2,1) + t132 * mrSges(2,2) - m(3) * t108 + t42 * mrSges(3,1) - mrSges(3,3) * t98 - t158 * t10 - t157 * t155 - t156 * t30 + (t150 * t71 + t124 + t145) * t41 - t80 * t29 + t113 * (t42 * pkin(2) - t108)) * g(1) (-m(5) * t49 - t158 * (t128 * t51 - t52 * t65) - t157 * (-t129 * t51 - t52 * t66) + t146 * (-t52 * pkin(9) + t49) + t142 * t52 + t154 * t51) * g(3) + (m(5) * t35 - t158 * (-t128 * t41 + t43 * t65) - t157 * (t129 * t41 + t43 * t66) + t146 * (t43 * pkin(9) - t35) - t142 * t43 - t154 * t41) * g(2) + (m(5) * t37 - t158 * (-t128 * t44 + t46 * t65) - t157 * (t129 * t44 + t46 * t66) + t146 * (t46 * pkin(9) - t37) - t142 * t46 - t154 * t44) * g(1) (t144 * t39 + t80 * t40) * g(3) + (t144 * t29 + t80 * t30) * g(2) + (t144 * t33 + t80 * t34) * g(1) (-(-t40 * t75 + t51 * t71) * mrSges(5,2) - m(7) * (-t40 * t53 - t51 * t54) - t150 * (-t40 * t71 - t51 * t75) + t149) * g(3) + (-(-t30 * t75 - t41 * t71) * mrSges(5,2) - m(7) * (-t30 * t53 + t41 * t54) - t150 * (-t30 * t71 + t41 * t75) + t147) * g(2) + (t16 * mrSges(5,2) - m(7) * (-t34 * t53 + t44 * t54) - t150 * t15 + t148) * g(1) (-t25 * t140 + t149) * g(3) + (t140 * t155 + t147) * g(2) + (-t13 * t140 + t148) * g(1) (-g(1) * t33 - g(2) * t29 - g(3) * t39) * m(7)];
taug  = t1(:);
