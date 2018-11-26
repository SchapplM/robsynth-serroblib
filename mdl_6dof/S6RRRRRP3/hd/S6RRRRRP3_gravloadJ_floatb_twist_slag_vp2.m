% Calculate Gravitation load on the joints for
% S6RRRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
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
% Datum: 2018-11-23 18:28
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRRRRP3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP3_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:27:45
% EndTime: 2018-11-23 18:27:46
% DurationCPUTime: 1.05s
% Computational Cost: add. (613->135), mult. (649->156), div. (0->0), fcn. (584->10), ass. (0->78)
t160 = mrSges(6,1) + mrSges(7,1);
t56 = qJ(4) + qJ(5);
t51 = cos(t56);
t61 = cos(qJ(4));
t161 = t61 * mrSges(5,1) + t160 * t51;
t159 = mrSges(6,2) + mrSges(7,2);
t49 = sin(t56);
t157 = t49 * t159;
t158 = -mrSges(5,3) - mrSges(6,3) - mrSges(7,3);
t60 = sin(qJ(1));
t63 = cos(qJ(1));
t144 = g(1) * t63 + g(2) * t60;
t57 = qJ(2) + qJ(3);
t50 = sin(t57);
t156 = t161 * t50;
t52 = cos(t57);
t155 = -t52 * mrSges(4,1) + (mrSges(4,2) + t158) * t50;
t53 = t61 * pkin(4);
t47 = t53 + pkin(3);
t64 = -pkin(10) - pkin(9);
t146 = t52 * t47 - t50 * t64;
t26 = pkin(5) * t51 + t53;
t22 = pkin(3) + t26;
t55 = -qJ(6) + t64;
t147 = t52 * t22 - t50 * t55;
t151 = t52 * pkin(3) + t50 * pkin(9);
t154 = -m(5) * t151 - m(6) * t146 - m(7) * t147;
t132 = m(6) * pkin(4);
t149 = mrSges(5,1) + t132;
t148 = mrSges(6,1) * t49 + t159 * t51;
t111 = t52 * t63;
t11 = -t49 * t111 + t51 * t60;
t12 = t51 * t111 + t49 * t60;
t143 = -t160 * t11 + t159 * t12;
t112 = t52 * t60;
t10 = -t51 * t112 + t49 * t63;
t9 = t49 * t112 + t51 * t63;
t142 = -t159 * t10 + t160 * t9;
t141 = -m(3) * pkin(7) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t59 = sin(qJ(2));
t126 = pkin(2) * t59;
t76 = -t47 * t50 - t52 * t64;
t78 = -t22 * t50 - t52 * t55;
t140 = -m(7) * (t78 - t126) - m(6) * (t76 - t126) - m(5) * (-pkin(3) * t50 - t126) + t156;
t139 = -m(6) * t76 - m(7) * t78 + t156;
t115 = t50 * t63;
t58 = sin(qJ(4));
t110 = t58 * mrSges(5,2);
t96 = t50 * t110;
t138 = t158 * t111 - t115 * t157 - t63 * t96;
t116 = t50 * t60;
t137 = t158 * t112 - t116 * t157 - t60 * t96;
t136 = m(4) + m(5) + m(6) + m(7);
t135 = t155 + (t110 + t157 - t161) * t52;
t125 = pkin(4) * t58;
t25 = pkin(5) * t49 + t125;
t134 = -m(6) * t125 - m(7) * t25;
t62 = cos(qJ(2));
t84 = t62 * mrSges(3,1) - t59 * mrSges(3,2);
t133 = m(3) * pkin(1) + mrSges(2,1) - t155 + t84;
t131 = m(7) * pkin(5);
t122 = g(3) * t50;
t54 = t62 * pkin(2);
t109 = t58 * t60;
t108 = t58 * t63;
t107 = t60 * t25;
t106 = t60 * t61;
t104 = t61 * t63;
t81 = mrSges(4,1) * t50 + mrSges(4,2) * t52;
t15 = -t52 * t108 + t106;
t13 = t52 * t109 + t104;
t65 = -pkin(8) - pkin(7);
t48 = t54 + pkin(1);
t39 = pkin(9) * t111;
t38 = pkin(9) * t112;
t16 = t52 * t104 + t109;
t14 = -t52 * t106 + t108;
t1 = [(-t109 * t132 - m(7) * t107 - t16 * mrSges(5,1) - t15 * mrSges(5,2) - t136 * (t63 * t48 - t60 * t65) - t160 * t12 - t159 * t11 + t141 * t60 + (-t133 + t154) * t63) * g(2) + (-t14 * mrSges(5,1) - t13 * mrSges(5,2) - t159 * t9 - t160 * t10 + (t136 * t65 + t134 + t141) * t63 + (m(4) * t48 - m(5) * (-t48 - t151) - m(6) * (-t48 - t146) - m(7) * (-t48 - t147) + t133) * t60) * g(1) (-m(5) * t38 + t140 * t60 + t137) * g(2) + (-m(5) * t39 + t140 * t63 + t138) * g(1) + (-t84 - m(4) * t54 - m(5) * (t54 + t151) - m(6) * (t54 + t146) - m(7) * (t54 + t147) + t135) * g(3) + t144 * (m(4) * t126 + mrSges(3,1) * t59 + mrSges(3,2) * t62 + t81) t144 * t81 + (-m(5) * (-pkin(3) * t116 + t38) + t139 * t60 + t137) * g(2) + (-m(5) * (-pkin(3) * t115 + t39) + t139 * t63 + t138) * g(1) + (t135 + t154) * g(3) (mrSges(5,1) * t58 + mrSges(7,1) * t49 + mrSges(5,2) * t61 - t134 + t148) * t122 + (-t14 * mrSges(5,2) - m(7) * (-t52 * t107 - t26 * t63) + t149 * t13 + t142) * g(2) + (t16 * mrSges(5,2) - m(7) * (-t25 * t111 + t26 * t60) - t149 * t15 + t143) * g(1) (-(-mrSges(7,1) - t131) * t49 + t148) * t122 + (t9 * t131 + t142) * g(2) + (-t11 * t131 + t143) * g(1) (g(3) * t52 - t144 * t50) * m(7)];
taug  = t1(:);
