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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 02:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

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
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP9_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP9_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 02:00:22
% EndTime: 2019-03-10 02:00:27
% DurationCPUTime: 1.62s
% Computational Cost: add. (875->156), mult. (1765->208), div. (0->0), fcn. (2087->12), ass. (0->73)
t64 = qJ(4) + qJ(5);
t60 = cos(t64);
t69 = cos(qJ(4));
t61 = t69 * pkin(4);
t47 = pkin(5) * t60 + t61;
t45 = pkin(3) + t47;
t58 = t61 + pkin(3);
t66 = sin(qJ(4));
t151 = -m(5) * pkin(3) - m(6) * t58 - m(7) * t45 - t69 * mrSges(5,1) + t66 * mrSges(5,2) - mrSges(4,1);
t72 = -pkin(11) - pkin(10);
t144 = m(6) * t72 + m(7) * (-qJ(6) + t72) + mrSges(4,2) - mrSges(6,3) - mrSges(7,3) - m(5) * pkin(10) - mrSges(5,3);
t150 = mrSges(6,1) + mrSges(7,1);
t149 = mrSges(6,2) + mrSges(7,2);
t120 = pkin(4) * t66;
t139 = mrSges(4,3) - mrSges(3,2);
t59 = sin(t64);
t46 = pkin(5) * t59 + t120;
t148 = -m(6) * t120 - m(7) * t46 - t66 * mrSges(5,1) - t69 * mrSges(5,2) - t139;
t117 = sin(qJ(1));
t68 = sin(qJ(2));
t71 = cos(qJ(2));
t101 = cos(pkin(6));
t118 = cos(qJ(1));
t90 = t101 * t118;
t42 = t117 * t71 + t68 * t90;
t67 = sin(qJ(3));
t70 = cos(qJ(3));
t65 = sin(pkin(6));
t97 = t65 * t118;
t28 = t42 * t70 - t67 * t97;
t41 = t117 * t68 - t71 * t90;
t147 = t28 * t66 - t41 * t69;
t146 = t28 * t69 + t41 * t66;
t145 = t28 * t59 - t41 * t60;
t10 = -t28 * t60 - t41 * t59;
t143 = t144 * t67 + t151 * t70 - mrSges(3,1);
t138 = -m(6) * pkin(4) - mrSges(5,1);
t109 = t65 * t71;
t110 = t65 * t68;
t40 = t101 * t67 + t110 * t70;
t25 = -t109 * t60 - t40 * t59;
t137 = -t149 * (t109 * t59 - t40 * t60) - t150 * t25;
t89 = t101 * t117;
t44 = t118 * t71 - t68 * t89;
t96 = t65 * t117;
t32 = t44 * t70 + t67 * t96;
t43 = t118 * t68 + t71 * t89;
t13 = -t32 * t59 + t43 * t60;
t14 = t32 * t60 + t43 * t59;
t136 = -t150 * t13 + t149 * t14;
t135 = -t149 * t10 + t145 * t150;
t134 = -m(4) - m(6) - m(7);
t133 = -t149 * t59 + t150 * t60 - t151;
t130 = m(6) * (pkin(9) + t120) + m(7) * (pkin(9) + t46) + t139;
t128 = -m(5) * pkin(9) + t148;
t125 = m(7) * pkin(5);
t112 = t59 * t70;
t111 = t60 * t70;
t106 = t70 * t71;
t102 = t118 * pkin(1) + pkin(8) * t96;
t100 = t44 * pkin(2) + t102;
t91 = -pkin(1) * t117 + pkin(8) * t97;
t15 = -t32 * t66 + t43 * t69;
t83 = t43 * pkin(9) + t100;
t82 = -t42 * pkin(2) + t91;
t27 = t42 * t67 + t70 * t97;
t76 = -t41 * pkin(9) + t82;
t39 = -t101 * t70 + t110 * t67;
t37 = t43 * pkin(2);
t35 = t41 * pkin(2);
t31 = t44 * t67 - t70 * t96;
t16 = t32 * t69 + t43 * t66;
t1 = [(-t118 * mrSges(2,1) + t117 * mrSges(2,2) - m(3) * t102 - t44 * mrSges(3,1) - mrSges(3,3) * t96 - m(4) * t83 - t32 * mrSges(4,1) - m(5) * (pkin(3) * t32 + t83) - t16 * mrSges(5,1) - t15 * mrSges(5,2) - m(6) * (t32 * t58 + t100) - m(7) * (t32 * t45 + t100) - t130 * t43 - t150 * t14 - t149 * t13 + t144 * t31) * g(2) + (t117 * mrSges(2,1) + t118 * mrSges(2,2) - m(3) * t91 + t42 * mrSges(3,1) - mrSges(3,3) * t97 - m(4) * t76 + t28 * mrSges(4,1) - m(5) * (-pkin(3) * t28 + t76) + t146 * mrSges(5,1) - t147 * mrSges(5,2) - m(6) * (-t28 * t58 + t82) - m(7) * (-t28 * t45 + t82) + t130 * t41 - t150 * t10 - t149 * t145 - t144 * t27) * g(1) (m(5) * t35 - t150 * (-t111 * t41 + t42 * t59) - t149 * (t112 * t41 + t42 * t60) + t134 * (t42 * pkin(9) - t35) + t128 * t42 - t143 * t41) * g(2) + (m(5) * t37 - t150 * (-t111 * t43 + t44 * t59) - t149 * (t112 * t43 + t44 * t60) + t134 * (t44 * pkin(9) - t37) + t128 * t44 - t143 * t43) * g(1) + ((-m(5) + t134) * (pkin(2) * t109 + pkin(9) * t110) + (-t150 * (t106 * t60 + t59 * t68) - t149 * (-t106 * t59 + t60 * t68) + t143 * t71 + t148 * t68) * t65) * g(3) (t133 * t39 + t144 * t40) * g(3) + (t133 * t27 + t144 * t28) * g(2) + (t133 * t31 + t144 * t32) * g(1) (-(t109 * t66 - t40 * t69) * mrSges(5,2) - m(7) * (-t109 * t47 - t40 * t46) + t138 * (-t109 * t69 - t40 * t66) + t137) * g(3) + (t146 * mrSges(5,2) - m(7) * (-t28 * t46 + t41 * t47) - t138 * t147 + t135) * g(2) + (t16 * mrSges(5,2) - m(7) * (-t32 * t46 + t43 * t47) + t138 * t15 + t136) * g(1) (-t125 * t25 + t137) * g(3) + (t125 * t145 + t135) * g(2) + (-t125 * t13 + t136) * g(1) (-g(1) * t31 - g(2) * t27 - g(3) * t39) * m(7)];
taug  = t1(:);
