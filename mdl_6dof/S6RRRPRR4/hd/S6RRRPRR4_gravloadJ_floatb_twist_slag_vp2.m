% Calculate Gravitation load on the joints for
% S6RRRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2018-11-23 17:53
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRRPRR4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR4_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:52:35
% EndTime: 2018-11-23 17:52:36
% DurationCPUTime: 0.92s
% Computational Cost: add. (597->137), mult. (584->146), div. (0->0), fcn. (517->12), ass. (0->74)
t53 = pkin(11) + qJ(5);
t47 = qJ(6) + t53;
t40 = sin(t47);
t45 = sin(t53);
t149 = -t45 * mrSges(6,2) - t40 * mrSges(7,2);
t148 = -mrSges(6,3) - mrSges(7,3);
t41 = cos(t47);
t46 = cos(t53);
t147 = t46 * mrSges(6,1) + t41 * mrSges(7,1);
t146 = m(5) * qJ(4);
t145 = -mrSges(5,3) + t148;
t55 = sin(pkin(11));
t105 = t55 * mrSges(5,2);
t54 = qJ(2) + qJ(3);
t48 = sin(t54);
t49 = cos(t54);
t144 = -t49 * t146 + (-t105 + t149) * t48;
t56 = cos(pkin(11));
t104 = t56 * mrSges(5,1);
t142 = t104 - t105;
t59 = sin(qJ(1));
t61 = cos(qJ(1));
t132 = g(1) * t61 + g(2) * t59;
t141 = (t104 + t147) * t48;
t140 = -t49 * mrSges(4,1) + (mrSges(4,2) + t148) * t48;
t42 = t56 * pkin(4) + pkin(3);
t57 = -pkin(9) - qJ(4);
t134 = t49 * t42 - t48 * t57;
t18 = pkin(5) * t46 + t42;
t52 = -pkin(10) + t57;
t135 = t49 * t18 - t48 * t52;
t139 = -m(6) * t134 - m(7) * t135;
t138 = -m(6) - m(7);
t136 = m(7) * pkin(5) + mrSges(6,1);
t120 = pkin(3) * t48;
t58 = sin(qJ(2));
t121 = pkin(2) * t58;
t75 = -t42 * t48 - t49 * t57;
t77 = -t18 * t48 - t49 * t52;
t131 = -m(7) * (t77 - t121) - m(6) * (t75 - t121) - m(5) * (-t120 - t121) + t141;
t62 = -pkin(8) - pkin(7);
t130 = -m(3) * pkin(7) + m(5) * t62 - t55 * mrSges(5,1) - t56 * mrSges(5,2) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t38 = t48 * mrSges(5,3);
t129 = -t38 + t140 + (-t142 - t147 - t149) * t49;
t128 = m(5) * t120 - m(6) * t75 - m(7) * t77 + t141;
t106 = t49 * t61;
t127 = t145 * t106 + t144 * t61;
t107 = t49 * t59;
t126 = t145 * t107 + t144 * t59;
t60 = cos(qJ(2));
t83 = t60 * mrSges(3,1) - t58 * mrSges(3,2);
t125 = -(m(5) * pkin(3) + t142) * t49 - mrSges(2,1) - m(3) * pkin(1) - t83 + t140;
t5 = t107 * t40 + t41 * t61;
t6 = -t107 * t41 + t40 * t61;
t123 = -t5 * mrSges(7,1) + t6 * mrSges(7,2);
t7 = -t106 * t40 + t41 * t59;
t8 = t106 * t41 + t40 * t59;
t122 = t7 * mrSges(7,1) - t8 * mrSges(7,2);
t119 = pkin(4) * t55;
t118 = pkin(5) * t45;
t115 = g(3) * t48;
t51 = t60 * pkin(2);
t35 = t48 * qJ(4);
t97 = t49 * pkin(3) + t35;
t80 = mrSges(4,1) * t48 + mrSges(4,2) * t49;
t79 = -mrSges(7,1) * t40 - mrSges(7,2) * t41;
t11 = -t106 * t45 + t46 * t59;
t9 = t107 * t45 + t46 * t61;
t44 = t51 + pkin(1);
t28 = t61 * t44;
t22 = t118 + t119;
t12 = t106 * t46 + t45 * t59;
t10 = -t107 * t46 + t45 * t61;
t1 = [(-m(5) * t28 - t12 * mrSges(6,1) - t8 * mrSges(7,1) - t11 * mrSges(6,2) - t7 * mrSges(7,2) + (-m(4) + t138) * (-t59 * t62 + t28) + (-m(6) * t119 - m(7) * t22 + t130) * t59 + (-(mrSges(5,3) + t146) * t48 + t125 + t139) * t61) * g(2) + (-t10 * mrSges(6,1) - t6 * mrSges(7,1) - t9 * mrSges(6,2) - t5 * mrSges(7,2) + (m(4) * t62 - m(6) * (-t62 + t119) - m(7) * (t22 - t62) + t130) * t61 + (m(4) * t44 - m(5) * (-t44 - t35) + t38 - m(6) * (-t44 - t134) - m(7) * (-t44 - t135) - t125) * t59) * g(1) (t131 * t59 + t126) * g(2) + (t131 * t61 + t127) * g(1) + (-t83 - m(4) * t51 - m(5) * (t51 + t97) - m(6) * (t51 + t134) - m(7) * (t51 + t135) + t129) * g(3) + t132 * (m(4) * t121 + mrSges(3,1) * t58 + mrSges(3,2) * t60 + t80) t132 * t80 + (t128 * t59 + t126) * g(2) + (t128 * t61 + t127) * g(1) + (-m(5) * t97 + t129 + t139) * g(3) (g(3) * t49 - t132 * t48) * (m(5) - t138) (m(7) * t118 + mrSges(6,1) * t45 + mrSges(6,2) * t46 - t79) * t115 + (-t10 * mrSges(6,2) + t136 * t9 - t123) * g(2) + (t12 * mrSges(6,2) - t136 * t11 - t122) * g(1), -g(1) * t122 - g(2) * t123 - t115 * t79];
taug  = t1(:);
