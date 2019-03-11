% Calculate Gravitation load on the joints for
% S6RRRRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-03-10 04:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRRR7_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR7_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR7_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR7_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR7_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR7_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 04:29:02
% EndTime: 2019-03-10 04:29:06
% DurationCPUTime: 1.68s
% Computational Cost: add. (967->161), mult. (1809->207), div. (0->0), fcn. (2146->14), ass. (0->77)
t57 = qJ(4) + qJ(5);
t52 = cos(t57);
t62 = cos(qJ(4));
t54 = t62 * pkin(4);
t37 = pkin(5) * t52 + t54;
t35 = pkin(3) + t37;
t53 = qJ(6) + t57;
t48 = sin(t53);
t49 = cos(t53);
t50 = t54 + pkin(3);
t51 = sin(t57);
t59 = sin(qJ(4));
t126 = m(5) * pkin(3) + m(6) * t50 + m(7) * t35 + t62 * mrSges(5,1) + t52 * mrSges(6,1) + t49 * mrSges(7,1) - t59 * mrSges(5,2) - t51 * mrSges(6,2) - t48 * mrSges(7,2) + mrSges(4,1);
t65 = -pkin(11) - pkin(10);
t146 = mrSges(4,2) + m(7) * (-pkin(12) + t65) - mrSges(7,3) + m(6) * t65 - mrSges(6,3) - m(5) * pkin(10) - mrSges(5,3);
t60 = sin(qJ(3));
t63 = cos(qJ(3));
t152 = -t126 * t63 + t146 * t60 - mrSges(3,1);
t150 = -m(4) - m(5);
t114 = pkin(4) * t59;
t36 = pkin(5) * t51 + t114;
t148 = m(7) * (pkin(9) + t36) + m(6) * (pkin(9) + t114);
t127 = m(6) + m(7) - t150;
t145 = pkin(2) * t127 - t152;
t110 = sin(qJ(1));
t61 = sin(qJ(2));
t64 = cos(qJ(2));
t111 = cos(qJ(1));
t97 = cos(pkin(6));
t87 = t97 * t111;
t32 = t110 * t64 + t61 * t87;
t58 = sin(pkin(6));
t93 = t58 * t111;
t20 = t32 * t63 - t60 * t93;
t31 = t110 * t61 - t64 * t87;
t144 = t20 * t59 - t31 * t62;
t143 = t20 * t62 + t31 * t59;
t142 = t20 * t51 - t31 * t52;
t141 = t20 * t52 + t31 * t51;
t140 = t20 * t48 - t31 * t49;
t139 = -t20 * t49 - t31 * t48;
t134 = mrSges(4,3) - mrSges(3,2);
t138 = -t59 * mrSges(5,1) - t51 * mrSges(6,1) - t48 * mrSges(7,1) - t62 * mrSges(5,2) - t52 * mrSges(6,2) - t49 * mrSges(7,2) - t134;
t133 = -m(6) * pkin(4) - mrSges(5,1);
t102 = t58 * t64;
t103 = t58 * t61;
t30 = t103 * t63 + t60 * t97;
t101 = (-t102 * t49 - t30 * t48) * mrSges(7,1) + (t102 * t48 - t30 * t49) * mrSges(7,2);
t77 = -t102 * t52 - t30 * t51;
t132 = -t77 * mrSges(6,1) - (t102 * t51 - t30 * t52) * mrSges(6,2) - t101;
t86 = t97 * t110;
t34 = t111 * t64 - t61 * t86;
t92 = t58 * t110;
t24 = t34 * t63 + t60 * t92;
t33 = t111 * t61 + t64 * t86;
t11 = -t24 * t51 + t33 * t52;
t5 = -t24 * t48 + t33 * t49;
t6 = t24 * t49 + t33 * t48;
t116 = t5 * mrSges(7,1) - t6 * mrSges(7,2);
t12 = t24 * t52 + t33 * t51;
t131 = -t11 * mrSges(6,1) + t12 * mrSges(6,2) - t116;
t117 = -t140 * mrSges(7,1) + t139 * mrSges(7,2);
t130 = mrSges(6,1) * t142 + t141 * mrSges(6,2) - t117;
t122 = t134 + t148;
t121 = t150 * pkin(9) + t138 - t148;
t118 = m(7) * pkin(5);
t98 = t111 * pkin(1) + pkin(8) * t92;
t96 = t34 * pkin(2) + t98;
t89 = -t32 * t60 - t63 * t93;
t88 = -pkin(1) * t110 + pkin(8) * t93;
t13 = -t24 * t59 + t33 * t62;
t79 = t33 * pkin(9) + t96;
t78 = -t32 * pkin(2) + t88;
t71 = -t31 * pkin(9) + t78;
t23 = t34 * t60 - t63 * t92;
t14 = t24 * t62 + t33 * t59;
t1 = [(-t111 * mrSges(2,1) + t110 * mrSges(2,2) - m(3) * t98 - t34 * mrSges(3,1) - mrSges(3,3) * t92 - m(4) * t79 - t24 * mrSges(4,1) - m(5) * (pkin(3) * t24 + t79) - t14 * mrSges(5,1) - t13 * mrSges(5,2) - m(6) * (t24 * t50 + t96) - t12 * mrSges(6,1) - t11 * mrSges(6,2) - m(7) * (t24 * t35 + t96) - t6 * mrSges(7,1) - t5 * mrSges(7,2) - t122 * t33 + t146 * t23) * g(2) + (t110 * mrSges(2,1) + t111 * mrSges(2,2) - m(3) * t88 + t32 * mrSges(3,1) - mrSges(3,3) * t93 - m(4) * t71 + t20 * mrSges(4,1) - m(5) * (-pkin(3) * t20 + t71) + t143 * mrSges(5,1) - t144 * mrSges(5,2) - m(6) * (-t20 * t50 + t78) + t141 * mrSges(6,1) - t142 * mrSges(6,2) - m(7) * (-t20 * t35 + t78) - t139 * mrSges(7,1) - t140 * mrSges(7,2) + t122 * t31 + t146 * t89) * g(1) (-t127 * (pkin(2) * t102 + pkin(9) * t103) + (t152 * t64 + (-m(6) * t114 - m(7) * t36 + t138) * t61) * t58) * g(3) + (t121 * t32 + t145 * t31) * g(2) + (t121 * t34 + t145 * t33) * g(1) (t146 * t30 - t126 * (-t103 * t60 + t63 * t97)) * g(3) + (-t126 * t89 + t146 * t20) * g(2) + (t126 * t23 + t146 * t24) * g(1) (-(t102 * t59 - t30 * t62) * mrSges(5,2) - m(7) * (-t102 * t37 - t30 * t36) + t133 * (-t102 * t62 - t30 * t59) + t132) * g(3) + (t143 * mrSges(5,2) - m(7) * (-t20 * t36 + t31 * t37) - t133 * t144 + t130) * g(2) + (t14 * mrSges(5,2) - m(7) * (-t24 * t36 + t33 * t37) + t133 * t13 + t131) * g(1) (-t118 * t77 + t132) * g(3) + (t118 * t142 + t130) * g(2) + (-t11 * t118 + t131) * g(1), -g(1) * t116 - g(2) * t117 - g(3) * t101];
taug  = t1(:);
