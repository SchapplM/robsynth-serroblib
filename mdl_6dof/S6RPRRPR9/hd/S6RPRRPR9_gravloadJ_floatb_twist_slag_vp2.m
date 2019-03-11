% Calculate Gravitation load on the joints for
% S6RPRRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
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
% Datum: 2019-03-09 05:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRPR9_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR9_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR9_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR9_gravloadJ_floatb_twist_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR9_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR9_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:26:56
% EndTime: 2019-03-09 05:27:00
% DurationCPUTime: 1.35s
% Computational Cost: add. (1086->140), mult. (2623->201), div. (0->0), fcn. (3289->16), ass. (0->68)
t110 = sin(pkin(7));
t63 = sin(pkin(6));
t104 = t63 * t110;
t112 = cos(pkin(7));
t70 = cos(qJ(1));
t109 = sin(pkin(12));
t127 = sin(qJ(1));
t111 = cos(pkin(12));
t113 = cos(pkin(6));
t92 = t113 * t111;
t83 = t127 * t109 - t70 * t92;
t146 = t70 * t104 + t83 * t112;
t105 = t63 * t112;
t37 = t70 * t105 - t83 * t110;
t66 = sin(qJ(4));
t123 = t37 * t66;
t128 = cos(qJ(3));
t90 = t113 * t109;
t47 = t111 * t127 + t70 * t90;
t67 = sin(qJ(3));
t24 = -t128 * t47 + t146 * t67;
t69 = cos(qJ(4));
t145 = t24 * t69 + t123;
t138 = -t24 * t66 + t37 * t69;
t62 = qJ(4) + pkin(13);
t59 = sin(t62);
t60 = cos(t62);
t144 = t24 * t59 - t37 * t60;
t143 = t24 * t60 + t37 * t59;
t65 = sin(qJ(6));
t68 = cos(qJ(6));
t136 = m(7) * pkin(5) + t68 * mrSges(7,1) - t65 * mrSges(7,2) + mrSges(6,1);
t133 = -m(7) * pkin(11) + mrSges(6,2) - mrSges(7,3);
t78 = t109 * t70 + t127 * t92;
t39 = t127 * t105 + t78 * t110;
t89 = t112 * t111;
t91 = t113 * t110;
t36 = t67 * t91 + (t109 * t128 + t67 * t89) * t63;
t46 = -t104 * t111 + t112 * t113;
t139 = -t36 * t66 + t46 * t69;
t135 = -t127 * t104 + t78 * t112;
t48 = t111 * t70 - t127 * t90;
t26 = t48 * t128 - t135 * t67;
t9 = -t26 * t66 + t39 * t69;
t21 = t146 * t128 + t47 * t67;
t137 = m(6) + m(7);
t93 = -m(5) * pkin(10) + mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t131 = -t65 * mrSges(7,1) - t68 * mrSges(7,2) + t93;
t132 = m(5) * pkin(3) + t69 * mrSges(5,1) - t66 * mrSges(5,2) - t133 * t59 + t136 * t60 + mrSges(4,1);
t121 = t39 * t66;
t118 = t63 * t70;
t106 = t63 * t127;
t114 = t70 * pkin(1) + qJ(2) * t106;
t97 = -pkin(1) * t127 + qJ(2) * t118;
t74 = -t47 * pkin(2) + t37 * pkin(9) + t97;
t73 = t48 * pkin(2) + t39 * pkin(9) + t114;
t25 = t135 * t128 + t48 * t67;
t58 = pkin(4) * t69 + pkin(3);
t64 = -qJ(5) - pkin(10);
t71 = pkin(4) * t121 - t25 * t64 + t26 * t58 + t73;
t35 = t63 * t109 * t67 + (-t63 * t89 - t91) * t128;
t16 = t36 * t60 + t46 * t59;
t10 = t26 * t69 + t121;
t8 = t26 * t60 + t39 * t59;
t7 = t26 * t59 - t39 * t60;
t2 = t25 * t65 + t68 * t8;
t1 = t25 * t68 - t65 * t8;
t3 = [(-t70 * mrSges(2,1) + t127 * mrSges(2,2) - m(3) * t114 - t48 * mrSges(3,1) + t78 * mrSges(3,2) - mrSges(3,3) * t106 - m(4) * t73 - t26 * mrSges(4,1) - t39 * mrSges(4,3) - m(5) * (t26 * pkin(3) + t73) - t10 * mrSges(5,1) - t9 * mrSges(5,2) - m(6) * t71 - t8 * mrSges(6,1) - m(7) * (t8 * pkin(5) + t71) - t2 * mrSges(7,1) - t1 * mrSges(7,2) + t133 * t7 + t93 * t25) * g(2) + (t127 * mrSges(2,1) + t70 * mrSges(2,2) - m(3) * t97 + t47 * mrSges(3,1) - t83 * mrSges(3,2) - mrSges(3,3) * t118 - m(4) * t74 - t24 * mrSges(4,1) - t37 * mrSges(4,3) - m(5) * (t24 * pkin(3) + t74) - t145 * mrSges(5,1) - t138 * mrSges(5,2) + t133 * t144 - t136 * t143 - t131 * t21 - t137 * (pkin(4) * t123 + t21 * t64 + t24 * t58 + t74)) * g(1) (-t113 * g(3) + (-t127 * g(1) + g(2) * t70) * t63) * (m(3) + m(4) + m(5) + t137) (-t137 * (-t35 * t58 - t36 * t64) + t131 * t36 + t132 * t35) * g(3) + (-t137 * (-t21 * t58 + t24 * t64) - t131 * t24 + t132 * t21) * g(2) + (-t137 * (-t25 * t58 - t26 * t64) + t131 * t26 + t132 * t25) * g(1) (-t139 * mrSges(5,1) - (-t36 * t69 - t46 * t66) * mrSges(5,2) + t133 * t16 - t136 * (-t36 * t59 + t46 * t60)) * g(3) + (t138 * mrSges(5,1) - t145 * mrSges(5,2) - t133 * t143 - t136 * t144) * g(2) + (-t9 * mrSges(5,1) + t10 * mrSges(5,2) + t133 * t8 + t136 * t7) * g(1) + (-g(1) * t9 + g(2) * t138 - g(3) * t139) * t137 * pkin(4), t137 * (-g(1) * t25 - g(2) * t21 - g(3) * t35) -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * ((t143 * t65 + t21 * t68) * mrSges(7,1) + (t143 * t68 - t21 * t65) * mrSges(7,2)) - g(3) * ((-t16 * t65 + t35 * t68) * mrSges(7,1) + (-t16 * t68 - t35 * t65) * mrSges(7,2))];
taug  = t3(:);
