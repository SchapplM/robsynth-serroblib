% Calculate Gravitation load on the joints for
% S6RPRRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2018-11-23 16:40
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RPRRRR11_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR11_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR11_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR11_gravloadJ_floatb_twist_slag_vp2: pkin has to be [13x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR11_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR11_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:39:46
% EndTime: 2018-11-23 16:39:48
% DurationCPUTime: 1.48s
% Computational Cost: add. (4120->151), mult. (4349->199), div. (0->0), fcn. (4343->24), ass. (0->92)
t133 = pkin(7) - qJ(3);
t126 = cos(t133);
t118 = t126 / 0.2e1;
t132 = pkin(7) + qJ(3);
t125 = cos(t132);
t109 = t118 - t125 / 0.2e1;
t135 = sin(pkin(6));
t102 = t109 * t135;
t143 = cos(qJ(1));
t116 = sin(t132) / 0.2e1;
t124 = sin(t133);
t156 = t116 - t124 / 0.2e1;
t142 = sin(qJ(1));
t130 = pkin(6) + pkin(13);
t110 = sin(t130) / 0.2e1;
t131 = pkin(6) - pkin(13);
t122 = sin(t131);
t54 = t110 - t122 / 0.2e1;
t69 = cos(pkin(13));
t48 = t142 * t69 + t143 * t54;
t75 = cos(qJ(3));
t134 = sin(pkin(13));
t111 = cos(t131) / 0.2e1;
t123 = cos(t130);
t96 = t111 + t123 / 0.2e1;
t89 = t142 * t134 - t143 * t96;
t161 = t156 * t89 - t48 * t75;
t27 = t102 * t143 + t161;
t136 = cos(pkin(7));
t112 = t136 * t135;
t68 = sin(pkin(7));
t37 = -t143 * t112 + t89 * t68;
t71 = sin(qJ(4));
t74 = cos(qJ(4));
t12 = -t27 * t74 + t37 * t71;
t169 = t27 * t71 + t37 * t74;
t160 = -m(5) - m(6);
t70 = sin(qJ(5));
t166 = mrSges(4,2) - mrSges(5,3) - m(7) * (t70 * pkin(5) + pkin(10));
t73 = cos(qJ(5));
t63 = pkin(5) * t73 + pkin(4);
t67 = qJ(5) + qJ(6);
t64 = sin(t67);
t65 = cos(t67);
t154 = m(6) * pkin(4) + m(7) * t63 + t73 * mrSges(6,1) + t65 * mrSges(7,1) - t70 * mrSges(6,2) - t64 * mrSges(7,2) + mrSges(5,1);
t148 = mrSges(5,2) + m(7) * (-pkin(12) - pkin(11)) - mrSges(7,3) - m(6) * pkin(11) - mrSges(6,3);
t86 = t134 * t143 + t142 * t96;
t77 = t142 * t112 + t86 * t68;
t155 = m(7) - t160;
t165 = pkin(3) * t155 - t148 * t71 + t154 * t74 + mrSges(4,1);
t164 = -t70 * mrSges(6,1) - t64 * mrSges(7,1) - t73 * mrSges(6,2) - t65 * mrSges(7,2) + t166;
t55 = t111 - t123 / 0.2e1;
t95 = t110 + t122 / 0.2e1;
t163 = t156 * t95 + t55 * t75;
t49 = -t142 * t54 + t143 * t69;
t162 = -t156 * t86 + t49 * t75;
t157 = -m(7) * pkin(5) - mrSges(6,1);
t117 = t125 / 0.2e1;
t101 = t118 + t117;
t72 = sin(qJ(3));
t98 = t116 + t124 / 0.2e1;
t93 = t98 * t135;
t23 = t101 * t89 + t143 * t93 + t48 * t72;
t149 = t160 * pkin(10) + t164;
t146 = (-t12 * t64 + t23 * t65) * mrSges(7,1) + (-t12 * t65 - t23 * t64) * mrSges(7,2);
t29 = t102 * t142 + t162;
t16 = t29 * t74 + t71 * t77;
t28 = t101 * t86 - t142 * t93 + t49 * t72;
t5 = -t16 * t64 + t28 * t65;
t6 = t16 * t65 + t28 * t64;
t145 = t5 * mrSges(7,1) - t6 * mrSges(7,2);
t137 = cos(pkin(6));
t33 = t109 * t137 + t163;
t47 = t136 * t137 - t68 * t95;
t18 = t33 * t74 + t47 * t71;
t32 = -t101 * t95 - t137 * t98 + t55 * t72;
t144 = (-t18 * t64 + t32 * t65) * mrSges(7,1) + (-t18 * t65 - t32 * t64) * mrSges(7,2);
t119 = t135 * t142;
t138 = t143 * pkin(1) + qJ(2) * t119;
t120 = t143 * t135;
t121 = -pkin(1) * t142 + qJ(2) * t120;
t7 = -t16 * t70 + t28 * t73;
t100 = t117 - t126 / 0.2e1;
t94 = t100 * t135;
t84 = -t48 * pkin(2) - pkin(9) * t37 + t121;
t82 = t27 * pkin(3) + t84;
t81 = t49 * pkin(2) + pkin(9) * t77 + t138;
t80 = t29 * pkin(3) + t81;
t78 = t28 * pkin(10) + t80;
t15 = t29 * t71 - t74 * t77;
t8 = t16 * t73 + t28 * t70;
t1 = [(-t143 * mrSges(2,1) + t142 * mrSges(2,2) - m(3) * t138 - t49 * mrSges(3,1) + t86 * mrSges(3,2) - mrSges(3,3) * t119 - m(4) * t81 - t29 * mrSges(4,1) - t77 * mrSges(4,3) - m(5) * t78 - t16 * mrSges(5,1) - m(6) * (t16 * pkin(4) + t78) - t8 * mrSges(6,1) - t7 * mrSges(6,2) - m(7) * (t16 * t63 + t80) - t6 * mrSges(7,1) - t5 * mrSges(7,2) + t166 * t28 + t148 * t15) * g(2) + (-m(3) * t121 - m(4) * t84 - m(7) * t82 + t142 * mrSges(2,1) + t48 * mrSges(3,1) - t27 * mrSges(4,1) + t143 * mrSges(2,2) - t89 * mrSges(3,2) - mrSges(3,3) * t120 + t37 * mrSges(4,3) + t160 * (-pkin(10) * t23 + t82) + t154 * t12 - t164 * t23 + t148 * t169) * g(1) (-g(1) * t119 + g(2) * t120 - g(3) * t137) * (m(3) + m(4) + t155) (t149 * (-t100 * t137 + t163) + t165 * t32) * g(3) + (t149 * (t143 * t94 - t161) + t165 * t23) * g(2) + (t149 * (-t142 * t94 + t162) + t165 * t28) * g(1) (t148 * t18 - t154 * (-t33 * t71 + t47 * t74)) * g(3) + (t148 * t12 - t154 * t169) * g(2) + (t148 * t16 + t154 * t15) * g(1) (-(-t18 * t73 - t32 * t70) * mrSges(6,2) - t144 + t157 * (-t18 * t70 + t32 * t73)) * g(3) + (-(-t12 * t73 - t23 * t70) * mrSges(6,2) - t146 + t157 * (-t12 * t70 + t23 * t73)) * g(2) + (t8 * mrSges(6,2) + t157 * t7 - t145) * g(1), -g(1) * t145 - g(2) * t146 - g(3) * t144];
taug  = t1(:);
