% Calculate Gravitation load on the joints for
% S6RRRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,theta5]';
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
% Datum: 2018-11-23 18:04
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRRRPP1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP1_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:04:17
% EndTime: 2018-11-23 18:04:19
% DurationCPUTime: 1.08s
% Computational Cost: add. (611->133), mult. (662->153), div. (0->0), fcn. (602->10), ass. (0->75)
t145 = -mrSges(6,1) - mrSges(7,1);
t122 = m(7) * pkin(5) - t145;
t143 = -mrSges(6,3) - mrSges(7,2);
t144 = mrSges(6,2) - mrSges(7,3);
t52 = sin(qJ(4));
t104 = t52 * mrSges(5,2);
t49 = qJ(4) + pkin(10);
t44 = sin(t49);
t50 = qJ(2) + qJ(3);
t46 = sin(t50);
t113 = t44 * t46;
t142 = -mrSges(6,2) * t113 - t46 * t104;
t134 = m(6) + m(7);
t135 = pkin(4) * t134 + mrSges(5,1);
t141 = -mrSges(5,3) + t143;
t47 = cos(t50);
t133 = t47 * pkin(3) + t46 * pkin(9);
t140 = m(5) * t133;
t138 = t143 * t46;
t137 = -t47 * mrSges(4,1) + (mrSges(4,2) - mrSges(5,3)) * t46;
t54 = sin(qJ(1));
t55 = cos(qJ(4));
t100 = t54 * t55;
t57 = cos(qJ(1));
t102 = t52 * t57;
t7 = -t47 * t102 + t100;
t130 = g(1) * t57 + g(2) * t54;
t114 = mrSges(5,1) * t55;
t42 = pkin(4) * t55 + pkin(3);
t45 = cos(t49);
t92 = qJ(6) * t44;
t136 = (t114 - m(7) * (-t42 - t92) + t44 * mrSges(7,3) + t122 * t45) * t46;
t129 = -m(3) * pkin(7) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t105 = t47 * t57;
t128 = t105 * t141 + t142 * t57;
t51 = -qJ(5) - pkin(9);
t107 = t47 * t51;
t53 = sin(qJ(2));
t120 = pkin(2) * t53;
t73 = -t42 * t46 - t107;
t127 = -m(7) * (-t107 - t120) - m(6) * (t73 - t120) - m(5) * (-pkin(3) * t46 - t120) + t136;
t126 = -m(6) * t73 + m(7) * t107 + t136;
t106 = t47 * t54;
t125 = t106 * t141 + t142 * t54;
t56 = cos(qJ(2));
t77 = t56 * mrSges(3,1) - t53 * mrSges(3,2);
t124 = m(3) * pkin(1) + mrSges(2,1) - t137 + t77;
t111 = t45 * t47;
t123 = t145 * t111 + t137 + t138 + (t144 * t44 + t104 - t114) * t47;
t121 = m(7) * qJ(6) - t144;
t48 = t56 * pkin(2);
t109 = t46 * t54;
t108 = t46 * t57;
t19 = t47 * t42;
t103 = t52 * t54;
t101 = t54 * t45;
t99 = t55 * t57;
t98 = t57 * t44;
t58 = -pkin(8) - pkin(7);
t97 = t57 * t58;
t83 = -t46 * t51 + t19;
t43 = t48 + pkin(1);
t82 = t57 * t43 - t54 * t58;
t74 = mrSges(4,1) * t46 + mrSges(4,2) * t47;
t69 = pkin(5) * t111 + t47 * t92 + t83;
t5 = t103 * t47 + t99;
t30 = pkin(9) * t105;
t29 = pkin(9) * t106;
t8 = t47 * t99 + t103;
t6 = -t100 * t47 + t102;
t4 = t105 * t45 + t44 * t54;
t3 = t47 * t98 - t101;
t2 = t101 * t47 - t98;
t1 = t106 * t44 + t45 * t57;
t9 = [(-t8 * mrSges(5,1) - t7 * mrSges(5,2) + (-m(4) - m(5)) * t82 - t134 * (pkin(4) * t103 + t42 * t105 - t108 * t51 + t82) - t122 * t4 - t121 * t3 + t143 * t108 + t129 * t54 + (-t124 - t140) * t57) * g(2) + (m(5) * t97 - t6 * mrSges(5,1) - t5 * mrSges(5,2) - t134 * (pkin(4) * t102 + t51 * t109 - t97) + t122 * t2 + t121 * t1 + (m(4) * t58 + t129) * t57 + (m(4) * t43 - m(5) * (-t43 - t133) - t134 * (-t43 - t19) + t124 - t138) * t54) * g(1) (-m(5) * t29 + t127 * t54 + t125) * g(2) + (-m(5) * t30 + t127 * t57 + t128) * g(1) + (-t77 - m(4) * t48 - m(5) * (t48 + t133) - m(6) * (t48 + t83) - m(7) * (t48 + t69) + t123) * g(3) + t130 * (m(4) * t120 + mrSges(3,1) * t53 + mrSges(3,2) * t56 + t74) t130 * t74 + (-m(5) * (-pkin(3) * t109 + t29) + t126 * t54 + t125) * g(2) + (-m(5) * (-pkin(3) * t108 + t30) + t126 * t57 + t128) * g(1) + (-m(6) * t83 - m(7) * t69 + t123 - t140) * g(3) (mrSges(5,2) * t55 - t121 * t45 + t122 * t44 + t135 * t52) * g(3) * t46 + (-t6 * mrSges(5,2) + t122 * t1 - t121 * t2 + t135 * t5) * g(2) + (t8 * mrSges(5,2) - t121 * t4 + t122 * t3 - t135 * t7) * g(1) (g(3) * t47 - t130 * t46) * t134 (-g(1) * t3 - g(2) * t1 - g(3) * t113) * m(7)];
taug  = t9(:);
