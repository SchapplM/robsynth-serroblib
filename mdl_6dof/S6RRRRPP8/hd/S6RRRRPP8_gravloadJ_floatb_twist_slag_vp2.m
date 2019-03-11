% Calculate Gravitation load on the joints for
% S6RRRRPP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4]';
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
% Datum: 2019-03-09 21:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRPP8_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP8_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP8_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP8_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP8_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP8_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:30:20
% EndTime: 2019-03-09 21:30:24
% DurationCPUTime: 1.29s
% Computational Cost: add. (832->151), mult. (2072->203), div. (0->0), fcn. (2509->10), ass. (0->83)
t74 = sin(qJ(4));
t77 = cos(qJ(4));
t158 = pkin(4) * t77 + qJ(5) * t74;
t155 = mrSges(7,3) - mrSges(6,2) - mrSges(5,3);
t157 = m(7) * qJ(6) + t155;
t144 = mrSges(5,2) - mrSges(7,2) - mrSges(6,3);
t94 = m(7) * pkin(5) + mrSges(5,1) + mrSges(6,1) + mrSges(7,1);
t133 = cos(qJ(3));
t75 = sin(qJ(3));
t154 = -t133 * pkin(3) - pkin(10) * t75;
t134 = cos(qJ(1));
t73 = sin(pkin(6));
t114 = t73 * t134;
t132 = sin(qJ(1));
t76 = sin(qJ(2));
t78 = cos(qJ(2));
t119 = cos(pkin(6));
t98 = t119 * t134;
t56 = t132 * t78 + t76 * t98;
t29 = -t75 * t114 + t133 * t56;
t55 = t132 * t76 - t78 * t98;
t5 = t29 * t74 - t55 * t77;
t153 = t29 * t77 + t55 * t74;
t152 = m(6) + m(7);
t151 = -mrSges(4,1) * t133 + mrSges(4,2) * t75 - mrSges(3,1);
t150 = mrSges(3,2) - mrSges(4,3);
t113 = t73 * t133;
t103 = -t134 * t113 - t56 * t75;
t148 = t158 * t103;
t112 = t73 * t132;
t97 = t119 * t132;
t58 = t134 * t78 - t76 * t97;
t32 = -t112 * t133 + t58 * t75;
t147 = t158 * t32;
t129 = t73 * t76;
t53 = t119 * t133 - t129 * t75;
t146 = t158 * t53;
t145 = t152 * pkin(4) + t94;
t143 = t144 * t74 - t77 * t94 - mrSges(4,1);
t142 = mrSges(4,2) - m(7) * (pkin(10) - qJ(6)) + t155;
t140 = -t157 * t75 - t151;
t137 = pkin(10) * t32;
t135 = t103 * pkin(10);
t128 = t73 * t78;
t123 = pkin(2) * t128 + pkin(9) * t129;
t122 = t134 * pkin(1) + pkin(8) * t112;
t118 = t75 * t128;
t111 = t74 * t133;
t110 = t77 * t133;
t109 = t78 * t133;
t108 = -t55 * pkin(2) + pkin(9) * t56;
t57 = t134 * t76 + t78 * t97;
t107 = -t57 * pkin(2) + pkin(9) * t58;
t22 = t103 * pkin(3);
t106 = pkin(10) * t29 + t22;
t24 = t32 * pkin(3);
t33 = t112 * t75 + t133 * t58;
t105 = pkin(10) * t33 - t24;
t48 = t53 * pkin(3);
t54 = t113 * t76 + t119 * t75;
t104 = pkin(10) * t54 + t48;
t101 = t73 * t109;
t102 = pkin(3) * t101 + pkin(10) * t118 + t123;
t100 = -pkin(1) * t132 + pkin(8) * t114;
t93 = t154 * t55 + t108;
t92 = t154 * t57 + t107;
t91 = t58 * pkin(2) + pkin(9) * t57 + t122;
t90 = t33 * pkin(3) + t91;
t88 = -t152 * qJ(5) + t144;
t84 = -t56 * pkin(2) - t55 * pkin(9) + t100;
t81 = -pkin(3) * t29 + t84;
t10 = t33 * t77 + t57 * t74;
t9 = t33 * t74 - t57 * t77;
t80 = t10 * pkin(4) + qJ(5) * t9 + t90;
t79 = -pkin(4) * t153 - qJ(5) * t5 + t81;
t36 = (t109 * t77 + t74 * t76) * t73;
t35 = t101 * t74 - t129 * t77;
t26 = t128 * t77 + t54 * t74;
t16 = -t110 * t57 + t58 * t74;
t15 = -t111 * t57 - t58 * t77;
t14 = -t110 * t55 + t56 * t74;
t13 = -t111 * t55 - t56 * t77;
t1 = [(-t134 * mrSges(2,1) + t132 * mrSges(2,2) - m(3) * t122 - t58 * mrSges(3,1) - mrSges(3,3) * t112 - m(4) * t91 - t33 * mrSges(4,1) - m(5) * (t90 + t137) - m(6) * (t80 + t137) - m(7) * t80 + t150 * t57 + t144 * t9 - t94 * t10 + t142 * t32) * g(2) + (t132 * mrSges(2,1) + t134 * mrSges(2,2) - m(3) * t100 + t56 * mrSges(3,1) - mrSges(3,3) * t114 - m(4) * t84 + t29 * mrSges(4,1) - m(5) * (t81 + t135) - m(6) * (t79 + t135) - m(7) * t79 - t150 * t55 + t94 * t153 - t144 * t5 + t142 * t103) * g(1) (-m(4) * t123 - m(5) * t102 - t152 * (t36 * pkin(4) + t35 * qJ(5) + t102) + (t150 * t76 + t151 * t78) * t73 - t94 * t36 + t144 * t35 + t157 * t118) * g(3) + (-m(4) * t108 - m(5) * t93 - t152 * (t14 * pkin(4) + qJ(5) * t13 + t93) + t150 * t56 - t94 * t14 + t144 * t13 + t140 * t55) * g(2) + (-m(4) * t107 - m(5) * t92 - t152 * (t16 * pkin(4) + qJ(5) * t15 + t92) + t150 * t58 - t94 * t16 + t144 * t15 + t140 * t57) * g(1) (-m(5) * t104 - m(6) * (t104 + t146) - m(7) * (t48 + t146) + t142 * t54 + t143 * t53) * g(3) + (-m(5) * t106 - m(6) * (t106 + t148) - m(7) * (t22 + t148) + t142 * t29 + t143 * t103) * g(2) + (-m(5) * t105 - m(6) * (t105 - t147) - m(7) * (-t24 - t147) + t142 * t33 - t143 * t32) * g(1) (t88 * (-t74 * t128 + t54 * t77) + t145 * t26) * g(3) + (t145 * t5 + t88 * t153) * g(2) + (t10 * t88 + t145 * t9) * g(1), t152 * (-g(1) * t9 - g(2) * t5 - g(3) * t26) (g(1) * t32 - g(2) * t103 - g(3) * t53) * m(7)];
taug  = t1(:);
