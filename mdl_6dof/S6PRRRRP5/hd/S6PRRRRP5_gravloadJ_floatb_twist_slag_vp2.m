% Calculate Gravitation load on the joints for
% S6PRRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-03-09 00:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRRRP5_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP5_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP5_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP5_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP5_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRP5_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:20:14
% EndTime: 2019-03-09 00:20:18
% DurationCPUTime: 1.44s
% Computational Cost: add. (1174->158), mult. (3223->241), div. (0->0), fcn. (4076->14), ass. (0->87)
t151 = mrSges(6,1) + mrSges(7,1);
t146 = -mrSges(6,2) - mrSges(7,2);
t73 = cos(qJ(5));
t68 = pkin(5) * t73 + pkin(4);
t150 = m(6) * pkin(4) + m(7) * t68 + mrSges(5,1);
t149 = m(6) * pkin(11) - m(7) * (-qJ(6) - pkin(11)) - mrSges(5,2) + mrSges(6,3) + mrSges(7,3);
t135 = m(7) * pkin(5);
t147 = mrSges(4,2) - mrSges(5,3);
t123 = cos(pkin(12));
t124 = cos(pkin(7));
t125 = cos(pkin(6));
t102 = t125 * t123;
t120 = sin(pkin(12));
t132 = sin(qJ(2));
t134 = cos(qJ(2));
t80 = -t102 * t134 + t120 * t132;
t121 = sin(pkin(7));
t122 = sin(pkin(6));
t98 = t122 * t121;
t145 = t123 * t98 + t80 * t124;
t100 = t125 * t120;
t81 = t100 * t134 + t123 * t132;
t97 = t122 * t120;
t144 = -t121 * t97 + t81 * t124;
t143 = -t121 * mrSges(4,3) + mrSges(3,2);
t99 = t124 * t122;
t142 = t125 * t121 + t134 * t99;
t141 = -m(5) - m(6) - m(7);
t70 = sin(qJ(5));
t140 = t146 * t70 + t151 * t73 + t150;
t139 = -t135 - t151;
t71 = sin(qJ(4));
t74 = cos(qJ(4));
t138 = t149 * t71 + t150 * t74 + mrSges(4,1);
t137 = -m(7) * (pkin(5) * t70 + pkin(10)) + t147;
t133 = cos(qJ(3));
t58 = t102 * t132 + t120 * t134;
t72 = sin(qJ(3));
t30 = t133 * t58 - t145 * t72;
t131 = t30 * t70;
t59 = -t100 * t132 + t123 * t134;
t32 = t59 * t133 - t144 * t72;
t130 = t32 * t70;
t105 = t122 * t132;
t47 = t133 * t105 + t142 * t72;
t129 = t47 * t70;
t128 = t70 * t74;
t127 = t73 * t74;
t106 = t134 * t122;
t88 = t132 * t98;
t126 = pkin(2) * t106 + pkin(9) * t88;
t90 = t132 * t99;
t54 = t106 * t133 - t72 * t90;
t119 = t54 * pkin(3) + t126;
t112 = pkin(9) * t121;
t111 = t71 * t121;
t110 = t72 * t124;
t109 = t74 * t121;
t107 = t124 * t133;
t53 = t106 * t72 + t133 * t90;
t103 = pkin(10) * t53 + t119;
t96 = -t80 * pkin(2) + t112 * t58;
t95 = -t81 * pkin(2) + t112 * t59;
t38 = -t110 * t58 - t133 * t80;
t94 = t38 * pkin(3) + t96;
t40 = -t110 * t59 - t133 * t81;
t93 = t40 * pkin(3) + t95;
t37 = t107 * t58 - t72 * t80;
t84 = t37 * pkin(10) + t94;
t39 = t107 * t59 - t72 * t81;
t83 = t39 * pkin(10) + t93;
t79 = t124 * t125 - t134 * t98;
t76 = t121 * t81 + t124 * t97;
t75 = t121 * t80 - t123 * t99;
t46 = t105 * t72 - t142 * t133;
t42 = t54 * t74 + t71 * t88;
t34 = t47 * t74 + t71 * t79;
t33 = t47 * t71 - t74 * t79;
t31 = t144 * t133 + t59 * t72;
t29 = t145 * t133 + t58 * t72;
t24 = t111 * t59 + t40 * t74;
t22 = t111 * t58 + t38 * t74;
t18 = t32 * t74 + t71 * t76;
t17 = t32 * t71 - t74 * t76;
t16 = t30 * t74 + t71 * t75;
t15 = t30 * t71 - t74 * t75;
t1 = [(-m(2) - m(3) - m(4) + t141) * g(3) (-mrSges(3,1) * t106 + mrSges(3,2) * t105 - m(4) * t126 - t54 * mrSges(4,1) - mrSges(4,3) * t88 - m(5) * t103 - t42 * mrSges(5,1) - m(6) * (pkin(4) * t42 + t103) - m(7) * (t42 * t68 + t119) + t137 * t53 - t151 * (t42 * t73 + t53 * t70) + t146 * (-t42 * t70 + t53 * t73) - t149 * (t54 * t71 - t74 * t88)) * g(3) + (mrSges(3,1) * t80 - m(4) * t96 - t38 * mrSges(4,1) - m(5) * t84 - t22 * mrSges(5,1) - m(6) * (t22 * pkin(4) + t84) - m(7) * (t22 * t68 + t94) + t146 * (-t22 * t70 + t37 * t73) + t143 * t58 + t137 * t37 - t151 * (t22 * t73 + t37 * t70) - t149 * (-t109 * t58 + t38 * t71)) * g(2) + (mrSges(3,1) * t81 - m(4) * t95 - t40 * mrSges(4,1) - m(5) * t83 - t24 * mrSges(5,1) - m(6) * (t24 * pkin(4) + t83) - m(7) * (t24 * t68 + t93) + t143 * t59 + t137 * t39 - t151 * (t24 * t73 + t39 * t70) + t146 * (-t24 * t70 + t39 * t73) - t149 * (-t109 * t59 + t40 * t71)) * g(1) (-t129 * t135 + t147 * t47 - t151 * (-t127 * t46 + t129) + t146 * (t128 * t46 + t47 * t73) + t141 * (-t46 * pkin(3) + t47 * pkin(10)) + t138 * t46) * g(3) + (-t131 * t135 - t151 * (-t127 * t29 + t131) + t146 * (t128 * t29 + t30 * t73) + t147 * t30 + t141 * (-t29 * pkin(3) + t30 * pkin(10)) + t138 * t29) * g(2) + (-t130 * t135 - t151 * (-t127 * t31 + t130) + t146 * (t128 * t31 + t32 * t73) + t147 * t32 + t141 * (-t31 * pkin(3) + t32 * pkin(10)) + t138 * t31) * g(1) (t140 * t33 - t149 * t34) * g(3) + (t140 * t15 - t149 * t16) * g(2) + (t140 * t17 - t149 * t18) * g(1) (t146 * (-t34 * t73 - t46 * t70) + t139 * (-t34 * t70 + t46 * t73)) * g(3) + (t146 * (-t16 * t73 - t29 * t70) + t139 * (-t16 * t70 + t29 * t73)) * g(2) + (t146 * (-t18 * t73 - t31 * t70) + t139 * (-t18 * t70 + t31 * t73)) * g(1) (-g(1) * t17 - g(2) * t15 - g(3) * t33) * m(7)];
taug  = t1(:);
