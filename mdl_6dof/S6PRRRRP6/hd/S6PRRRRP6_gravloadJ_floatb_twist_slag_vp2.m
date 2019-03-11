% Calculate Gravitation load on the joints for
% S6PRRRRP6
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
% Datum: 2019-03-09 00:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRRRP6_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP6_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP6_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP6_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP6_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRP6_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:28:39
% EndTime: 2019-03-09 00:28:41
% DurationCPUTime: 1.17s
% Computational Cost: add. (1324->159), mult. (3680->248), div. (0->0), fcn. (4689->14), ass. (0->90)
t122 = m(7) * pkin(5) + mrSges(6,1) + mrSges(7,1);
t121 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t156 = -m(6) - m(7);
t164 = mrSges(4,2) - mrSges(5,3);
t163 = mrSges(6,3) + mrSges(7,2);
t91 = sin(qJ(4));
t94 = cos(qJ(4));
t162 = t94 * mrSges(5,1) - mrSges(5,2) * t91 + mrSges(4,1);
t138 = sin(pkin(6));
t139 = cos(pkin(12));
t109 = t139 * t138;
t140 = cos(pkin(7));
t89 = sin(pkin(7));
t141 = cos(pkin(6));
t112 = t141 * t139;
t137 = sin(pkin(12));
t152 = sin(qJ(2));
t154 = cos(qJ(2));
t97 = -t112 * t154 + t137 * t152;
t161 = t89 * t109 + t97 * t140;
t108 = t138 * t137;
t111 = t141 * t137;
t98 = t111 * t154 + t139 * t152;
t160 = t89 * t108 - t98 * t140;
t110 = t140 * t138;
t159 = t154 * t110 + t141 * t89;
t158 = mrSges(5,2) - t163;
t90 = sin(qJ(5));
t93 = cos(qJ(5));
t157 = t121 * t90 - t122 * t93 - mrSges(5,1);
t155 = pkin(4) * t94;
t153 = cos(qJ(3));
t79 = t112 * t152 + t137 * t154;
t92 = sin(qJ(3));
t42 = t161 * t153 + t79 * t92;
t151 = t42 * t91;
t80 = -t111 * t152 + t139 * t154;
t44 = -t160 * t153 + t80 * t92;
t150 = t44 * t91;
t118 = t138 * t152;
t64 = t118 * t92 - t159 * t153;
t149 = t64 * t91;
t148 = t89 * t91;
t147 = t89 * t94;
t146 = t90 * t94;
t145 = t93 * t94;
t125 = t92 * t140;
t53 = -t125 * t79 - t153 * t97;
t76 = t97 * pkin(2);
t144 = t53 * pkin(3) - t76;
t55 = -t125 * t80 - t153 * t98;
t77 = t98 * pkin(2);
t143 = t55 * pkin(3) - t77;
t107 = t89 * t118;
t119 = t154 * t138;
t142 = pkin(2) * t119 + pkin(9) * t107;
t136 = -m(5) + t156;
t102 = t152 * t110;
t75 = -t102 * t92 + t119 * t153;
t133 = t75 * pkin(3) + t142;
t132 = -m(4) + t136;
t43 = t153 * t79 - t161 * t92;
t131 = -t42 * pkin(3) + pkin(10) * t43;
t45 = t80 * t153 + t160 * t92;
t130 = -t44 * pkin(3) + pkin(10) * t45;
t65 = t153 * t118 + t159 * t92;
t129 = -t64 * pkin(3) + pkin(10) * t65;
t120 = t140 * t153;
t104 = pkin(10) * t136 + t164;
t103 = pkin(11) * t156 + t158;
t100 = mrSges(3,2) + (pkin(9) * t132 - mrSges(4,3)) * t89;
t78 = -t119 * t89 + t140 * t141;
t74 = t102 * t153 + t119 * t92;
t67 = t108 * t140 + t89 * t98;
t66 = -t109 * t140 + t89 * t97;
t58 = t107 * t91 + t75 * t94;
t54 = t120 * t80 - t92 * t98;
t52 = t120 * t79 - t92 * t97;
t47 = t65 * t94 + t78 * t91;
t46 = -t65 * t91 + t78 * t94;
t28 = t148 * t80 + t55 * t94;
t26 = t148 * t79 + t53 * t94;
t20 = t45 * t94 + t67 * t91;
t19 = -t45 * t91 + t67 * t94;
t18 = t43 * t94 + t66 * t91;
t17 = -t43 * t91 + t66 * t94;
t15 = t47 * t90 - t64 * t93;
t3 = t20 * t90 - t44 * t93;
t1 = t18 * t90 - t42 * t93;
t2 = [(-m(2) - m(3) + t132) * g(3) (-mrSges(3,1) * t119 + mrSges(3,2) * t118 - m(4) * t142 - t75 * mrSges(4,1) - mrSges(4,3) * t107 - m(5) * t133 - t58 * mrSges(5,1) + t104 * t74 - t122 * (t58 * t93 + t74 * t90) + t121 * (t58 * t90 - t74 * t93) + t103 * (-t107 * t94 + t75 * t91) + t156 * (t58 * pkin(4) + t133)) * g(3) + (t97 * mrSges(3,1) + m(4) * t76 - t53 * mrSges(4,1) - m(5) * t144 - t26 * mrSges(5,1) + t121 * (t26 * t90 - t52 * t93) + t100 * t79 + t104 * t52 - t122 * (t26 * t93 + t52 * t90) + t103 * (-t147 * t79 + t53 * t91) + t156 * (t26 * pkin(4) + t144)) * g(2) + (t98 * mrSges(3,1) + m(4) * t77 - t55 * mrSges(4,1) - m(5) * t143 - t28 * mrSges(5,1) + t100 * t80 + t104 * t54 - t122 * (t28 * t93 + t54 * t90) + t121 * (t28 * t90 - t54 * t93) + t103 * (-t147 * t80 + t55 * t91) + t156 * (t28 * pkin(4) + t143)) * g(1) (-m(5) * t129 + t164 * t65 + t162 * t64 + t156 * (-pkin(11) * t149 - t64 * t155 + t129) - t122 * (-t145 * t64 + t65 * t90) + t121 * (-t146 * t64 - t65 * t93) + t163 * t149) * g(3) + (-m(5) * t131 + t156 * (-pkin(11) * t151 - t42 * t155 + t131) - t122 * (-t145 * t42 + t43 * t90) + t121 * (-t146 * t42 - t43 * t93) + t164 * t43 + t162 * t42 + t163 * t151) * g(2) + (-m(5) * t130 + t156 * (-pkin(11) * t150 - t44 * t155 + t130) - t122 * (-t145 * t44 + t45 * t90) + t121 * (-t146 * t44 - t45 * t93) + t164 * t45 + t162 * t44 + t163 * t150) * g(1) (t156 * (t46 * pkin(4) + pkin(11) * t47) + t158 * t47 + t157 * t46) * g(3) + (t156 * (t17 * pkin(4) + pkin(11) * t18) + t158 * t18 + t157 * t17) * g(2) + (t156 * (t19 * pkin(4) + pkin(11) * t20) + t158 * t20 + t157 * t19) * g(1) (t121 * (t47 * t93 + t64 * t90) + t122 * t15) * g(3) + (t121 * (t18 * t93 + t42 * t90) + t122 * t1) * g(2) + (t121 * (t20 * t93 + t44 * t90) + t122 * t3) * g(1) (-g(1) * t3 - g(2) * t1 - g(3) * t15) * m(7)];
taug  = t2(:);
