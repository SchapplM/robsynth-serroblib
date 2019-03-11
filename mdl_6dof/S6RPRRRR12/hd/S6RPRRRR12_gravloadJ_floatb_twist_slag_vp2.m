% Calculate Gravitation load on the joints for
% S6RPRRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 08:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRRR12_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR12_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR12_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RPRRRR12_gravloadJ_floatb_twist_slag_vp2: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR12_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR12_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:50:37
% EndTime: 2019-03-09 07:50:42
% DurationCPUTime: 2.16s
% Computational Cost: add. (2462->166), mult. (6956->251), div. (0->0), fcn. (9078->18), ass. (0->90)
t157 = cos(qJ(4));
t144 = cos(pkin(8));
t142 = sin(pkin(6));
t145 = cos(pkin(7));
t123 = t145 * t142;
t141 = sin(pkin(7));
t159 = cos(qJ(1));
t143 = cos(pkin(14));
t146 = cos(pkin(6));
t126 = t146 * t143;
t140 = sin(pkin(14));
t156 = sin(qJ(1));
t70 = -t126 * t159 + t140 * t156;
t108 = -t159 * t123 + t70 * t141;
t79 = sin(pkin(8));
t153 = t108 * t79;
t122 = t142 * t141;
t158 = cos(qJ(3));
t116 = t158 * t122;
t137 = t70 * t145;
t155 = sin(qJ(3));
t124 = t146 * t140;
t71 = t124 * t159 + t143 * t156;
t54 = t116 * t159 + t137 * t158 + t155 * t71;
t177 = t144 * t54 - t153;
t115 = t155 * t122;
t53 = -t159 * t115 - t137 * t155 + t158 * t71;
t82 = sin(qJ(4));
t22 = -t53 * t157 + t177 * t82;
t38 = t108 * t144 + t54 * t79;
t81 = sin(qJ(5));
t84 = cos(qJ(5));
t183 = t22 * t81 + t38 * t84;
t182 = t22 * t84 - t38 * t81;
t171 = m(6) + m(7);
t120 = -pkin(12) * t171 + mrSges(5,2) - mrSges(6,3);
t80 = sin(qJ(6));
t83 = cos(qJ(6));
t179 = -t80 * mrSges(7,1) - t83 * mrSges(7,2) + t120;
t178 = -m(5) - t171;
t168 = m(7) * pkin(5) + t83 * mrSges(7,1) - t80 * mrSges(7,2) + mrSges(6,1);
t133 = -m(7) * pkin(13) + mrSges(6,2) - mrSges(7,3);
t176 = t53 * t82;
t106 = t126 * t156 + t140 * t159;
t173 = t106 * t141 + t156 * t123;
t100 = t106 * t145;
t107 = -t124 * t156 + t143 * t159;
t90 = t100 * t158 + t107 * t155 - t116 * t156;
t175 = t144 * t173 + t90 * t79;
t172 = pkin(4) * t171 - t133 * t81 + t168 * t84 + mrSges(5,1);
t169 = t90 * t144 - t173 * t79;
t104 = -t122 * t143 + t145 * t146;
t121 = t142 * t140;
t165 = t143 * t123 + t146 * t141;
t94 = t155 * t121 - t165 * t158;
t167 = t104 * t79 - t94 * t144;
t152 = t79 * t81;
t151 = t79 * t84;
t129 = t142 * t156;
t147 = t159 * pkin(1) + qJ(2) * t129;
t135 = t82 * t144;
t130 = t159 * t142;
t132 = -pkin(1) * t156 + qJ(2) * t130;
t131 = t144 * t157;
t111 = mrSges(4,2) + (pkin(11) * t178 - mrSges(5,3)) * t79;
t109 = -t71 * pkin(2) - t108 * pkin(10) + t132;
t101 = -t53 * pkin(3) - pkin(11) * t38 + t109;
t97 = t107 * pkin(2) + t173 * pkin(10) + t147;
t56 = -t100 * t155 + t107 * t158 + t115 * t156;
t87 = t56 * pkin(3) + t175 * pkin(11) + t97;
t24 = t56 * t157 - t169 * t82;
t86 = t24 * pkin(4) + t87;
t65 = t158 * t121 + t165 * t155;
t64 = t94 * pkin(3);
t51 = t90 * pkin(3);
t49 = t54 * pkin(3);
t48 = t104 * t144 + t79 * t94;
t42 = -t135 * t65 - t157 * t94;
t35 = t65 * t157 + t167 * t82;
t34 = -t167 * t157 + t65 * t82;
t30 = -t135 * t56 - t157 * t90;
t28 = -t135 * t53 - t157 * t54;
t23 = t169 * t157 + t56 * t82;
t19 = t157 * t177 + t176;
t14 = t35 * t84 + t48 * t81;
t8 = t175 * t81 + t24 * t84;
t7 = -t175 * t84 + t24 * t81;
t2 = t23 * t80 + t8 * t83;
t1 = t23 * t83 - t8 * t80;
t3 = [(-t159 * mrSges(2,1) + t156 * mrSges(2,2) - m(3) * t147 - t107 * mrSges(3,1) + t106 * mrSges(3,2) - mrSges(3,3) * t129 - m(4) * t97 - t56 * mrSges(4,1) + t90 * mrSges(4,2) - t173 * mrSges(4,3) - m(5) * t87 - t24 * mrSges(5,1) - t175 * mrSges(5,3) - m(6) * t86 - t8 * mrSges(6,1) - m(7) * (t8 * pkin(5) + t86) - t2 * mrSges(7,1) - t1 * mrSges(7,2) + t133 * t7 + t120 * t23) * g(2) + (t156 * mrSges(2,1) + t159 * mrSges(2,2) - m(3) * t132 + t71 * mrSges(3,1) - t70 * mrSges(3,2) - mrSges(3,3) * t130 - m(4) * t109 + t53 * mrSges(4,1) - t54 * mrSges(4,2) + t108 * mrSges(4,3) - m(5) * t101 - t22 * mrSges(5,1) + t38 * mrSges(5,3) + t133 * t183 - t168 * t182 + t179 * (-t131 * t54 + t153 * t157 - t176) + t171 * (-t22 * pkin(4) - t101)) * g(1) (-g(1) * t129 + g(2) * t130 - g(3) * t146) * (m(3) + m(4) - t178) (t94 * mrSges(4,1) + m(5) * t64 - t42 * mrSges(5,1) + t111 * t65 - t168 * (t152 * t65 + t42 * t84) + t179 * (t131 * t65 - t82 * t94) + t133 * (-t151 * t65 + t42 * t81) - t171 * (t42 * pkin(4) - t64)) * g(3) + (t54 * mrSges(4,1) + m(5) * t49 - t28 * mrSges(5,1) + t133 * (-t151 * t53 + t28 * t81) + t111 * t53 - t168 * (t152 * t53 + t28 * t84) + t179 * (t131 * t53 - t54 * t82) - t171 * (t28 * pkin(4) - t49)) * g(2) + (t90 * mrSges(4,1) + m(5) * t51 - t30 * mrSges(5,1) + t111 * t56 - t168 * (t152 * t56 + t30 * t84) + t179 * (t131 * t56 - t82 * t90) + t133 * (-t151 * t56 + t30 * t81) - t171 * (t30 * pkin(4) - t51)) * g(1) (t172 * t34 + t179 * t35) * g(3) + (t172 * t19 - t179 * t22) * g(2) + (t172 * t23 + t179 * t24) * g(1) (t133 * t14 - t168 * (-t35 * t81 + t48 * t84)) * g(3) + (-t133 * t182 - t168 * t183) * g(2) + (t133 * t8 + t168 * t7) * g(1), -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * ((t182 * t80 + t19 * t83) * mrSges(7,1) + (t182 * t83 - t19 * t80) * mrSges(7,2)) - g(3) * ((-t14 * t80 + t34 * t83) * mrSges(7,1) + (-t14 * t83 - t34 * t80) * mrSges(7,2))];
taug  = t3(:);
