% Calculate Gravitation load on the joints for
% S6RRRPRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-03-09 17:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPRP7_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP7_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP7_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP7_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP7_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP7_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:03:56
% EndTime: 2019-03-09 17:04:00
% DurationCPUTime: 1.56s
% Computational Cost: add. (876->143), mult. (1564->195), div. (0->0), fcn. (1833->12), ass. (0->72)
t161 = m(7) * qJ(6) - mrSges(6,2) + mrSges(7,3);
t75 = sin(qJ(5));
t162 = t161 * t75;
t105 = m(7) * pkin(5) + mrSges(6,1) + mrSges(7,1);
t148 = m(6) + m(7);
t146 = -mrSges(6,3) - mrSges(7,2);
t72 = qJ(3) + pkin(11);
t70 = cos(t72);
t160 = t105 * t70;
t152 = mrSges(5,2) + t146;
t144 = -m(4) * pkin(9) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t159 = -t105 * t75 + t144;
t158 = m(5) + t148;
t69 = sin(t72);
t76 = sin(qJ(3));
t80 = cos(qJ(3));
t155 = -m(4) * pkin(2) - t80 * mrSges(4,1) - t70 * mrSges(5,1) + t76 * mrSges(4,2) + t69 * mrSges(5,2) - mrSges(3,1);
t120 = cos(pkin(6));
t73 = sin(pkin(6));
t77 = sin(qJ(2));
t130 = t73 * t77;
t154 = t120 * t80 - t76 * t130;
t128 = t73 * t80;
t78 = sin(qJ(1));
t109 = t78 * t120;
t139 = cos(qJ(1));
t81 = cos(qJ(2));
t50 = -t77 * t109 + t139 * t81;
t23 = t78 * t128 - t50 * t76;
t79 = cos(qJ(5));
t153 = t105 * t79 + mrSges(5,1) + t162;
t113 = t73 * t139;
t100 = t120 * t139;
t48 = t77 * t100 + t78 * t81;
t18 = -t69 * t113 + t48 * t70;
t47 = -t81 * t100 + t77 * t78;
t1 = t18 * t75 - t47 * t79;
t151 = t18 * t79 + t47 * t75;
t150 = t161 * t79 + t159;
t142 = pkin(4) * t70;
t149 = -t148 * (-pkin(10) * t69 - t142) + t79 * t160 + t70 * t162 - t155 - t146 * t69;
t87 = t80 * t113 + t48 * t76;
t83 = t87 * pkin(3);
t129 = t73 * t78;
t127 = t73 * t81;
t124 = t79 * t81;
t121 = t139 * pkin(1) + pkin(8) * t129;
t119 = t69 * t127;
t118 = t76 * t129;
t115 = t75 * t127;
t112 = -pkin(1) * t78 + pkin(8) * t113;
t17 = -t70 * t113 - t48 * t69;
t63 = t76 * t113;
t110 = -t48 * t80 + t63;
t103 = t23 * pkin(3);
t49 = t81 * t109 + t139 * t77;
t68 = pkin(3) * t80 + pkin(2);
t74 = -qJ(4) - pkin(9);
t101 = pkin(3) * t118 - t49 * t74 + t50 * t68 + t121;
t94 = t154 * pkin(3);
t92 = pkin(3) * t63 + t47 * t74 - t48 * t68 + t112;
t90 = -t148 * pkin(10) + t152;
t53 = t68 * t127;
t37 = t120 * t69 + t70 * t130;
t36 = t120 * t70 - t69 * t130;
t24 = t50 * t80 + t118;
t22 = t69 * t129 + t50 * t70;
t21 = -t70 * t129 + t50 * t69;
t15 = t73 * t124 + t37 * t75;
t6 = t22 * t79 + t49 * t75;
t5 = t22 * t75 - t49 * t79;
t2 = [(-t139 * mrSges(2,1) - m(3) * t121 - t50 * mrSges(3,1) - m(4) * (pkin(2) * t50 + t121) - t24 * mrSges(4,1) - t23 * mrSges(4,2) - m(5) * t101 - t22 * mrSges(5,1) + (-mrSges(3,3) * t73 + mrSges(2,2)) * t78 - t105 * t6 - t161 * t5 + t144 * t49 + t90 * t21 - t148 * (t22 * pkin(4) + t101)) * g(2) + (t78 * mrSges(2,1) + t139 * mrSges(2,2) - m(3) * t112 + t48 * mrSges(3,1) - mrSges(3,3) * t113 - m(4) * (-pkin(2) * t48 + t112) - t110 * mrSges(4,1) - t87 * mrSges(4,2) - m(5) * t92 + t18 * mrSges(5,1) + t105 * t151 + t161 * t1 - t144 * t47 + t90 * t17 + t148 * (pkin(4) * t18 - t92)) * g(1) (-t158 * (-t47 * t68 - t48 * t74) + t150 * t48 + t149 * t47) * g(2) + (-t158 * (-t49 * t68 - t50 * t74) + t150 * t50 + t149 * t49) * g(1) + (-m(5) * t53 - t148 * (pkin(10) * t119 + t127 * t142 - t74 * t130 + t53) - t161 * (t70 * t115 - t79 * t130) + t146 * t119 + (-t124 * t160 + t155 * t81 + (m(5) * t74 + t159) * t77) * t73) * g(3) (-t154 * mrSges(4,1) - (-t120 * t76 - t77 * t128) * mrSges(4,2) - m(5) * t94 - t148 * (t36 * pkin(4) + pkin(10) * t37 + t94) + t152 * t37 - t153 * t36) * g(3) + (m(5) * t83 + t87 * mrSges(4,1) - t110 * mrSges(4,2) + t152 * t18 - t153 * t17 - t148 * (t17 * pkin(4) + t18 * pkin(10) - t83)) * g(2) + (-m(5) * t103 - mrSges(4,1) * t23 + mrSges(4,2) * t24 - t148 * (-t21 * pkin(4) + pkin(10) * t22 + t103) + t152 * t22 + t153 * t21) * g(1), t158 * (-g(1) * t49 - g(2) * t47 + g(3) * t127) (-t161 * (t37 * t79 - t115) + t105 * t15) * g(3) + (t105 * t1 - t151 * t161) * g(2) + (t105 * t5 - t161 * t6) * g(1) (-g(1) * t5 - g(2) * t1 - g(3) * t15) * m(7)];
taug  = t2(:);
