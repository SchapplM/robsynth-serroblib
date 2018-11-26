% Calculate Gravitation load on the joints for
% S6RRRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2018-11-23 17:55
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRRPRR7_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR7_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR7_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR7_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR7_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR7_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:54:54
% EndTime: 2018-11-23 17:54:55
% DurationCPUTime: 1.21s
% Computational Cost: add. (1541->170), mult. (1470->216), div. (0->0), fcn. (1423->18), ass. (0->86)
t157 = mrSges(6,2) - mrSges(7,3);
t84 = sin(qJ(6));
t88 = cos(qJ(6));
t156 = t88 * mrSges(7,1) - t84 * mrSges(7,2) + mrSges(6,1);
t158 = m(7) * pkin(5) + t156;
t114 = -m(7) * pkin(11) + t157;
t154 = m(6) + m(7);
t153 = m(5) * pkin(3) + mrSges(4,1);
t83 = -qJ(4) - pkin(9);
t93 = -m(4) * pkin(9) + m(5) * t83 + mrSges(3,2) - mrSges(4,3) - mrSges(5,3) - mrSges(6,3);
t152 = t84 * mrSges(7,1) + t88 * mrSges(7,2) - t93;
t129 = pkin(6) + qJ(2);
t112 = cos(t129) / 0.2e1;
t130 = pkin(6) - qJ(2);
t117 = cos(t130);
t60 = t112 - t117 / 0.2e1;
t80 = qJ(3) + pkin(12);
t76 = qJ(5) + t80;
t71 = sin(t76);
t72 = cos(t76);
t82 = cos(pkin(6));
t34 = t60 * t71 + t72 * t82;
t35 = -t60 * t72 + t71 * t82;
t151 = -t156 * t34 + t157 * t35;
t81 = sin(pkin(6));
t87 = sin(qJ(1));
t141 = t81 * t87;
t115 = sin(t129);
t110 = t115 / 0.2e1;
t116 = sin(t130);
t101 = t110 - t116 / 0.2e1;
t142 = cos(qJ(1));
t90 = cos(qJ(2));
t125 = t142 * t90;
t46 = -t101 * t87 + t125;
t17 = -t141 * t72 + t46 * t71;
t18 = t141 * t71 + t46 * t72;
t150 = t156 * t17 + t157 * t18;
t126 = t81 * t142;
t139 = t87 * t90;
t43 = t101 * t142 + t139;
t13 = -t72 * t126 - t43 * t71;
t14 = -t71 * t126 + t43 * t72;
t149 = -t156 * t13 + t157 * t14;
t89 = cos(qJ(3));
t77 = t89 * pkin(3);
t73 = t77 + pkin(2);
t74 = sin(t80);
t75 = cos(t80);
t85 = sin(qJ(3));
t148 = m(4) * pkin(2) + m(5) * t73 + t89 * mrSges(4,1) + t75 * mrSges(5,1) - t85 * mrSges(4,2) - t74 * mrSges(5,2) - t114 * t71 + t158 * t72 + mrSges(3,1);
t144 = pkin(3) * t85;
t61 = pkin(4) * t74 + t144;
t62 = pkin(4) * t75 + t77;
t135 = t62 * t141 - t46 * t61;
t132 = t60 * t61 + t82 * t62;
t131 = t142 * pkin(1) + pkin(8) * t141;
t123 = -t87 * pkin(1) + pkin(8) * t126;
t122 = t13 * pkin(5) + t14 * pkin(11);
t121 = -t17 * pkin(5) + pkin(11) * t18;
t120 = t34 * pkin(5) + pkin(11) * t35;
t119 = t74 * t126 - t43 * t75;
t64 = t85 * t126;
t118 = -t43 * t89 + t64;
t86 = sin(qJ(2));
t92 = t117 / 0.2e1 + t112;
t45 = t142 * t86 + t87 * t92;
t58 = pkin(2) + t62;
t79 = -pkin(10) + t83;
t113 = t61 * t141 - t45 * t79 + t46 * t58 + t131;
t111 = t116 / 0.2e1;
t103 = -t126 * t62 - t43 * t61;
t102 = t111 - t115 / 0.2e1;
t21 = t141 * t89 - t46 * t85;
t96 = t126 * t75 + t43 * t74;
t95 = t126 * t89 + t43 * t85;
t59 = t110 + t111;
t47 = t102 * t87 + t125;
t44 = -t102 * t142 + t139;
t42 = -t142 * t92 + t86 * t87;
t22 = t141 * t85 + t46 * t89;
t20 = t141 * t74 + t46 * t75;
t19 = t141 * t75 - t46 * t74;
t2 = t18 * t88 + t45 * t84;
t1 = -t18 * t84 + t45 * t88;
t3 = [(-t142 * mrSges(2,1) - m(3) * t131 - t46 * mrSges(3,1) - m(4) * (pkin(2) * t46 + t131) - t22 * mrSges(4,1) - t21 * mrSges(4,2) - m(5) * (t46 * t73 + t131) - t20 * mrSges(5,1) - t19 * mrSges(5,2) - m(6) * t113 - t18 * mrSges(6,1) - m(7) * (pkin(5) * t18 + t113) - t2 * mrSges(7,1) - t1 * mrSges(7,2) + (mrSges(2,2) + (-m(5) * t144 - mrSges(3,3)) * t81) * t87 + t114 * t17 + t93 * t45) * g(2) + (t87 * mrSges(2,1) + t142 * mrSges(2,2) - m(3) * t123 + t43 * mrSges(3,1) - mrSges(3,3) * t126 - m(4) * (-pkin(2) * t43 + t123) - t118 * mrSges(4,1) - t95 * mrSges(4,2) - m(5) * (pkin(3) * t64 - t43 * t73 + t123) - t119 * mrSges(5,1) - t96 * mrSges(5,2) + t114 * t13 + t158 * t14 + t152 * t42 + t154 * (-t61 * t126 - t42 * t79 + t43 * t58 - t123)) * g(1) (-t154 * (t59 * t58 + t60 * t79) + t152 * t60 - t148 * t59) * g(3) + (-t154 * (-t42 * t58 - t44 * t79) - t152 * t44 + t148 * t42) * g(2) + (-t154 * (-t45 * t58 - t47 * t79) - t152 * t47 + t148 * t45) * g(1) (-(t60 * t89 - t82 * t85) * mrSges(4,2) - (t60 * t74 + t75 * t82) * mrSges(5,1) - (t60 * t75 - t74 * t82) * mrSges(5,2) - m(6) * t132 - m(7) * (t120 + t132) - t153 * (t60 * t85 + t82 * t89) + t151) * g(3) + (-t118 * mrSges(4,2) + t96 * mrSges(5,1) - t119 * mrSges(5,2) - m(6) * t103 - m(7) * (t103 + t122) + t153 * t95 + t149) * g(2) + (mrSges(4,2) * t22 - t19 * mrSges(5,1) + t20 * mrSges(5,2) - m(6) * t135 - m(7) * (t121 + t135) - t153 * t21 + t150) * g(1) (m(5) + t154) * (-g(1) * t45 - g(2) * t42 + g(3) * t59) (-m(7) * t120 + t151) * g(3) + (-m(7) * t122 + t149) * g(2) + (-m(7) * t121 + t150) * g(1), -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * ((-t14 * t84 + t42 * t88) * mrSges(7,1) + (-t14 * t88 - t42 * t84) * mrSges(7,2)) - g(3) * ((-t35 * t84 - t59 * t88) * mrSges(7,1) + (-t35 * t88 + t59 * t84) * mrSges(7,2))];
taug  = t3(:);
