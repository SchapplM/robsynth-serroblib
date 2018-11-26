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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 18:10
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

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
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP8_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP8_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:09:50
% EndTime: 2018-11-23 18:09:52
% DurationCPUTime: 1.47s
% Computational Cost: add. (1957->156), mult. (2447->199), div. (0->0), fcn. (2509->14), ass. (0->90)
t168 = -mrSges(5,2) + mrSges(7,2) + mrSges(6,3);
t103 = m(7) * pkin(5) + mrSges(5,1) + mrSges(6,1) + mrSges(7,1);
t79 = sin(qJ(4));
t82 = cos(qJ(4));
t165 = pkin(4) * t82 + qJ(5) * t79;
t164 = mrSges(4,2) + mrSges(7,3) - mrSges(6,2) - mrSges(5,3);
t163 = t103 * t82 + t168 * t79 + mrSges(4,1);
t142 = cos(qJ(3));
t80 = sin(qJ(3));
t162 = t142 * pkin(3) + pkin(10) * t80;
t131 = sin(pkin(6));
t143 = cos(qJ(1));
t109 = t143 * t131;
t141 = sin(qJ(1));
t83 = cos(qJ(2));
t122 = t141 * t83;
t130 = pkin(6) - qJ(2);
t113 = sin(t130);
t129 = pkin(6) + qJ(2);
t112 = sin(t129);
t75 = t112 / 0.2e1;
t98 = t75 - t113 / 0.2e1;
t56 = t143 * t98 + t122;
t30 = -t80 * t109 + t142 * t56;
t81 = sin(qJ(2));
t107 = cos(t129) / 0.2e1;
t114 = cos(t130);
t86 = t114 / 0.2e1 + t107;
t55 = t141 * t81 - t143 * t86;
t5 = t30 * t79 - t55 * t82;
t161 = t30 * t82 + t55 * t79;
t160 = m(6) + m(7);
t159 = mrSges(3,2) - mrSges(4,3);
t115 = -t142 * t109 - t56 * t80;
t158 = t165 * t115;
t108 = t131 * t141;
t123 = t143 * t83;
t59 = -t141 * t98 + t123;
t33 = -t108 * t142 + t59 * t80;
t157 = t165 * t33;
t69 = t107 - t114 / 0.2e1;
t78 = cos(pkin(6));
t53 = t142 * t78 + t69 * t80;
t156 = t165 * t53;
t155 = t160 * pkin(4) + t103;
t152 = -m(7) * (pkin(10) - qJ(6)) + t164;
t149 = -mrSges(4,1) * t142 - mrSges(3,1) + (m(7) * qJ(6) + t164) * t80;
t146 = pkin(10) * t33;
t144 = t115 * pkin(10);
t134 = t143 * pkin(1) + pkin(8) * t108;
t125 = t79 * t142;
t124 = t82 * t142;
t106 = t113 / 0.2e1;
t85 = t106 - t112 / 0.2e1;
t57 = -t143 * t85 + t122;
t121 = -t55 * pkin(2) + pkin(9) * t57;
t58 = t141 * t86 + t143 * t81;
t60 = t141 * t85 + t123;
t120 = -t58 * pkin(2) + pkin(9) * t60;
t68 = t75 + t106;
t119 = t68 * pkin(2) - pkin(9) * t69;
t24 = t115 * pkin(3);
t118 = pkin(10) * t30 + t24;
t26 = t33 * pkin(3);
t34 = t108 * t80 + t142 * t59;
t117 = pkin(10) * t34 - t26;
t48 = t53 * pkin(3);
t54 = -t142 * t69 + t78 * t80;
t116 = pkin(10) * t54 + t48;
t111 = -pkin(1) * t141 + pkin(8) * t109;
t102 = -t162 * t55 + t121;
t101 = -t162 * t58 + t120;
t100 = t59 * pkin(2) + pkin(9) * t58 + t134;
t99 = t162 * t68 + t119;
t97 = t34 * pkin(3) + t100;
t95 = -t160 * qJ(5) - t168;
t92 = -t56 * pkin(2) - t55 * pkin(9) + t111;
t88 = -pkin(3) * t30 + t92;
t10 = t34 * t82 + t58 * t79;
t9 = t34 * t79 - t58 * t82;
t87 = t10 * pkin(4) + qJ(5) * t9 + t97;
t84 = -pkin(4) * t161 - qJ(5) * t5 + t88;
t36 = t124 * t68 - t69 * t79;
t35 = t125 * t68 + t69 * t82;
t18 = t54 * t79 + t68 * t82;
t16 = -t124 * t58 + t60 * t79;
t15 = -t125 * t58 - t60 * t82;
t14 = -t124 * t55 + t57 * t79;
t13 = -t125 * t55 - t57 * t82;
t1 = [(-mrSges(2,1) * t143 + mrSges(2,2) * t141 - m(3) * t134 - t59 * mrSges(3,1) - mrSges(3,3) * t108 - m(4) * t100 - t34 * mrSges(4,1) - m(5) * (t97 + t146) - m(6) * (t87 + t146) - m(7) * t87 + t159 * t58 - t168 * t9 - t103 * t10 + t152 * t33) * g(2) + (mrSges(2,1) * t141 + mrSges(2,2) * t143 - m(3) * t111 + t56 * mrSges(3,1) - mrSges(3,3) * t109 - m(4) * t92 + t30 * mrSges(4,1) - m(5) * (t88 + t144) - m(6) * (t84 + t144) - m(7) * t84 - t159 * t55 + t103 * t161 + t168 * t5 + t152 * t115) * g(1) (-m(4) * t119 - m(5) * t99 - t160 * (t36 * pkin(4) + qJ(5) * t35 + t99) - t159 * t69 - t103 * t36 - t168 * t35 + t149 * t68) * g(3) + (-m(4) * t121 - m(5) * t102 - t160 * (t14 * pkin(4) + qJ(5) * t13 + t102) + t159 * t57 - t103 * t14 - t168 * t13 - t149 * t55) * g(2) + (-m(4) * t120 - m(5) * t101 - t160 * (t16 * pkin(4) + qJ(5) * t15 + t101) + t159 * t60 - t103 * t16 - t168 * t15 - t149 * t58) * g(1) (-m(5) * t116 - m(6) * (t116 + t156) - m(7) * (t48 + t156) + t152 * t54 - t163 * t53) * g(3) + (-m(5) * t118 - m(6) * (t118 + t158) - m(7) * (t24 + t158) + t152 * t30 - t163 * t115) * g(2) + (-m(5) * t117 - m(6) * (t117 - t157) - m(7) * (-t26 - t157) + t152 * t34 + t163 * t33) * g(1) (t95 * (t54 * t82 - t68 * t79) + t155 * t18) * g(3) + (t155 * t5 + t95 * t161) * g(2) + (t10 * t95 + t155 * t9) * g(1), t160 * (-g(1) * t9 - g(2) * t5 - g(3) * t18) (g(1) * t33 - g(2) * t115 - g(3) * t53) * m(7)];
taug  = t1(:);
