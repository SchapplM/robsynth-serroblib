% Calculate Gravitation load on the joints for
% S6PRRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1,theta5]';
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
% Datum: 2018-11-23 15:25
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6PRRRPR5_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR5_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR5_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR5_gravloadJ_floatb_twist_slag_vp2: pkin has to be [13x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR5_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR5_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:24:44
% EndTime: 2018-11-23 15:24:46
% DurationCPUTime: 1.35s
% Computational Cost: add. (3309->144), mult. (3309->209), div. (0->0), fcn. (3195->24), ass. (0->80)
t81 = sin(qJ(6));
t84 = cos(qJ(6));
t158 = m(7) * pkin(5) + t84 * mrSges(7,1) - t81 * mrSges(7,2) + mrSges(6,1);
t152 = -m(7) * pkin(11) + mrSges(6,2) - mrSges(7,3);
t149 = m(6) + m(7);
t157 = t149 * pkin(4);
t133 = cos(pkin(6));
t78 = sin(pkin(7));
t79 = cos(pkin(7));
t128 = pkin(6) + qJ(2);
t111 = sin(t128) / 0.2e1;
t129 = pkin(6) - qJ(2);
t120 = sin(t129);
t98 = t111 + t120 / 0.2e1;
t56 = t133 * t79 - t98 * t78;
t82 = sin(qJ(4));
t85 = cos(qJ(4));
t126 = pkin(7) + qJ(3);
t110 = sin(t126) / 0.2e1;
t127 = pkin(7) - qJ(3);
t119 = sin(t127);
t66 = t110 - t119 / 0.2e1;
t112 = cos(t126) / 0.2e1;
t121 = cos(t127);
t68 = t112 - t121 / 0.2e1;
t113 = cos(t128) / 0.2e1;
t122 = cos(t129);
t69 = t113 - t122 / 0.2e1;
t86 = cos(qJ(3));
t92 = -t133 * t68 + t98 * t66 - t69 * t86;
t156 = t56 * t85 - t92 * t82;
t130 = sin(pkin(12));
t131 = sin(pkin(6));
t106 = t131 * t130;
t100 = t122 / 0.2e1 + t113;
t132 = cos(pkin(12));
t148 = sin(qJ(2));
t91 = t130 * t100 + t132 * t148;
t43 = t79 * t106 + t91 * t78;
t67 = t111 - t120 / 0.2e1;
t87 = cos(qJ(2));
t59 = -t130 * t67 + t132 * t87;
t88 = -t68 * t106 + t59 * t86 - t91 * t66;
t155 = t43 * t85 - t88 * t82;
t107 = t132 * t131;
t90 = -t132 * t100 + t130 * t148;
t42 = -t79 * t107 + t90 * t78;
t57 = t130 * t87 + t132 * t67;
t89 = t68 * t107 + t57 * t86 - t90 * t66;
t154 = t42 * t85 - t89 * t82;
t125 = m(4) + m(5) + t149;
t153 = pkin(2) * t125 + mrSges(3,1);
t102 = -m(5) * pkin(3) - t85 * mrSges(5,1) + t82 * mrSges(5,2) - mrSges(4,1);
t77 = qJ(4) + pkin(13);
t75 = sin(t77);
t76 = cos(t77);
t151 = -t152 * t75 + t158 * t76 - t102;
t96 = -m(5) * pkin(10) - t81 * mrSges(7,1) - t84 * mrSges(7,2) + mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t141 = t75 * t78;
t140 = t76 * t78;
t99 = t121 / 0.2e1 + t112;
t97 = t110 + t119 / 0.2e1;
t95 = t97 * t131;
t93 = -mrSges(3,2) + (t85 * mrSges(5,2) + mrSges(4,3) + (mrSges(5,1) + t157) * t82 + t125 * pkin(9)) * t78;
t83 = sin(qJ(3));
t80 = -qJ(5) - pkin(10);
t74 = pkin(4) * t85 + pkin(3);
t41 = t69 * t66 + t98 * t86;
t40 = -t69 * t99 + t98 * t83;
t33 = -t133 * t97 - t69 * t83 - t98 * t99;
t31 = -t59 * t66 - t91 * t86;
t30 = t59 * t99 - t91 * t83;
t29 = -t57 * t66 - t90 * t86;
t28 = t57 * t99 - t90 * t83;
t18 = -t130 * t95 + t59 * t83 + t91 * t99;
t15 = t132 * t95 + t57 * t83 + t90 * t99;
t10 = t56 * t75 + t76 * t92;
t4 = t43 * t75 + t76 * t88;
t2 = t42 * t75 + t76 * t89;
t1 = [(-m(2) - m(3) - t125) * g(3) (t152 * (t69 * t140 + t41 * t75) + t102 * t41 - t158 * (-t69 * t141 + t41 * t76) + t96 * t40 + t93 * t69 - t153 * t98 - t149 * (-t40 * t80 + t41 * t74)) * g(3) + (t152 * (-t140 * t57 + t29 * t75) - t158 * (t141 * t57 + t29 * t76) + t102 * t29 + t96 * t28 - t93 * t57 + t153 * t90 - t149 * (-t28 * t80 + t29 * t74)) * g(2) + (t152 * (-t140 * t59 + t31 * t75) - t158 * (t141 * t59 + t31 * t76) + t102 * t31 + t96 * t30 - t93 * t59 + t153 * t91 - t149 * (-t30 * t80 + t31 * t74)) * g(1) (-t149 * (-t33 * t74 - t80 * t92) + t96 * t92 + t151 * t33) * g(3) + (-t149 * (-t15 * t74 - t80 * t89) + t96 * t89 + t151 * t15) * g(2) + (-t149 * (-t18 * t74 - t80 * t88) + t96 * t88 + t151 * t18) * g(1) (-t156 * mrSges(5,1) - (-t56 * t82 - t85 * t92) * mrSges(5,2) - t158 * (t56 * t76 - t75 * t92) + t152 * t10) * g(3) + (-t154 * mrSges(5,1) - (-t42 * t82 - t85 * t89) * mrSges(5,2) + t152 * t2 - t158 * (t42 * t76 - t75 * t89)) * g(2) + (-t155 * mrSges(5,1) - (-t43 * t82 - t85 * t88) * mrSges(5,2) + t152 * t4 - t158 * (t43 * t76 - t75 * t88)) * g(1) + (-g(1) * t155 - g(2) * t154 - g(3) * t156) * t157, t149 * (-g(1) * t18 - g(2) * t15 - g(3) * t33) -g(1) * ((t18 * t84 - t4 * t81) * mrSges(7,1) + (-t18 * t81 - t4 * t84) * mrSges(7,2)) - g(2) * ((t15 * t84 - t2 * t81) * mrSges(7,1) + (-t15 * t81 - t2 * t84) * mrSges(7,2)) - g(3) * ((-t10 * t81 + t33 * t84) * mrSges(7,1) + (-t10 * t84 - t33 * t81) * mrSges(7,2))];
taug  = t1(:);
