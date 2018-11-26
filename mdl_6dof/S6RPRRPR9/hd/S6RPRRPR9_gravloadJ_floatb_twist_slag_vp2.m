% Calculate Gravitation load on the joints for
% S6RPRRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
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
% Datum: 2018-11-23 16:21
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RPRRPR9_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR9_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR9_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR9_gravloadJ_floatb_twist_slag_vp2: pkin has to be [13x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR9_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR9_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:20:45
% EndTime: 2018-11-23 16:20:46
% DurationCPUTime: 1.51s
% Computational Cost: add. (3207->154), mult. (3330->206), div. (0->0), fcn. (3289->24), ass. (0->81)
t76 = sin(pkin(6));
t86 = cos(qJ(1));
t141 = t76 * t86;
t75 = sin(pkin(7));
t78 = cos(pkin(7));
t135 = sin(pkin(12));
t150 = sin(qJ(1));
t131 = pkin(6) - pkin(12);
t111 = cos(t131) / 0.2e1;
t130 = pkin(6) + pkin(12);
t123 = cos(t130);
t99 = t111 + t123 / 0.2e1;
t96 = t150 * t135 - t86 * t99;
t41 = -t78 * t141 + t96 * t75;
t81 = sin(qJ(4));
t146 = t41 * t81;
t133 = pkin(7) + qJ(3);
t115 = sin(t133) / 0.2e1;
t134 = pkin(7) - qJ(3);
t125 = sin(t134);
t155 = t115 - t125 / 0.2e1;
t110 = sin(t130) / 0.2e1;
t122 = sin(t131);
t60 = t110 - t122 / 0.2e1;
t77 = cos(pkin(12));
t54 = t150 * t77 + t86 * t60;
t116 = cos(t133) / 0.2e1;
t126 = cos(t134);
t62 = t116 - t126 / 0.2e1;
t85 = cos(qJ(3));
t25 = -t62 * t141 + t155 * t96 - t54 * t85;
t84 = cos(qJ(4));
t168 = -t25 * t84 + t146;
t161 = t25 * t81 + t41 * t84;
t74 = qJ(4) + pkin(13);
t71 = sin(t74);
t72 = cos(t74);
t4 = -t25 * t72 + t41 * t71;
t167 = t25 * t71 + t41 * t72;
t80 = sin(qJ(6));
t83 = cos(qJ(6));
t157 = m(7) * pkin(5) + t83 * mrSges(7,1) - t80 * mrSges(7,2) + mrSges(6,1);
t156 = -m(7) * pkin(11) + mrSges(6,2) - mrSges(7,3);
t127 = t76 * t150;
t92 = t86 * t135 + t150 * t99;
t42 = t78 * t127 + t92 * t75;
t136 = cos(pkin(6));
t61 = t111 - t123 / 0.2e1;
t98 = t110 + t122 / 0.2e1;
t159 = -t136 * t62 + t155 * t98 + t61 * t85;
t53 = t136 * t78 - t98 * t75;
t162 = -t159 * t81 + t53 * t84;
t55 = -t150 * t60 + t86 * t77;
t160 = -t62 * t127 - t155 * t92 + t55 * t85;
t9 = -t160 * t81 + t42 * t84;
t158 = m(6) + m(7);
t112 = -m(5) * pkin(10) + mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t152 = -mrSges(7,1) * t80 - mrSges(7,2) * t83 + t112;
t102 = t126 / 0.2e1 + t116;
t82 = sin(qJ(3));
t100 = t115 + t125 / 0.2e1;
t97 = t76 * t100;
t21 = t96 * t102 + t54 * t82 + t86 * t97;
t153 = m(5) * pkin(3) + mrSges(5,1) * t84 - mrSges(5,2) * t81 - t156 * t71 + t157 * t72 + mrSges(4,1);
t144 = t42 * t81;
t137 = t86 * pkin(1) + qJ(2) * t127;
t117 = -t150 * pkin(1) + qJ(2) * t141;
t90 = -t54 * pkin(2) - t41 * pkin(9) + t117;
t89 = t55 * pkin(2) + t42 * pkin(9) + t137;
t26 = t92 * t102 - t150 * t97 + t55 * t82;
t70 = pkin(4) * t84 + pkin(3);
t79 = -qJ(5) - pkin(10);
t87 = pkin(4) * t144 + t160 * t70 - t26 * t79 + t89;
t31 = -t136 * t100 - t98 * t102 + t61 * t82;
t12 = t159 * t72 + t53 * t71;
t10 = t160 * t84 + t144;
t8 = t160 * t72 + t42 * t71;
t7 = t160 * t71 - t42 * t72;
t2 = t26 * t80 + t8 * t83;
t1 = t26 * t83 - t8 * t80;
t3 = [(-t86 * mrSges(2,1) + t150 * mrSges(2,2) - m(3) * t137 - t55 * mrSges(3,1) + t92 * mrSges(3,2) - mrSges(3,3) * t127 - m(4) * t89 - t160 * mrSges(4,1) - t42 * mrSges(4,3) - m(5) * (pkin(3) * t160 + t89) - t10 * mrSges(5,1) - t9 * mrSges(5,2) - m(6) * t87 - t8 * mrSges(6,1) - m(7) * (t8 * pkin(5) + t87) - t2 * mrSges(7,1) - t1 * mrSges(7,2) + t156 * t7 + t112 * t26) * g(2) + (t150 * mrSges(2,1) + t86 * mrSges(2,2) - m(3) * t117 + t54 * mrSges(3,1) - t96 * mrSges(3,2) - mrSges(3,3) * t141 - m(4) * t90 - t25 * mrSges(4,1) + t41 * mrSges(4,3) - m(5) * (t25 * pkin(3) + t90) + t168 * mrSges(5,1) + t161 * mrSges(5,2) + t156 * t167 + t157 * t4 - t152 * t21 - t158 * (-pkin(4) * t146 + t21 * t79 + t25 * t70 + t90)) * g(1) (-t136 * g(3) + (-t150 * g(1) + t86 * g(2)) * t76) * (m(3) + m(4) + m(5) + t158) (-t158 * (-t159 * t79 - t31 * t70) + t152 * t159 + t153 * t31) * g(3) + (-t158 * (-t21 * t70 + t25 * t79) - t152 * t25 + t153 * t21) * g(2) + (-t158 * (-t160 * t79 - t26 * t70) + t152 * t160 + t153 * t26) * g(1) (-t162 * mrSges(5,1) - (-t159 * t84 - t53 * t81) * mrSges(5,2) + t156 * t12 - t157 * (-t159 * t71 + t53 * t72)) * g(3) + (-t161 * mrSges(5,1) + t168 * mrSges(5,2) + t156 * t4 - t157 * t167) * g(2) + (-t9 * mrSges(5,1) + t10 * mrSges(5,2) + t156 * t8 + t157 * t7) * g(1) + (-g(1) * t9 - g(2) * t161 - g(3) * t162) * t158 * pkin(4), t158 * (-g(1) * t26 - g(2) * t21 - g(3) * t31) -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * ((t21 * t83 - t4 * t80) * mrSges(7,1) + (-t21 * t80 - t4 * t83) * mrSges(7,2)) - g(3) * ((-t12 * t80 + t31 * t83) * mrSges(7,1) + (-t12 * t83 - t31 * t80) * mrSges(7,2))];
taug  = t3(:);
