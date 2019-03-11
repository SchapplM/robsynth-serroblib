% Calculate Gravitation load on the joints for
% S6PRRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1]';
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
% Datum: 2019-03-08 22:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRRPP2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPP2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPP2_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPP2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:48:46
% EndTime: 2019-03-08 22:48:49
% DurationCPUTime: 0.85s
% Computational Cost: add. (586->112), mult. (1498->158), div. (0->0), fcn. (1793->10), ass. (0->70)
t59 = sin(qJ(4));
t62 = cos(qJ(4));
t126 = pkin(4) * t62 + qJ(5) * t59;
t123 = mrSges(7,3) - mrSges(6,2) - mrSges(5,3);
t125 = m(7) * qJ(6) + t123;
t74 = m(7) * pkin(5) + mrSges(5,1) + mrSges(6,1) + mrSges(7,1);
t113 = mrSges(5,2) - mrSges(7,2) - mrSges(6,3);
t60 = sin(qJ(3));
t63 = cos(qJ(3));
t122 = -pkin(3) * t63 - pkin(9) * t60;
t121 = m(6) + m(7);
t120 = -t63 * mrSges(4,1) + t60 * mrSges(4,2) - mrSges(3,1);
t119 = mrSges(3,2) - mrSges(4,3);
t57 = sin(pkin(10));
t61 = sin(qJ(2));
t64 = cos(qJ(2));
t91 = cos(pkin(10));
t92 = cos(pkin(6));
t71 = t92 * t91;
t43 = t57 * t64 + t61 * t71;
t58 = sin(pkin(6));
t80 = t58 * t91;
t20 = -t43 * t60 - t63 * t80;
t117 = t126 * t20;
t103 = t58 * t63;
t81 = t57 * t92;
t45 = -t61 * t81 + t64 * t91;
t22 = t103 * t57 - t45 * t60;
t116 = t126 * t22;
t104 = t58 * t61;
t46 = -t104 * t60 + t63 * t92;
t115 = t126 * t46;
t114 = t121 * pkin(4) + t74;
t112 = t113 * t59 - t74 * t62 - mrSges(4,1);
t111 = mrSges(4,2) - m(7) * (pkin(9) - qJ(6)) + t123;
t109 = -t125 * t60 - t120;
t102 = t58 * t64;
t101 = t59 * t63;
t97 = t62 * t63;
t96 = t63 * t64;
t95 = pkin(2) * t102 + pkin(8) * t104;
t90 = t60 * t102;
t89 = t58 * t96;
t42 = t57 * t61 - t64 * t71;
t86 = -t42 * pkin(2) + pkin(8) * t43;
t44 = t61 * t91 + t64 * t81;
t85 = -t44 * pkin(2) + pkin(8) * t45;
t17 = t20 * pkin(3);
t21 = t43 * t63 - t60 * t80;
t84 = pkin(9) * t21 + t17;
t18 = t22 * pkin(3);
t23 = t57 * t58 * t60 + t45 * t63;
t83 = pkin(9) * t23 + t18;
t41 = t46 * pkin(3);
t47 = t103 * t61 + t60 * t92;
t82 = pkin(9) * t47 + t41;
t79 = pkin(3) * t89 + pkin(9) * t90 + t95;
t73 = t122 * t42 + t86;
t72 = t122 * t44 + t85;
t70 = -t121 * qJ(5) + t113;
t28 = (t59 * t61 + t62 * t96) * t58;
t27 = -t104 * t62 + t59 * t89;
t24 = t102 * t62 + t47 * t59;
t12 = -t44 * t97 + t45 * t59;
t11 = -t101 * t44 - t45 * t62;
t10 = -t42 * t97 + t43 * t59;
t9 = -t101 * t42 - t43 * t62;
t5 = t23 * t59 - t44 * t62;
t3 = t21 * t59 - t42 * t62;
t1 = [(-m(2) - m(3) - m(4) - m(5) - t121) * g(3) (-m(4) * t95 - m(5) * t79 - t121 * (t28 * pkin(4) + t27 * qJ(5) + t79) + (t119 * t61 + t120 * t64) * t58 + t125 * t90 - t74 * t28 + t113 * t27) * g(3) + (-m(4) * t86 - m(5) * t73 - t121 * (t10 * pkin(4) + qJ(5) * t9 + t73) + t119 * t43 + t113 * t9 - t74 * t10 + t109 * t42) * g(2) + (-m(4) * t85 - m(5) * t72 - t121 * (t12 * pkin(4) + qJ(5) * t11 + t72) + t119 * t45 - t74 * t12 + t113 * t11 + t109 * t44) * g(1) (-m(5) * t82 - m(6) * (t82 + t115) - m(7) * (t41 + t115) + t111 * t47 + t112 * t46) * g(3) + (-m(5) * t84 - m(6) * (t84 + t117) - m(7) * (t17 + t117) + t111 * t21 + t112 * t20) * g(2) + (-m(5) * t83 - m(6) * (t83 + t116) - m(7) * (t18 + t116) + t111 * t23 + t112 * t22) * g(1) (t70 * (-t59 * t102 + t47 * t62) + t114 * t24) * g(3) + (t70 * (t21 * t62 + t42 * t59) + t114 * t3) * g(2) + (t70 * (t23 * t62 + t44 * t59) + t114 * t5) * g(1), t121 * (-g(1) * t5 - g(2) * t3 - g(3) * t24) (-g(1) * t22 - g(2) * t20 - g(3) * t46) * m(7)];
taug  = t1(:);
