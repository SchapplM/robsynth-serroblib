% Calculate Gravitation load on the joints for
% S6PRRRPP3
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
% Datum: 2019-03-08 22:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRRPP3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPP3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPP3_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPP3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:54:28
% EndTime: 2019-03-08 22:54:29
% DurationCPUTime: 0.86s
% Computational Cost: add. (591->117), mult. (1511->165), div. (0->0), fcn. (1811->10), ass. (0->76)
t59 = sin(qJ(4));
t62 = cos(qJ(4));
t123 = pkin(4) * t62 + qJ(5) * t59;
t71 = m(7) * qJ(6) + mrSges(5,1) - mrSges(6,2) + mrSges(7,3);
t122 = -mrSges(7,1) - mrSges(6,1) - mrSges(5,3);
t113 = mrSges(5,2) - mrSges(7,2) - mrSges(6,3);
t121 = m(6) + m(7);
t60 = sin(qJ(3));
t63 = cos(qJ(3));
t120 = -t63 * mrSges(4,1) + t60 * mrSges(4,2) - mrSges(3,1);
t118 = mrSges(3,2) - mrSges(4,3);
t57 = sin(pkin(10));
t61 = sin(qJ(2));
t64 = cos(qJ(2));
t91 = cos(pkin(10));
t92 = cos(pkin(6));
t72 = t92 * t91;
t43 = t57 * t64 + t61 * t72;
t58 = sin(pkin(6));
t80 = t58 * t91;
t20 = -t43 * t60 - t63 * t80;
t117 = t123 * t20;
t102 = t58 * t63;
t81 = t57 * t92;
t45 = -t61 * t81 + t64 * t91;
t22 = t102 * t57 - t45 * t60;
t116 = t123 * t22;
t103 = t58 * t61;
t46 = -t103 * t60 + t63 * t92;
t115 = t123 * t46;
t114 = t121 * pkin(4) + t71;
t112 = t113 * t59 - t71 * t62 - mrSges(4,1);
t111 = mrSges(4,2) - m(7) * (pkin(5) + pkin(9)) + t122;
t110 = -t122 * t60 - t120;
t107 = pkin(3) * t63;
t42 = t57 * t61 - t64 * t72;
t105 = t42 * t60;
t44 = t61 * t91 + t64 * t81;
t104 = t44 * t60;
t101 = t58 * t64;
t100 = t59 * t63;
t96 = t62 * t63;
t95 = t63 * t64;
t94 = pkin(2) * t101 + pkin(8) * t103;
t90 = t60 * t101;
t89 = t58 * t95;
t87 = -t42 * pkin(2) + t43 * pkin(8);
t86 = -t44 * pkin(2) + pkin(8) * t45;
t17 = t20 * pkin(3);
t21 = t43 * t63 - t60 * t80;
t85 = pkin(9) * t21 + t17;
t18 = t22 * pkin(3);
t23 = t57 * t58 * t60 + t45 * t63;
t84 = pkin(9) * t23 + t18;
t41 = t46 * pkin(3);
t47 = t102 * t61 + t60 * t92;
t83 = pkin(9) * t47 + t41;
t79 = pkin(3) * t89 + pkin(9) * t90 + t94;
t74 = -pkin(9) * t105 - t42 * t107 + t87;
t73 = -pkin(9) * t104 - t44 * t107 + t86;
t70 = -t121 * qJ(5) + t113;
t10 = -t42 * t96 + t43 * t59;
t9 = -t100 * t42 - t43 * t62;
t66 = t10 * pkin(4) + qJ(5) * t9 + t74;
t11 = -t100 * t44 - t45 * t62;
t12 = -t44 * t96 + t45 * t59;
t65 = t12 * pkin(4) + qJ(5) * t11 + t73;
t28 = (t59 * t61 + t62 * t95) * t58;
t27 = -t103 * t62 + t59 * t89;
t25 = -t101 * t59 + t47 * t62;
t24 = t101 * t62 + t47 * t59;
t6 = t23 * t62 + t44 * t59;
t5 = t23 * t59 - t44 * t62;
t4 = t21 * t62 + t42 * t59;
t3 = t21 * t59 - t42 * t62;
t1 = [(-m(2) - m(3) - m(4) - m(5) - t121) * g(3) (-m(4) * t94 - m(5) * t79 - t121 * (t28 * pkin(4) + t27 * qJ(5) + t79) + (t118 * t61 + t120 * t64) * t58 + (-m(7) * pkin(5) + t122) * t90 - t71 * t28 + t113 * t27) * g(3) + (-m(4) * t87 - m(5) * t74 - m(6) * t66 - m(7) * (-pkin(5) * t105 + t66) + t118 * t43 + t113 * t9 - t71 * t10 + t110 * t42) * g(2) + (-m(4) * t86 - m(5) * t73 - m(6) * t65 - m(7) * (-pkin(5) * t104 + t65) + t118 * t45 - t71 * t12 + t113 * t11 + t110 * t44) * g(1) (-m(5) * t83 - m(6) * (t83 + t115) - m(7) * (t41 + t115) + t111 * t47 + t112 * t46) * g(3) + (-m(5) * t85 - m(6) * (t85 + t117) - m(7) * (t17 + t117) + t111 * t21 + t112 * t20) * g(2) + (-m(5) * t84 - m(6) * (t84 + t116) - m(7) * (t18 + t116) + t111 * t23 + t112 * t22) * g(1) (t114 * t24 + t25 * t70) * g(3) + (t114 * t3 + t4 * t70) * g(2) + (t114 * t5 + t6 * t70) * g(1), t121 * (-g(1) * t5 - g(2) * t3 - g(3) * t24) (-g(1) * t6 - g(2) * t4 - g(3) * t25) * m(7)];
taug  = t1(:);
