% Calculate Gravitation load on the joints for
% S6RPRRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 07:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRRR11_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR11_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR11_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR11_gravloadJ_floatb_twist_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR11_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR11_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:38:18
% EndTime: 2019-03-09 07:38:21
% DurationCPUTime: 1.32s
% Computational Cost: add. (1315->132), mult. (3414->184), div. (0->0), fcn. (4343->16), ass. (0->70)
t110 = cos(pkin(7));
t115 = cos(qJ(1));
t106 = sin(pkin(13));
t113 = sin(qJ(1));
t109 = cos(pkin(13));
t111 = cos(pkin(6));
t95 = t111 * t109;
t81 = t113 * t106 - t115 * t95;
t107 = sin(pkin(7));
t108 = sin(pkin(6));
t91 = t108 * t107;
t142 = t81 * t110 + t115 * t91;
t114 = cos(qJ(3));
t93 = t111 * t106;
t41 = t113 * t109 + t115 * t93;
t59 = sin(qJ(3));
t26 = -t114 * t41 + t142 * t59;
t58 = sin(qJ(4));
t61 = cos(qJ(4));
t92 = t110 * t108;
t69 = -t81 * t107 + t115 * t92;
t141 = t26 * t58 - t69 * t61;
t140 = t26 * t61 + t69 * t58;
t134 = -m(5) - m(6);
t57 = sin(qJ(5));
t137 = mrSges(4,2) - mrSges(5,3) - m(7) * (pkin(5) * t57 + pkin(10));
t60 = cos(qJ(5));
t52 = pkin(5) * t60 + pkin(4);
t56 = qJ(5) + qJ(6);
t53 = sin(t56);
t54 = cos(t56);
t126 = m(6) * pkin(4) + m(7) * t52 + t60 * mrSges(6,1) + t54 * mrSges(7,1) - t57 * mrSges(6,2) - t53 * mrSges(7,2) + mrSges(5,1);
t120 = mrSges(5,2) + m(7) * (-pkin(12) - pkin(11)) - mrSges(7,3) - m(6) * pkin(11) - mrSges(6,3);
t75 = t106 * t115 + t113 * t95;
t63 = t75 * t107 + t113 * t92;
t127 = m(7) - t134;
t136 = pkin(3) * t127 - t120 * t58 + t126 * t61 + mrSges(4,1);
t135 = -t57 * mrSges(6,1) - t53 * mrSges(7,1) - t60 * mrSges(6,2) - t54 * mrSges(7,2) + t137;
t23 = t142 * t114 + t41 * t59;
t131 = t75 * t110 - t113 * t91;
t130 = t111 * t107 + t109 * t92;
t128 = -m(7) * pkin(5) - mrSges(6,1);
t121 = t134 * pkin(10) + t135;
t118 = (t140 * t53 + t23 * t54) * mrSges(7,1) + (t140 * t54 - t23 * t53) * mrSges(7,2);
t42 = t109 * t115 - t113 * t93;
t28 = t42 * t114 - t131 * t59;
t16 = t28 * t61 + t63 * t58;
t27 = t131 * t114 + t42 * t59;
t5 = -t16 * t53 + t27 * t54;
t6 = t16 * t54 + t27 * t53;
t117 = t5 * mrSges(7,1) - t6 * mrSges(7,2);
t90 = t108 * t106;
t33 = t114 * t90 + t130 * t59;
t40 = -t109 * t91 + t111 * t110;
t22 = t33 * t61 + t40 * t58;
t32 = -t130 * t114 + t59 * t90;
t116 = (-t22 * t53 + t32 * t54) * mrSges(7,1) + (-t22 * t54 - t32 * t53) * mrSges(7,2);
t99 = t108 * t113;
t112 = t115 * pkin(1) + qJ(2) * t99;
t100 = t115 * t108;
t102 = -pkin(1) * t113 + qJ(2) * t100;
t7 = -t16 * t57 + t27 * t60;
t70 = -t41 * pkin(2) + t69 * pkin(9) + t102;
t68 = t26 * pkin(3) + t70;
t67 = t42 * pkin(2) + t63 * pkin(9) + t112;
t66 = t28 * pkin(3) + t67;
t64 = t27 * pkin(10) + t66;
t15 = t28 * t58 - t63 * t61;
t8 = t16 * t60 + t27 * t57;
t1 = [(-mrSges(2,1) * t115 + t113 * mrSges(2,2) - m(3) * t112 - t42 * mrSges(3,1) + t75 * mrSges(3,2) - mrSges(3,3) * t99 - m(4) * t67 - t28 * mrSges(4,1) - t63 * mrSges(4,3) - m(5) * t64 - t16 * mrSges(5,1) - m(6) * (t16 * pkin(4) + t64) - t8 * mrSges(6,1) - t7 * mrSges(6,2) - m(7) * (t16 * t52 + t66) - t6 * mrSges(7,1) - t5 * mrSges(7,2) + t137 * t27 + t120 * t15) * g(2) + (-m(3) * t102 - m(4) * t70 - m(7) * t68 + t113 * mrSges(2,1) + t41 * mrSges(3,1) - t26 * mrSges(4,1) + mrSges(2,2) * t115 - t81 * mrSges(3,2) - mrSges(3,3) * t100 - t69 * mrSges(4,3) + t134 * (-pkin(10) * t23 + t68) - t126 * t140 - t135 * t23 + t120 * t141) * g(1) (-g(1) * t99 + g(2) * t100 - g(3) * t111) * (m(3) + m(4) + t127) (t121 * t33 + t136 * t32) * g(3) + (-t121 * t26 + t136 * t23) * g(2) + (t121 * t28 + t136 * t27) * g(1) (t120 * t22 - t126 * (-t33 * t58 + t40 * t61)) * g(3) + (-t120 * t140 - t126 * t141) * g(2) + (t120 * t16 + t126 * t15) * g(1) (-(-t22 * t60 - t32 * t57) * mrSges(6,2) - t116 + t128 * (-t22 * t57 + t32 * t60)) * g(3) + (-(t140 * t60 - t23 * t57) * mrSges(6,2) - t118 + t128 * (t140 * t57 + t23 * t60)) * g(2) + (t8 * mrSges(6,2) + t128 * t7 - t117) * g(1), -g(1) * t117 - g(2) * t118 - g(3) * t116];
taug  = t1(:);
