% Calculate Gravitation load on the joints for
% S6RRRPRR8
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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPRR8_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR8_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR8_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR8_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR8_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR8_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:44:09
% EndTime: 2019-03-09 18:44:13
% DurationCPUTime: 1.58s
% Computational Cost: add. (869->142), mult. (1456->193), div. (0->0), fcn. (1686->14), ass. (0->63)
t61 = cos(qJ(5));
t46 = pkin(5) * t61 + pkin(4);
t54 = qJ(5) + qJ(6);
t50 = sin(t54);
t51 = cos(t54);
t57 = sin(qJ(5));
t121 = m(6) * pkin(4) + m(7) * t46 + t61 * mrSges(6,1) + t51 * mrSges(7,1) - t57 * mrSges(6,2) - t50 * mrSges(7,2) + mrSges(5,1);
t120 = mrSges(5,2) + m(7) * (-pkin(11) - pkin(10)) - mrSges(7,3) - m(6) * pkin(10) - mrSges(6,3);
t53 = qJ(3) + pkin(12);
t48 = sin(t53);
t49 = cos(t53);
t58 = sin(qJ(3));
t62 = cos(qJ(3));
t119 = m(4) * pkin(2) + t62 * mrSges(4,1) - t58 * mrSges(4,2) - t120 * t48 + t121 * t49 + mrSges(3,1);
t78 = m(4) * pkin(9) - mrSges(3,2) + mrSges(4,3) + mrSges(5,3);
t136 = t50 * mrSges(7,1) + mrSges(6,2) * t61 + t51 * mrSges(7,2) + t78;
t132 = m(7) * pkin(5);
t98 = t57 * t132;
t118 = -t57 * mrSges(6,1) - t136 - t98;
t122 = m(5) + m(6) + m(7);
t55 = sin(pkin(6));
t59 = sin(qJ(2));
t106 = t55 * t59;
t99 = cos(pkin(6));
t131 = -t58 * t106 + t99 * t62;
t104 = t55 * t62;
t110 = cos(qJ(1));
t63 = cos(qJ(2));
t60 = sin(qJ(1));
t88 = t60 * t99;
t32 = t110 * t63 - t59 * t88;
t17 = t60 * t104 - t32 * t58;
t123 = mrSges(6,1) + t132;
t83 = t99 * t110;
t30 = t59 * t83 + t60 * t63;
t92 = t55 * t110;
t72 = t30 * t58 + t62 * t92;
t12 = t30 * t49 - t48 * t92;
t29 = t59 * t60 - t63 * t83;
t115 = (-t12 * t50 + t29 * t51) * mrSges(7,1) + (-t12 * t51 - t29 * t50) * mrSges(7,2);
t105 = t55 * t60;
t16 = t48 * t105 + t32 * t49;
t31 = t110 * t59 + t63 * t88;
t5 = -t16 * t50 + t31 * t51;
t6 = t16 * t51 + t31 * t50;
t114 = t5 * mrSges(7,1) - t6 * mrSges(7,2);
t103 = t55 * t63;
t24 = t49 * t106 + t99 * t48;
t111 = (-t51 * t103 - t24 * t50) * mrSges(7,1) + (t50 * t103 - t24 * t51) * mrSges(7,2);
t100 = t110 * pkin(1) + pkin(8) * t105;
t97 = t58 * t105;
t90 = -t60 * pkin(1) + pkin(8) * t92;
t11 = -t30 * t48 - t49 * t92;
t41 = t58 * t92;
t89 = -t30 * t62 + t41;
t47 = pkin(3) * t62 + pkin(2);
t56 = -qJ(4) - pkin(9);
t84 = pkin(3) * t97 - t31 * t56 + t32 * t47 + t100;
t7 = -t16 * t57 + t31 * t61;
t18 = t32 * t62 + t97;
t15 = -t49 * t105 + t32 * t48;
t8 = t16 * t61 + t31 * t57;
t1 = [(-t110 * mrSges(2,1) - m(3) * t100 - t32 * mrSges(3,1) - m(4) * (pkin(2) * t32 + t100) - t18 * mrSges(4,1) - t17 * mrSges(4,2) - m(5) * t84 - t16 * mrSges(5,1) - m(6) * (pkin(4) * t16 + t84) - t8 * mrSges(6,1) - t7 * mrSges(6,2) - m(7) * (t16 * t46 + t84) - t6 * mrSges(7,1) - t5 * mrSges(7,2) + (-mrSges(3,3) * t55 + mrSges(2,2)) * t60 + (-t78 - t98) * t31 + t120 * t15) * g(2) + (t60 * mrSges(2,1) + t110 * mrSges(2,2) - m(3) * t90 + t30 * mrSges(3,1) - mrSges(3,3) * t92 - m(4) * (-pkin(2) * t30 + t90) - t89 * mrSges(4,1) - t72 * mrSges(4,2) + t121 * t12 + (t123 * t57 + t136) * t29 + t120 * t11 + t122 * (-pkin(3) * t41 - t29 * t56 + t30 * t47 - t90)) * g(1) (-t122 * (-t29 * t47 - t30 * t56) + t118 * t30 + t119 * t29) * g(2) + (-t122 * (-t31 * t47 - t32 * t56) + t118 * t32 + t119 * t31) * g(1) + (-t122 * t47 * t103 + (-t119 * t63 + (t122 * t56 + t118) * t59) * t55) * g(3) (-t131 * mrSges(4,1) - (-t59 * t104 - t99 * t58) * mrSges(4,2) + t120 * t24 - t121 * (-t48 * t106 + t99 * t49)) * g(3) + (t72 * mrSges(4,1) - t89 * mrSges(4,2) - t121 * t11 + t12 * t120) * g(2) + (-mrSges(4,1) * t17 + mrSges(4,2) * t18 + t120 * t16 + t121 * t15) * g(1) + (-g(1) * t17 + g(2) * t72 - g(3) * t131) * t122 * pkin(3), t122 * (-g(1) * t31 - g(2) * t29 + g(3) * t103) (-(t57 * t103 - t24 * t61) * mrSges(6,2) - t111 - t123 * (-t61 * t103 - t24 * t57)) * g(3) + (-(-t12 * t61 - t29 * t57) * mrSges(6,2) - t115 - t123 * (-t12 * t57 + t29 * t61)) * g(2) + (t8 * mrSges(6,2) - t123 * t7 - t114) * g(1), -g(1) * t114 - g(2) * t115 - g(3) * t111];
taug  = t1(:);
