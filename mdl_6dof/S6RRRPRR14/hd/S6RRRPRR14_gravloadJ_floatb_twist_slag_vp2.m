% Calculate Gravitation load on the joints for
% S6RRRPRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
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
% Datum: 2019-03-09 20:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPRR14_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR14_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR14_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR14_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR14_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR14_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:12:24
% EndTime: 2019-03-09 20:12:27
% DurationCPUTime: 1.32s
% Computational Cost: add. (721->143), mult. (1627->186), div. (0->0), fcn. (1919->12), ass. (0->81)
t147 = m(7) * pkin(5);
t130 = m(5) + m(6) + m(7);
t134 = mrSges(4,2) - mrSges(5,3);
t69 = -qJ(4) * t130 + t134;
t53 = sin(qJ(5));
t91 = t53 * t147;
t146 = t69 - t91;
t145 = m(7) * (-pkin(11) - pkin(10)) - mrSges(7,3) - mrSges(4,1) + mrSges(5,2);
t51 = qJ(5) + qJ(6);
t48 = sin(t51);
t49 = cos(t51);
t57 = cos(qJ(5));
t143 = -t57 * mrSges(6,1) - t49 * mrSges(7,1) + t53 * mrSges(6,2) + t48 * mrSges(7,2);
t141 = t53 * mrSges(6,1);
t100 = t57 * mrSges(6,2);
t74 = t48 * mrSges(7,1) + t49 * mrSges(7,2);
t140 = -t100 - t74;
t129 = mrSges(3,2) - mrSges(5,1) - mrSges(4,3);
t139 = t129 + t143;
t127 = -m(6) * pkin(10) - mrSges(6,3);
t64 = t127 + t145;
t54 = sin(qJ(3));
t58 = cos(qJ(3));
t138 = -mrSges(3,1) + t145 * t58 + (-t74 - t91 + t134) * t54;
t47 = pkin(5) * t57 + pkin(4);
t137 = -m(6) * pkin(4) - m(7) * t47;
t136 = pkin(3) * t130 - t64;
t55 = sin(qJ(2));
t56 = sin(qJ(1));
t59 = cos(qJ(2));
t110 = cos(qJ(1));
t93 = cos(pkin(6));
t77 = t93 * t110;
t32 = t55 * t56 - t59 * t77;
t109 = t32 * t58;
t94 = qJ(4) * t54;
t133 = -pkin(3) * t109 - t32 * t94;
t81 = t56 * t93;
t34 = t110 * t55 + t59 * t81;
t108 = t34 * t58;
t132 = -pkin(3) * t108 - t34 * t94;
t131 = -mrSges(6,1) - t147;
t126 = -pkin(9) * (m(4) + t130) + t129 + t137;
t122 = t140 - t141 + t146;
t119 = t58 * mrSges(6,3) - t138 + (t100 + t141) * t54;
t118 = -m(6) * (pkin(4) + pkin(9)) - m(7) * (pkin(9) + t47) + t139;
t33 = t55 * t77 + t56 * t59;
t52 = sin(pkin(6));
t86 = t52 * t110;
t15 = t33 * t54 + t58 * t86;
t115 = (t15 * t49 - t32 * t48) * mrSges(7,1) + (-t15 * t48 - t32 * t49) * mrSges(7,2);
t105 = t52 * t58;
t35 = t110 * t59 - t55 * t81;
t19 = -t105 * t56 + t35 * t54;
t5 = t19 * t49 - t34 * t48;
t6 = t19 * t48 + t34 * t49;
t114 = t5 * mrSges(7,1) - t6 * mrSges(7,2);
t104 = t52 * t59;
t107 = t52 * t55;
t30 = t107 * t54 - t58 * t93;
t111 = (t104 * t48 + t30 * t49) * mrSges(7,1) + (t104 * t49 - t30 * t48) * mrSges(7,2);
t106 = t52 * t56;
t102 = t53 * t59;
t99 = t57 * t59;
t97 = t58 * t59;
t96 = pkin(2) * t104 + pkin(9) * t107;
t95 = t110 * pkin(1) + pkin(8) * t106;
t26 = t32 * pkin(2);
t90 = -t26 + t133;
t28 = t34 * pkin(2);
t89 = -t28 + t132;
t88 = t35 * pkin(2) + t95;
t84 = -t56 * pkin(1) + pkin(8) * t86;
t83 = t33 * pkin(9) - t26;
t82 = t35 * pkin(9) - t28;
t16 = t33 * t58 - t54 * t86;
t78 = -t33 * pkin(2) + t84;
t7 = t19 * t57 - t34 * t53;
t20 = t106 * t54 + t35 * t58;
t8 = t19 * t53 + t34 * t57;
t1 = [(-t110 * mrSges(2,1) - m(3) * t95 - t35 * mrSges(3,1) - m(4) * t88 - t8 * mrSges(6,1) - t7 * mrSges(6,2) - t6 * mrSges(7,1) - t5 * mrSges(7,2) + (-mrSges(3,3) * t52 + mrSges(2,2)) * t56 + t146 * t19 + t64 * t20 + t126 * t34 - t130 * (t20 * pkin(3) + t88)) * g(2) + (t56 * mrSges(2,1) + t110 * mrSges(2,2) - m(3) * t84 + t33 * mrSges(3,1) - mrSges(3,3) * t86 - m(4) * t78 - t64 * t16 - (t131 * t53 + t140 + t69) * t15 + (-t126 - t143) * t32 + t130 * (pkin(3) * t16 - t78)) * g(1) (-m(4) * t96 - t130 * (t52 * pkin(3) * t97 + t94 * t104 + t96) + (t127 * t97 + (-t102 * mrSges(6,1) - t99 * mrSges(6,2)) * t54 + t138 * t59 + (t137 + t139) * t55) * t52) * g(3) + (-m(4) * t83 - m(5) * (t83 + t133) - m(6) * (-pkin(10) * t109 + t90) - m(7) * t90 + t118 * t33 + t119 * t32) * g(2) + (-m(4) * t82 - m(5) * (t82 + t132) - m(6) * (-pkin(10) * t108 + t89) - m(7) * t89 + t118 * t35 + t119 * t34) * g(1) (t122 * (t105 * t55 + t54 * t93) + t136 * t30) * g(3) + (t122 * t16 + t136 * t15) * g(2) + (t122 * t20 + t136 * t19) * g(1), t130 * (-g(1) * t19 - g(2) * t15 - g(3) * t30) (-(-t30 * t53 + t52 * t99) * mrSges(6,2) - t111 + t131 * (t102 * t52 + t30 * t57)) * g(3) + (-(-t15 * t53 - t32 * t57) * mrSges(6,2) - t115 + t131 * (t15 * t57 - t32 * t53)) * g(2) + (t8 * mrSges(6,2) + t131 * t7 - t114) * g(1), -g(1) * t114 - g(2) * t115 - g(3) * t111];
taug  = t1(:);
