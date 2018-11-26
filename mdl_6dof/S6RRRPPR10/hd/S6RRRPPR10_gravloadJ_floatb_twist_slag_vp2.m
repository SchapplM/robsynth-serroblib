% Calculate Gravitation load on the joints for
% S6RRRPPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta5]';
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
% Datum: 2018-11-23 17:40
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRRPPR10_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR10_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR10_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR10_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR10_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR10_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:40:16
% EndTime: 2018-11-23 17:40:17
% DurationCPUTime: 1.26s
% Computational Cost: add. (1467->127), mult. (1771->151), div. (0->0), fcn. (1760->16), ass. (0->75)
t49 = sin(pkin(11));
t51 = cos(pkin(11));
t141 = -t51 * mrSges(6,2) - (m(7) * pkin(5) + mrSges(6,1)) * t49 + mrSges(4,2) - mrSges(5,3);
t129 = m(6) + m(7);
t124 = m(5) + t129;
t139 = t124 * qJ(4);
t65 = m(7) * (-pkin(10) - qJ(5)) - mrSges(7,3) - m(6) * qJ(5) - mrSges(4,1) + mrSges(5,2) - mrSges(6,3);
t48 = pkin(11) + qJ(6);
t45 = sin(t48);
t46 = cos(t48);
t73 = t45 * mrSges(7,1) + t46 * mrSges(7,2);
t137 = -t141 + t73;
t53 = sin(qJ(3));
t56 = cos(qJ(3));
t136 = pkin(3) * t56 + qJ(4) * t53;
t114 = t137 * t53 - t65 * t56 + mrSges(3,1);
t119 = -t51 * mrSges(6,1) + t49 * mrSges(6,2) - mrSges(5,1) + mrSges(3,2) - mrSges(4,3);
t44 = pkin(5) * t51 + pkin(4);
t74 = t46 * mrSges(7,1) - t45 * mrSges(7,2);
t116 = -m(6) * (pkin(4) + pkin(9)) + t119 - m(7) * (pkin(9) + t44) - t74;
t131 = pkin(3) * t124 - t65;
t109 = cos(qJ(1));
t54 = sin(qJ(2));
t55 = sin(qJ(1));
t96 = pkin(6) + qJ(2);
t78 = cos(t96) / 0.2e1;
t97 = pkin(6) - qJ(2);
t82 = cos(t97);
t63 = t82 / 0.2e1 + t78;
t24 = -t109 * t63 + t54 * t55;
t127 = t136 * t24;
t27 = t109 * t54 + t55 * t63;
t126 = t136 * t27;
t80 = sin(t96);
t76 = t80 / 0.2e1;
t81 = sin(t97);
t77 = t81 / 0.2e1;
t35 = t76 + t77;
t125 = t136 * t35;
t123 = -m(6) * pkin(4) - m(7) * t44 - pkin(9) * (m(4) + t124) + t119;
t120 = -t137 - t139;
t50 = sin(pkin(6));
t107 = t50 * t55;
t57 = cos(qJ(2));
t104 = t55 * t57;
t102 = t109 * pkin(1) + pkin(8) * t107;
t99 = cos(pkin(6));
t70 = t76 - t81 / 0.2e1;
t89 = t109 * t57;
t28 = -t55 * t70 + t89;
t93 = t28 * pkin(2) + t102;
t90 = t50 * t109;
t87 = -pkin(1) * t55 + pkin(8) * t90;
t18 = t24 * pkin(2);
t62 = t77 - t80 / 0.2e1;
t26 = -t109 * t62 + t104;
t86 = pkin(9) * t26 - t18;
t20 = t27 * pkin(2);
t29 = t55 * t62 + t89;
t85 = pkin(9) * t29 - t20;
t34 = t35 * pkin(2);
t36 = t78 - t82 / 0.2e1;
t84 = -pkin(9) * t36 + t34;
t25 = t109 * t70 + t104;
t8 = t25 * t56 - t53 * t90;
t79 = -t25 * pkin(2) + t87;
t7 = t25 * t53 + t56 * t90;
t59 = -t139 + t141;
t23 = -t36 * t56 + t53 * t99;
t22 = -t36 * t53 - t56 * t99;
t12 = t107 * t53 + t28 * t56;
t11 = -t107 * t56 + t28 * t53;
t2 = t11 * t45 + t27 * t46;
t1 = t11 * t46 - t27 * t45;
t3 = [(-t109 * mrSges(2,1) - m(3) * t102 - t28 * mrSges(3,1) - m(4) * t93 - t2 * mrSges(7,1) - t1 * mrSges(7,2) + (-mrSges(3,3) * t50 + mrSges(2,2)) * t55 + t65 * t12 + t59 * t11 + t123 * t27 - t124 * (t12 * pkin(3) + t93)) * g(2) + (t55 * mrSges(2,1) + t109 * mrSges(2,2) - m(3) * t87 + t25 * mrSges(3,1) - mrSges(3,3) * t90 - m(4) * t79 - t65 * t8 - (t59 - t73) * t7 + (-t123 + t74) * t24 + t124 * (pkin(3) * t8 - t79)) * g(1) (-m(4) * t84 - m(5) * (t84 + t125) - t129 * (t34 + t125) - t116 * t36 - t114 * t35) * g(3) + (-m(4) * t86 - m(5) * (t86 - t127) - t129 * (-t18 - t127) + t116 * t26 + t114 * t24) * g(2) + (-m(4) * t85 - m(5) * (t85 - t126) - t129 * (-t20 - t126) + t116 * t29 + t114 * t27) * g(1) (t120 * t23 + t131 * t22) * g(3) + (t120 * t8 + t131 * t7) * g(2) + (t131 * t11 + t120 * t12) * g(1), t124 * (-g(1) * t11 - g(2) * t7 - g(3) * t22) t129 * (-g(1) * t12 - g(2) * t8 - g(3) * t23) -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * ((-t24 * t45 + t46 * t7) * mrSges(7,1) + (-t24 * t46 - t45 * t7) * mrSges(7,2)) - g(3) * ((t22 * t46 + t35 * t45) * mrSges(7,1) + (-t22 * t45 + t35 * t46) * mrSges(7,2))];
taug  = t3(:);
