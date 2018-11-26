% Calculate Gravitation load on the joints for
% S6RRRRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2018-11-23 18:43
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRRRRR7_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR7_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR7_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR7_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR7_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR7_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:42:42
% EndTime: 2018-11-23 18:42:44
% DurationCPUTime: 1.51s
% Computational Cost: add. (1930->156), mult. (2130->196), div. (0->0), fcn. (2146->18), ass. (0->78)
t63 = qJ(4) + qJ(5);
t58 = cos(t63);
t69 = cos(qJ(4));
t60 = t69 * pkin(4);
t44 = pkin(5) * t58 + t60;
t40 = pkin(3) + t44;
t59 = qJ(6) + t63;
t54 = sin(t59);
t55 = cos(t59);
t56 = t60 + pkin(3);
t57 = sin(t63);
t65 = sin(qJ(4));
t132 = m(5) * pkin(3) + m(6) * t56 + m(7) * t40 + t69 * mrSges(5,1) + t58 * mrSges(6,1) + t55 * mrSges(7,1) - t65 * mrSges(5,2) - t57 * mrSges(6,2) - t54 * mrSges(7,2) + mrSges(4,1);
t72 = -pkin(11) - pkin(10);
t76 = mrSges(4,2) + m(7) * (-pkin(12) + t72) - mrSges(7,3) + m(6) * t72 - mrSges(6,3) - m(5) * pkin(10) - mrSges(5,3);
t147 = -m(4) - m(5);
t109 = m(6) + m(7) - t147;
t66 = sin(qJ(3));
t70 = cos(qJ(3));
t142 = pkin(2) * t109 + t132 * t70 - t76 * t66 + mrSges(3,1);
t143 = -t57 * mrSges(6,1) - t54 * mrSges(7,1) - t69 * mrSges(5,2) - t58 * mrSges(6,2) - t55 * mrSges(7,2);
t122 = pkin(4) * t65;
t140 = mrSges(3,2) - mrSges(4,3);
t43 = pkin(5) * t57 + t122;
t128 = -t65 * mrSges(5,1) + t140 + t147 * pkin(9) - m(6) * (pkin(9) + t122) - m(7) * (pkin(9) + t43) + t143;
t139 = m(6) * pkin(4) + mrSges(5,1);
t112 = pkin(6) - qJ(2);
t101 = cos(t112);
t111 = pkin(6) + qJ(2);
t95 = cos(t111) / 0.2e1;
t42 = t95 - t101 / 0.2e1;
t64 = cos(pkin(6));
t30 = -t42 * t70 + t64 * t66;
t99 = sin(t111);
t93 = t99 / 0.2e1;
t100 = sin(t112);
t94 = t100 / 0.2e1;
t41 = t93 + t94;
t116 = (-t30 * t54 - t41 * t55) * mrSges(7,1) + (-t30 * t55 + t41 * t54) * mrSges(7,2);
t87 = -t30 * t57 - t41 * t58;
t138 = -t87 * mrSges(6,1) - (-t30 * t58 + t41 * t57) * mrSges(6,2) - t116;
t113 = sin(pkin(6));
t68 = sin(qJ(1));
t102 = t68 * t113;
t119 = cos(qJ(1));
t71 = cos(qJ(2));
t106 = t119 * t71;
t84 = t93 - t100 / 0.2e1;
t35 = -t68 * t84 + t106;
t24 = t66 * t102 + t35 * t70;
t67 = sin(qJ(2));
t78 = t101 / 0.2e1 + t95;
t34 = t119 * t67 + t68 * t78;
t11 = -t24 * t57 + t34 * t58;
t12 = t24 * t58 + t34 * t57;
t5 = -t24 * t54 + t34 * t55;
t6 = t24 * t55 + t34 * t54;
t124 = t5 * mrSges(7,1) - t6 * mrSges(7,2);
t137 = -t11 * mrSges(6,1) + t12 * mrSges(6,2) - t124;
t118 = t68 * t71;
t32 = t119 * t84 + t118;
t96 = t119 * t113;
t20 = t32 * t70 - t66 * t96;
t31 = -t119 * t78 + t67 * t68;
t125 = (-t20 * t54 + t31 * t55) * mrSges(7,1) + (-t20 * t55 - t31 * t54) * mrSges(7,2);
t89 = -t20 * t57 + t31 * t58;
t136 = -t89 * mrSges(6,1) - (-t20 * t58 - t31 * t57) * mrSges(6,2) - t125;
t133 = m(7) * t43 + t109 * pkin(9) - t140;
t126 = m(7) * pkin(5);
t114 = t119 * pkin(1) + pkin(8) * t102;
t110 = t35 * pkin(2) + t114;
t103 = -t68 * pkin(1) + pkin(8) * t96;
t19 = -t32 * t66 - t70 * t96;
t13 = -t24 * t65 + t34 * t69;
t85 = t94 - t99 / 0.2e1;
t23 = -t70 * t102 + t35 * t66;
t14 = t24 * t69 + t34 * t65;
t1 = [(-t119 * mrSges(2,1) - m(3) * t114 - t35 * mrSges(3,1) - m(4) * t110 - t24 * mrSges(4,1) - m(5) * (pkin(3) * t24 + t110) - t14 * mrSges(5,1) - t13 * mrSges(5,2) - m(6) * (t24 * t56 + t110) - t12 * mrSges(6,1) - t11 * mrSges(6,2) - m(7) * (t24 * t40 + t110) - t6 * mrSges(7,1) - t5 * mrSges(7,2) + (-t113 * mrSges(3,3) + mrSges(2,2)) * t68 + (-m(6) * t122 - t133) * t34 + t76 * t23) * g(2) + (t68 * mrSges(2,1) + t119 * mrSges(2,2) - m(3) * t103 + t32 * mrSges(3,1) - mrSges(3,3) * t96 + t132 * t20 + (t139 * t65 + t133 - t143) * t31 + t76 * t19 + t109 * (t32 * pkin(2) - t103)) * g(1) (-t128 * t42 - t142 * t41) * g(3) + (t128 * (-t119 * t85 + t118) + t142 * t31) * g(2) + (t128 * (t68 * t85 + t106) + t142 * t34) * g(1) (t76 * t30 - t132 * (t42 * t66 + t64 * t70)) * g(3) + (-t132 * t19 + t76 * t20) * g(2) + (t132 * t23 + t76 * t24) * g(1) (-(-t30 * t69 + t41 * t65) * mrSges(5,2) - m(7) * (-t30 * t43 - t41 * t44) - t139 * (-t30 * t65 - t41 * t69) + t138) * g(3) + (-(-t20 * t69 - t31 * t65) * mrSges(5,2) - m(7) * (-t20 * t43 + t31 * t44) - t139 * (-t20 * t65 + t31 * t69) + t136) * g(2) + (t14 * mrSges(5,2) - m(7) * (-t24 * t43 + t34 * t44) - t139 * t13 + t137) * g(1) (-t87 * t126 + t138) * g(3) + (-t89 * t126 + t136) * g(2) + (-t11 * t126 + t137) * g(1), -g(1) * t124 - g(2) * t125 - g(3) * t116];
taug  = t1(:);
