% Calculate Gravitation load on the joints for
% S6RRPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3,theta4]';
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
% Datum: 2018-11-23 16:50
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRPPRR3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPPRR3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR3_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:50:04
% EndTime: 2018-11-23 16:50:05
% DurationCPUTime: 0.92s
% Computational Cost: add. (1464->132), mult. (1164->170), div. (0->0), fcn. (1066->22), ass. (0->74)
t72 = sin(qJ(6));
t75 = cos(qJ(6));
t121 = m(7) * pkin(5) + mrSges(7,1) * t75 - mrSges(7,2) * t72 + mrSges(6,1);
t99 = -m(7) * pkin(10) + mrSges(6,2) - mrSges(7,3);
t123 = -m(4) - m(5);
t122 = m(6) + m(7);
t67 = sin(pkin(12));
t69 = cos(pkin(12));
t81 = m(5) * pkin(3) + t69 * mrSges(5,1) - t67 * mrSges(5,2) + mrSges(4,1);
t91 = -m(5) * qJ(4) + mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t120 = mrSges(7,1) * t72 + mrSges(7,2) * t75 - t91;
t65 = pkin(6) + qJ(2);
t52 = cos(t65) / 0.2e1;
t66 = pkin(6) - qJ(2);
t62 = cos(t66);
t38 = t62 / 0.2e1 + t52;
t119 = m(5) + t122;
t104 = m(4) + t119;
t76 = cos(qJ(2));
t118 = t104 * (pkin(2) * t76 + pkin(1)) + m(3) * pkin(1) + t76 * mrSges(3,1) + mrSges(2,1);
t63 = pkin(12) + qJ(5);
t55 = sin(t63);
t59 = cos(t63);
t117 = t121 * t59 - t55 * t99 + t81;
t115 = sin(t65) / 0.2e1;
t58 = sin(t66);
t113 = pkin(2) * t58;
t68 = sin(pkin(6));
t74 = sin(qJ(1));
t112 = t68 * t74;
t77 = cos(qJ(1));
t111 = t68 * t77;
t73 = sin(qJ(2));
t110 = t73 * t77;
t109 = t74 * t73;
t48 = pkin(2) * t115;
t108 = t113 / 0.2e1 + t48;
t64 = qJ(2) + pkin(11);
t107 = pkin(4) * t67 * t68;
t56 = sin(t64);
t101 = pkin(6) - t64;
t87 = cos(t101) / 0.2e1;
t100 = pkin(6) + t64;
t93 = cos(t100);
t80 = t93 / 0.2e1 + t87;
t22 = t56 * t77 + t74 * t80;
t53 = pkin(4) * t69 + pkin(3);
t71 = -pkin(9) - qJ(4);
t86 = sin(t100) / 0.2e1;
t92 = sin(t101);
t33 = t86 - t92 / 0.2e1;
t60 = cos(t64);
t94 = -t33 * t74 + t60 * t77;
t105 = t107 * t74 - t22 * t71 + t53 * t94;
t95 = t33 * t77 + t60 * t74;
t4 = -t111 * t55 + t59 * t95;
t3 = -t111 * t59 - t55 * t95;
t36 = t38 * pkin(2);
t98 = -pkin(2) * t109 + t36 * t77;
t89 = -pkin(2) * t110 - t74 * t36;
t37 = t115 - t58 / 0.2e1;
t78 = t37 * mrSges(3,1) + mrSges(2,2) + t104 * (-t113 / 0.2e1 + t48 - t68 * (pkin(8) + qJ(3))) + (-m(3) * pkin(8) - mrSges(5,1) * t67 - mrSges(5,2) * t69 - mrSges(3,3) - mrSges(4,3)) * t68;
t70 = cos(pkin(6));
t35 = t87 - t93 / 0.2e1;
t34 = t92 / 0.2e1 + t86;
t26 = -t38 * t74 - t110;
t25 = -t38 * t77 + t109;
t19 = t56 * t74 - t77 * t80;
t18 = t35 * t59 + t55 * t70;
t8 = t112 * t55 + t59 * t94;
t7 = -t112 * t59 + t55 * t94;
t2 = t22 * t72 + t75 * t8;
t1 = t22 * t75 - t72 * t8;
t5 = [(-t26 * mrSges(3,2) - m(6) * t105 - t8 * mrSges(6,1) - m(7) * (pkin(5) * t8 + t105) - t2 * mrSges(7,1) - t1 * mrSges(7,2) + t99 * t7 - t81 * t94 + t78 * t74 + t91 * t22 - t118 * t77) * g(2) + (-t25 * mrSges(3,2) + t99 * t3 + t81 * t95 + t118 * t74 + t78 * t77 + t121 * t4 + t120 * t19 + t122 * (-t107 * t77 - t19 * t71 + t53 * t95)) * g(1) (-(t115 + t58 / 0.2e1) * mrSges(3,1) - (t52 - t62 / 0.2e1) * mrSges(3,2) + t123 * t108 - t122 * (t34 * t53 - t35 * t71 + t108) - t120 * t35 - t117 * t34) * g(3) + (t25 * mrSges(3,1) - (-t37 * t77 - t74 * t76) * mrSges(3,2) + t123 * t98 - t122 * (-t19 * t53 - t71 * t95 + t98) - t120 * t95 + t117 * t19) * g(2) + (-t26 * mrSges(3,1) - (t37 * t74 - t76 * t77) * mrSges(3,2) + t123 * t89 - t122 * (-t22 * t53 - t71 * t94 + t89) - t120 * t94 + t117 * t22) * g(1) (-t70 * g(3) + (-g(1) * t74 + g(2) * t77) * t68) * t104, t119 * (-g(1) * t22 - g(2) * t19 + g(3) * t34) (t99 * t18 - t121 * (-t35 * t55 + t59 * t70)) * g(3) + (-t121 * t3 + t4 * t99) * g(2) + (t121 * t7 + t8 * t99) * g(1), -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * ((t19 * t75 - t4 * t72) * mrSges(7,1) + (-t19 * t72 - t4 * t75) * mrSges(7,2)) - g(3) * ((-t18 * t72 - t34 * t75) * mrSges(7,1) + (-t18 * t75 + t34 * t72) * mrSges(7,2))];
taug  = t5(:);
