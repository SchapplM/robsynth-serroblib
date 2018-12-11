% Calculate potential energy for
% S6RRPRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-12-10 18:39
% Revision: bb42a8b95257d9bc83910d26e849f5825122f662 (2018-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function U = S6RRPRRR14_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(14,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRPRRR14_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR14_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_energypot_floatb_twist_slag_vp2: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR14_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR14_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-12-10 18:09:45
% EndTime: 2018-12-10 18:09:46
% DurationCPUTime: 1.22s
% Computational Cost: add. (3330->149), mult. (3365->176), div. (0->0), fcn. (3324->30), ass. (0->89)
t133 = m(6) + m(7);
t135 = t133 + m(4) + m(5);
t74 = pkin(7) + pkin(14);
t63 = sin(t74) / 0.2e1;
t75 = pkin(7) - pkin(14);
t67 = sin(t75);
t50 = t63 + t67 / 0.2e1;
t64 = cos(t75) / 0.2e1;
t68 = cos(t74);
t52 = t64 + t68 / 0.2e1;
t76 = pkin(6) + qJ(2);
t65 = sin(t76) / 0.2e1;
t77 = pkin(6) - qJ(2);
t69 = sin(t77);
t55 = t65 + t69 / 0.2e1;
t66 = cos(t77) / 0.2e1;
t70 = cos(t76);
t59 = t66 - t70 / 0.2e1;
t78 = sin(pkin(14));
t85 = cos(pkin(6));
t32 = t50 * t85 + t52 * t55 - t59 * t78;
t80 = sin(pkin(7));
t84 = cos(pkin(7));
t44 = -t55 * t80 + t84 * t85;
t79 = sin(pkin(8));
t83 = cos(pkin(8));
t23 = -t32 * t79 + t44 * t83;
t81 = sin(pkin(6));
t90 = sin(qJ(1));
t125 = t81 * t90;
t58 = t66 + t70 / 0.2e1;
t89 = sin(qJ(2));
t95 = cos(qJ(1));
t47 = -t90 * t58 - t89 * t95;
t56 = t65 - t69 / 0.2e1;
t94 = cos(qJ(2));
t48 = -t90 * t56 + t94 * t95;
t28 = t125 * t50 + t47 * t52 - t48 * t78;
t117 = t84 * t125;
t39 = -t47 * t80 + t117;
t19 = -t28 * t79 + t39 * t83;
t124 = t81 * t95;
t45 = t58 * t95 - t90 * t89;
t46 = t56 * t95 + t90 * t94;
t26 = -t124 * t50 + t45 * t52 - t46 * t78;
t38 = -t124 * t84 - t45 * t80;
t18 = -t26 * t79 + t38 * t83;
t123 = t95 * pkin(1) + pkin(10) * t125;
t122 = qJ(3) * t80;
t121 = qJ(3) * t84;
t120 = pkin(8) - qJ(4);
t119 = pkin(8) + qJ(4);
t118 = pkin(9) + r_base(3);
t116 = t85 * pkin(10) + t118;
t114 = cos(t119);
t113 = sin(t120);
t112 = -m(7) * pkin(13) + mrSges(6,2) - mrSges(7,3);
t111 = cos(t120) / 0.2e1;
t110 = sin(t119) / 0.2e1;
t109 = -m(1) - m(2) - m(3) - t135;
t72 = t90 * pkin(1);
t108 = t46 * pkin(2) - t122 * t45 + t72;
t107 = t48 * pkin(2) + qJ(3) * t117 - t122 * t47 + t123;
t86 = sin(qJ(6));
t91 = cos(qJ(6));
t106 = -m(7) * pkin(5) - mrSges(7,1) * t91 + mrSges(7,2) * t86 - mrSges(6,1);
t105 = t59 * pkin(2) + t85 * t121 - t122 * t55 + t116;
t104 = t111 + t114 / 0.2e1;
t103 = t110 + t113 / 0.2e1;
t102 = -mrSges(7,1) * t86 - mrSges(7,2) * t91 - pkin(12) * t133 + mrSges(5,2) - mrSges(6,3);
t51 = t63 - t67 / 0.2e1;
t53 = t64 - t68 / 0.2e1;
t82 = cos(pkin(14));
t27 = -t124 * t53 + t45 * t51 + t46 * t82;
t101 = t27 * pkin(3) + t18 * pkin(11) + t108;
t29 = t125 * t53 + t47 * t51 + t48 * t82;
t99 = t29 * pkin(3) + t19 * pkin(11) + t107;
t33 = t51 * t55 + t53 * t85 + t59 * t82;
t97 = t33 * pkin(3) + t23 * pkin(11) + t105;
t93 = cos(qJ(4));
t92 = cos(qJ(5));
t88 = sin(qJ(4));
t87 = sin(qJ(5));
t57 = t111 - t114 / 0.2e1;
t54 = t110 - t113 / 0.2e1;
t15 = t32 * t54 + t33 * t93 + t44 * t57;
t12 = t28 * t54 + t29 * t93 + t39 * t57;
t10 = t26 * t54 + t27 * t93 + t38 * t57;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t118 - mrSges(2,3) - m(3) * t116 - t59 * mrSges(3,1) - t55 * mrSges(3,2) - t85 * mrSges(3,3) - m(4) * t105 - t33 * mrSges(4,1) - t32 * mrSges(4,2) - t44 * mrSges(4,3) - m(5) * t97 - t15 * mrSges(5,1) - t23 * mrSges(5,3) + t112 * (t15 * t87 - t23 * t92) + t106 * (t15 * t92 + t23 * t87) + t102 * (-t103 * t44 - t104 * t32 + t33 * t88) + t133 * (-t15 * pkin(4) - t97)) * g(3) + (-t90 * mrSges(2,1) - m(3) * t72 - m(4) * t108 - m(5) * t101 - t45 * mrSges(3,2) - t46 * mrSges(3,1) - t38 * mrSges(4,3) - t26 * mrSges(4,2) - t27 * mrSges(4,1) - t18 * mrSges(5,3) - t10 * mrSges(5,1) - mrSges(1,2) + t109 * r_base(2) + t106 * (t10 * t92 + t18 * t87) + t102 * (-t38 * t103 - t26 * t104 + t27 * t88) + t112 * (t10 * t87 - t18 * t92) + (-mrSges(2,2) + (m(3) * pkin(10) + mrSges(3,3) + t135 * (pkin(10) + t121)) * t81) * t95 + t133 * (-t10 * pkin(4) - t101)) * g(2) + (-t95 * mrSges(2,1) - m(4) * t107 - m(5) * t99 - m(3) * t123 - t47 * mrSges(3,2) - t48 * mrSges(3,1) - t39 * mrSges(4,3) - t28 * mrSges(4,2) - t29 * mrSges(4,1) - t12 * mrSges(5,1) - t19 * mrSges(5,3) - mrSges(1,1) + t109 * r_base(1) + (-t81 * mrSges(3,3) + mrSges(2,2)) * t90 + t112 * (t12 * t87 - t19 * t92) + t106 * (t12 * t92 + t19 * t87) + t102 * (-t103 * t39 - t28 * t104 + t29 * t88) + t133 * (-t12 * pkin(4) - t99)) * g(1);
U  = t1;
