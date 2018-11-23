% Calculate potential energy for
% S6RRRRRR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2018-11-23 11:28
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function U = S6RRRRRR10_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(14,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRRRR10_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR10_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRRRRR10_energypot_floatb_twist_slag_vp2: pkin has to be [14x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR10_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR10_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 10:29:11
% EndTime: 2018-11-23 10:29:13
% DurationCPUTime: 1.19s
% Computational Cost: add. (3330->147), mult. (3365->173), div. (0->0), fcn. (3324->30), ass. (0->87)
t135 = m(6) + m(7);
t137 = t135 + m(4) + m(5);
t80 = sin(pkin(6));
t89 = sin(qJ(1));
t124 = t80 * t89;
t77 = pkin(6) - qJ(2);
t66 = cos(t77) / 0.2e1;
t76 = pkin(6) + qJ(2);
t70 = cos(t76);
t58 = t66 + t70 / 0.2e1;
t88 = sin(qJ(2));
t95 = cos(qJ(1));
t47 = -t89 * t58 - t88 * t95;
t79 = sin(pkin(7));
t82 = cos(pkin(7));
t39 = t82 * t124 - t47 * t79;
t64 = sin(t76) / 0.2e1;
t68 = sin(t77);
t53 = t64 + t68 / 0.2e1;
t83 = cos(pkin(6));
t44 = -t53 * t79 + t82 * t83;
t74 = pkin(7) + qJ(3);
t63 = sin(t74) / 0.2e1;
t75 = pkin(7) - qJ(3);
t67 = sin(t75);
t51 = t63 + t67 / 0.2e1;
t65 = cos(t75) / 0.2e1;
t69 = cos(t74);
t56 = t65 + t69 / 0.2e1;
t59 = t66 - t70 / 0.2e1;
t87 = sin(qJ(3));
t32 = t51 * t83 + t53 * t56 - t59 * t87;
t78 = sin(pkin(8));
t81 = cos(pkin(8));
t23 = -t32 * t78 + t44 * t81;
t54 = t64 - t68 / 0.2e1;
t94 = cos(qJ(2));
t48 = -t89 * t54 + t94 * t95;
t28 = t124 * t51 + t47 * t56 - t48 * t87;
t19 = -t28 * t78 + t39 * t81;
t123 = t80 * t95;
t45 = t58 * t95 - t89 * t88;
t46 = t54 * t95 + t89 * t94;
t26 = -t123 * t51 + t45 * t56 - t46 * t87;
t127 = t45 * t79;
t38 = -t123 * t82 - t127;
t18 = -t26 * t78 + t38 * t81;
t121 = t95 * pkin(1) + pkin(10) * t124;
t120 = pkin(8) - qJ(4);
t119 = pkin(8) + qJ(4);
t118 = pkin(9) + r_base(3);
t115 = t83 * pkin(10) + t118;
t114 = cos(t119);
t113 = sin(t120);
t112 = -m(7) * pkin(14) + mrSges(6,2) - mrSges(7,3);
t111 = cos(t120) / 0.2e1;
t110 = sin(t119) / 0.2e1;
t72 = t89 * pkin(1);
t109 = t46 * pkin(2) - pkin(11) * t127 + t72;
t108 = -m(1) - m(2) - m(3) - t137;
t107 = t48 * pkin(2) + t39 * pkin(11) + t121;
t84 = sin(qJ(6));
t90 = cos(qJ(6));
t106 = -m(7) * pkin(5) - mrSges(7,1) * t90 + mrSges(7,2) * t84 - mrSges(6,1);
t105 = t59 * pkin(2) + t44 * pkin(11) + t115;
t104 = t111 + t114 / 0.2e1;
t103 = t110 + t113 / 0.2e1;
t52 = t63 - t67 / 0.2e1;
t57 = t65 - t69 / 0.2e1;
t93 = cos(qJ(3));
t27 = -t123 * t57 + t45 * t52 + t46 * t93;
t102 = t27 * pkin(3) + t18 * pkin(12) + t109;
t101 = -mrSges(7,1) * t84 - mrSges(7,2) * t90 - pkin(13) * t135 + mrSges(5,2) - mrSges(6,3);
t29 = t124 * t57 + t47 * t52 + t48 * t93;
t99 = t29 * pkin(3) + t19 * pkin(12) + t107;
t33 = t52 * t53 + t57 * t83 + t59 * t93;
t97 = t33 * pkin(3) + t23 * pkin(12) + t105;
t92 = cos(qJ(4));
t91 = cos(qJ(5));
t86 = sin(qJ(4));
t85 = sin(qJ(5));
t55 = t111 - t114 / 0.2e1;
t50 = t110 - t113 / 0.2e1;
t15 = t32 * t50 + t33 * t92 + t44 * t55;
t12 = t28 * t50 + t29 * t92 + t39 * t55;
t10 = t26 * t50 + t27 * t92 + t38 * t55;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t118 - mrSges(2,3) - m(3) * t115 - t59 * mrSges(3,1) - t53 * mrSges(3,2) - t83 * mrSges(3,3) - m(4) * t105 - t33 * mrSges(4,1) - t32 * mrSges(4,2) - t44 * mrSges(4,3) - m(5) * t97 - t15 * mrSges(5,1) - t23 * mrSges(5,3) + t112 * (t15 * t85 - t23 * t91) + t106 * (t15 * t91 + t23 * t85) + t101 * (-t103 * t44 - t104 * t32 + t33 * t86) + t135 * (-t15 * pkin(4) - t97)) * g(3) + (-m(5) * t102 - t89 * mrSges(2,1) - m(3) * t72 - m(4) * t109 - t45 * mrSges(3,2) - t46 * mrSges(3,1) - t38 * mrSges(4,3) - t18 * mrSges(5,3) - t26 * mrSges(4,2) - t27 * mrSges(4,1) - t10 * mrSges(5,1) - mrSges(1,2) + t108 * r_base(2) + t106 * (t10 * t91 + t18 * t85) + t101 * (-t103 * t38 - t104 * t26 + t27 * t86) + t112 * (t10 * t85 - t18 * t91) + (-mrSges(2,2) + (m(3) * pkin(10) + mrSges(3,3) + t137 * (pkin(11) * t82 + pkin(10))) * t80) * t95 + t135 * (-t10 * pkin(4) - t102)) * g(2) + (-t95 * mrSges(2,1) - m(4) * t107 - m(5) * t99 - m(3) * t121 - t47 * mrSges(3,2) - t48 * mrSges(3,1) - t39 * mrSges(4,3) - t19 * mrSges(5,3) - t28 * mrSges(4,2) - t29 * mrSges(4,1) - t12 * mrSges(5,1) - mrSges(1,1) + t108 * r_base(1) + (-t80 * mrSges(3,3) + mrSges(2,2)) * t89 + t112 * (t12 * t85 - t19 * t91) + t106 * (t12 * t91 + t19 * t85) + t101 * (-t103 * t39 - t28 * t104 + t29 * t86) + t135 * (-t12 * pkin(4) - t99)) * g(1);
U  = t1;
