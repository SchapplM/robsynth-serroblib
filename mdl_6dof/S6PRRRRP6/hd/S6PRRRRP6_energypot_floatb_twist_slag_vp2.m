% Calculate potential energy for
% S6PRRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,theta1]';
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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRRRP6_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP6_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PRRRRP6_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP6_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP6_energypot_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP6_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRP6_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:28:37
% EndTime: 2019-03-09 00:28:38
% DurationCPUTime: 0.76s
% Computational Cost: add. (621->116), mult. (1546->150), div. (0->0), fcn. (1949->14), ass. (0->60)
t50 = sin(pkin(7));
t53 = cos(pkin(7));
t54 = cos(pkin(6));
t51 = sin(pkin(6));
t60 = cos(qJ(2));
t85 = t51 * t60;
t71 = -t50 * t85 + t54 * t53;
t49 = sin(pkin(12));
t52 = cos(pkin(12));
t58 = sin(qJ(2));
t82 = t54 * t60;
t37 = -t49 * t82 - t52 * t58;
t87 = t51 * t53;
t72 = -t37 * t50 + t49 * t87;
t98 = -m(1) - m(2);
t97 = -m(6) - m(7);
t96 = mrSges(4,2) - mrSges(5,3);
t95 = mrSges(5,2) - mrSges(6,3) - mrSges(7,2);
t94 = -m(7) * pkin(5) - mrSges(6,1) - mrSges(7,1);
t93 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t92 = cos(qJ(3));
t91 = cos(qJ(4));
t89 = t49 * t51;
t88 = t51 * t52;
t86 = t51 * t58;
t83 = t54 * t58;
t79 = qJ(1) + r_base(3);
t78 = t50 * t92;
t77 = t53 * t92;
t76 = t52 * pkin(1) + pkin(8) * t89 + r_base(1);
t75 = t54 * pkin(8) + t79;
t74 = t51 * t78;
t35 = -t49 * t58 + t52 * t82;
t73 = -t35 * t50 - t52 * t87;
t70 = t49 * pkin(1) - pkin(8) * t88 + r_base(2);
t38 = -t49 * t83 + t52 * t60;
t69 = t38 * pkin(2) + t72 * pkin(9) + t76;
t68 = pkin(2) * t86 + t71 * pkin(9) + t75;
t57 = sin(qJ(3));
t20 = -t37 * t77 + t38 * t57 - t49 * t74;
t21 = t38 * t92 + (t37 * t53 + t50 * t89) * t57;
t67 = t21 * pkin(3) + pkin(10) * t20 + t69;
t28 = -t54 * t78 + t57 * t86 - t77 * t85;
t29 = t54 * t50 * t57 + (t53 * t57 * t60 + t92 * t58) * t51;
t66 = t29 * pkin(3) + t28 * pkin(10) + t68;
t36 = t49 * t60 + t52 * t83;
t65 = t36 * pkin(2) + t73 * pkin(9) + t70;
t18 = -t35 * t77 + t36 * t57 + t52 * t74;
t19 = t36 * t92 + (t35 * t53 - t50 * t88) * t57;
t62 = t19 * pkin(3) + pkin(10) * t18 + t65;
t59 = cos(qJ(5));
t56 = sin(qJ(4));
t55 = sin(qJ(5));
t23 = t29 * t91 + t71 * t56;
t22 = t29 * t56 - t71 * t91;
t12 = t21 * t91 + t56 * t72;
t11 = t21 * t56 - t72 * t91;
t10 = t19 * t91 + t56 * t73;
t9 = t19 * t56 - t73 * t91;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t79 - mrSges(2,3) - m(3) * t75 - t54 * mrSges(3,3) - (t58 * mrSges(3,1) + t60 * mrSges(3,2)) * t51 - m(4) * t68 - t29 * mrSges(4,1) - t71 * mrSges(4,3) - m(5) * t66 - t23 * mrSges(5,1) + t97 * (t23 * pkin(4) + t22 * pkin(11) + t66) + t94 * (t23 * t59 + t28 * t55) + t93 * (t23 * t55 - t28 * t59) + t96 * t28 + t95 * t22) * g(3) + (-m(3) * t70 - m(4) * t65 - m(5) * t62 - t49 * mrSges(2,1) - t36 * mrSges(3,1) - t19 * mrSges(4,1) - t10 * mrSges(5,1) - t52 * mrSges(2,2) - t35 * mrSges(3,2) + mrSges(3,3) * t88 - t73 * mrSges(4,3) - mrSges(1,2) + t98 * r_base(2) + t97 * (t10 * pkin(4) + pkin(11) * t9 + t62) + t94 * (t10 * t59 + t18 * t55) + t96 * t18 + t93 * (t10 * t55 - t18 * t59) + t95 * t9) * g(2) + (-m(3) * t76 - m(4) * t69 - m(5) * t67 - t52 * mrSges(2,1) - t38 * mrSges(3,1) - t21 * mrSges(4,1) - t12 * mrSges(5,1) + t49 * mrSges(2,2) - t37 * mrSges(3,2) - mrSges(3,3) * t89 - t72 * mrSges(4,3) - mrSges(1,1) + t98 * r_base(1) + t97 * (t12 * pkin(4) + pkin(11) * t11 + t67) + t94 * (t12 * t59 + t20 * t55) + t93 * (t12 * t55 - t20 * t59) + t96 * t20 + t95 * t11) * g(1);
U  = t1;
