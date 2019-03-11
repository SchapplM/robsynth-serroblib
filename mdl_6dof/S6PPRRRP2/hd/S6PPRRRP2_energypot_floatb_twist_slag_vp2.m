% Calculate potential energy for
% S6PPRRRP2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-03-08 18:58
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PPRRRP2_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP2_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PPRRRP2_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRP2_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP2_energypot_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRP2_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRRP2_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:55:55
% EndTime: 2019-03-08 18:55:56
% DurationCPUTime: 0.80s
% Computational Cost: add. (621->116), mult. (1546->151), div. (0->0), fcn. (1949->14), ass. (0->61)
t51 = sin(pkin(7));
t55 = cos(pkin(7));
t56 = cos(pkin(6));
t52 = sin(pkin(6));
t53 = cos(pkin(12));
t87 = t52 * t53;
t71 = -t51 * t87 + t55 * t56;
t49 = sin(pkin(12));
t54 = cos(pkin(11));
t50 = sin(pkin(11));
t88 = t50 * t56;
t37 = -t49 * t54 - t53 * t88;
t85 = t52 * t55;
t72 = -t37 * t51 + t50 * t85;
t99 = -m(1) - m(2);
t98 = -m(6) - m(7);
t97 = mrSges(4,2) - mrSges(5,3);
t96 = mrSges(5,2) - mrSges(6,3) - mrSges(7,2);
t95 = -m(7) * pkin(5) - mrSges(6,1) - mrSges(7,1);
t94 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t93 = cos(qJ(3));
t92 = cos(qJ(4));
t90 = t49 * t52;
t89 = t50 * t52;
t86 = t52 * t54;
t84 = t54 * t56;
t82 = qJ(2) * t52;
t79 = qJ(1) + r_base(3);
t78 = t51 * t93;
t77 = t55 * t93;
t76 = t54 * pkin(1) + t50 * t82 + r_base(1);
t75 = t56 * qJ(2) + t79;
t74 = t52 * t78;
t35 = -t49 * t50 + t53 * t84;
t73 = -t35 * t51 - t54 * t85;
t70 = t50 * pkin(1) - t54 * t82 + r_base(2);
t38 = -t49 * t88 + t53 * t54;
t69 = t38 * pkin(2) + t72 * pkin(8) + t76;
t68 = pkin(2) * t90 + t71 * pkin(8) + t75;
t59 = sin(qJ(3));
t20 = -t37 * t77 + t38 * t59 - t50 * t74;
t21 = t38 * t93 + (t37 * t55 + t51 * t89) * t59;
t67 = t21 * pkin(3) + pkin(9) * t20 + t69;
t28 = -t56 * t78 + t59 * t90 - t77 * t87;
t29 = t56 * t51 * t59 + (t53 * t55 * t59 + t49 * t93) * t52;
t66 = t29 * pkin(3) + pkin(9) * t28 + t68;
t36 = t49 * t84 + t50 * t53;
t64 = t36 * pkin(2) + pkin(8) * t73 + t70;
t18 = -t35 * t77 + t36 * t59 + t54 * t74;
t19 = t36 * t93 + (t35 * t55 - t51 * t86) * t59;
t62 = t19 * pkin(3) + pkin(9) * t18 + t64;
t60 = cos(qJ(5));
t58 = sin(qJ(4));
t57 = sin(qJ(5));
t23 = t29 * t92 + t58 * t71;
t22 = t29 * t58 - t71 * t92;
t10 = t21 * t92 + t58 * t72;
t9 = t21 * t58 - t72 * t92;
t8 = t19 * t92 + t58 * t73;
t7 = t19 * t58 - t73 * t92;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t79 - mrSges(2,3) - m(3) * t75 - t56 * mrSges(3,3) - (t49 * mrSges(3,1) + t53 * mrSges(3,2)) * t52 - m(4) * t68 - t29 * mrSges(4,1) - t71 * mrSges(4,3) - m(5) * t66 - t23 * mrSges(5,1) + t98 * (t23 * pkin(4) + pkin(10) * t22 + t66) + t97 * t28 + t95 * (t23 * t60 + t28 * t57) + t94 * (t23 * t57 - t28 * t60) + t96 * t22) * g(3) + (-m(3) * t70 - m(4) * t64 - m(5) * t62 - t50 * mrSges(2,1) - t36 * mrSges(3,1) - t19 * mrSges(4,1) - t8 * mrSges(5,1) - t54 * mrSges(2,2) - t35 * mrSges(3,2) + mrSges(3,3) * t86 - t73 * mrSges(4,3) - mrSges(1,2) + t99 * r_base(2) + t98 * (t8 * pkin(4) + pkin(10) * t7 + t62) + t95 * (t18 * t57 + t60 * t8) + t97 * t18 + t94 * (-t18 * t60 + t57 * t8) + t96 * t7) * g(2) + (-m(3) * t76 - m(4) * t69 - m(5) * t67 - t54 * mrSges(2,1) - t38 * mrSges(3,1) - t21 * mrSges(4,1) - t10 * mrSges(5,1) + t50 * mrSges(2,2) - t37 * mrSges(3,2) - mrSges(3,3) * t89 - t72 * mrSges(4,3) - mrSges(1,1) + t99 * r_base(1) + t98 * (t10 * pkin(4) + pkin(10) * t9 + t67) + t95 * (t10 * t60 + t20 * t57) + t94 * (t10 * t57 - t20 * t60) + t97 * t20 + t96 * t9) * g(1);
U  = t1;
