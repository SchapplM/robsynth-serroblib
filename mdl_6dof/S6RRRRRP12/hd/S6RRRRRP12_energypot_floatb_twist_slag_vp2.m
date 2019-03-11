% Calculate potential energy for
% S6RRRRRP12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5]';
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
% Datum: 2019-03-10 03:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRRP12_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP12_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRRRP12_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP12_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP12_energypot_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP12_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP12_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 02:59:31
% EndTime: 2019-03-10 02:59:32
% DurationCPUTime: 0.80s
% Computational Cost: add. (621->116), mult. (1546->147), div. (0->0), fcn. (1949->14), ass. (0->61)
t52 = cos(pkin(6));
t57 = sin(qJ(1));
t59 = cos(qJ(2));
t83 = t57 * t59;
t56 = sin(qJ(2));
t60 = cos(qJ(1));
t85 = t56 * t60;
t37 = -t52 * t83 - t85;
t49 = sin(pkin(7));
t51 = cos(pkin(7));
t50 = sin(pkin(6));
t89 = t50 * t57;
t72 = -t37 * t49 + t51 * t89;
t88 = t50 * t59;
t71 = -t49 * t88 + t51 * t52;
t99 = -m(1) - m(2);
t98 = -m(6) - m(7);
t97 = mrSges(4,2) - mrSges(5,3);
t96 = mrSges(5,2) - mrSges(6,3) - mrSges(7,2);
t95 = -m(7) * pkin(5) - mrSges(6,1) - mrSges(7,1);
t94 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t93 = cos(qJ(3));
t92 = cos(qJ(4));
t90 = t50 * t56;
t87 = t50 * t60;
t84 = t57 * t56;
t82 = t59 * t60;
t81 = pkin(8) + r_base(3);
t78 = t49 * t93;
t77 = t51 * t93;
t76 = t52 * pkin(9) + t81;
t75 = t60 * pkin(1) + pkin(9) * t89 + r_base(1);
t74 = t50 * t78;
t35 = t52 * t82 - t84;
t73 = -t35 * t49 - t51 * t87;
t70 = t57 * pkin(1) - pkin(9) * t87 + r_base(2);
t38 = -t52 * t84 + t82;
t69 = t38 * pkin(2) + t72 * pkin(10) + t75;
t68 = pkin(2) * t90 + t71 * pkin(10) + t76;
t55 = sin(qJ(3));
t22 = -t37 * t77 + t38 * t55 - t57 * t74;
t23 = t38 * t93 + (t37 * t51 + t49 * t89) * t55;
t67 = t23 * pkin(3) + pkin(11) * t22 + t69;
t28 = -t52 * t78 + t55 * t90 - t77 * t88;
t29 = t52 * t49 * t55 + (t51 * t55 * t59 + t56 * t93) * t50;
t66 = t29 * pkin(3) + pkin(11) * t28 + t68;
t36 = t52 * t85 + t83;
t65 = t36 * pkin(2) + pkin(10) * t73 + t70;
t20 = -t35 * t77 + t36 * t55 + t60 * t74;
t21 = t36 * t93 + (t35 * t51 - t49 * t87) * t55;
t62 = t21 * pkin(3) + t20 * pkin(11) + t65;
t58 = cos(qJ(5));
t54 = sin(qJ(4));
t53 = sin(qJ(5));
t19 = t29 * t92 + t54 * t71;
t18 = t29 * t54 - t71 * t92;
t12 = t23 * t92 + t54 * t72;
t11 = t23 * t54 - t72 * t92;
t10 = t21 * t92 + t54 * t73;
t9 = t21 * t54 - t73 * t92;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t81 - mrSges(2,3) - m(3) * t76 - t52 * mrSges(3,3) - (t56 * mrSges(3,1) + t59 * mrSges(3,2)) * t50 - m(4) * t68 - t29 * mrSges(4,1) - t71 * mrSges(4,3) - m(5) * t66 - t19 * mrSges(5,1) + t98 * (t19 * pkin(4) + pkin(12) * t18 + t66) + t95 * (t19 * t58 + t28 * t53) + t94 * (t19 * t53 - t28 * t58) + t97 * t28 + t96 * t18) * g(3) + (-m(3) * t70 - m(4) * t65 - m(5) * t62 - t57 * mrSges(2,1) - t36 * mrSges(3,1) - t21 * mrSges(4,1) - t10 * mrSges(5,1) - t60 * mrSges(2,2) - t35 * mrSges(3,2) + mrSges(3,3) * t87 - t73 * mrSges(4,3) - mrSges(1,2) + t99 * r_base(2) + t98 * (t10 * pkin(4) + t9 * pkin(12) + t62) + t97 * t20 + t95 * (t10 * t58 + t20 * t53) + t94 * (t10 * t53 - t20 * t58) + t96 * t9) * g(2) + (-m(3) * t75 - m(4) * t69 - m(5) * t67 - t60 * mrSges(2,1) - t38 * mrSges(3,1) - t23 * mrSges(4,1) - t12 * mrSges(5,1) + t57 * mrSges(2,2) - t37 * mrSges(3,2) - mrSges(3,3) * t89 - t72 * mrSges(4,3) - mrSges(1,1) + t99 * r_base(1) + t98 * (t12 * pkin(4) + pkin(12) * t11 + t67) + t95 * (t12 * t58 + t22 * t53) + t94 * (t12 * t53 - t22 * t58) + t97 * t22 + t96 * t11) * g(1);
U  = t1;
