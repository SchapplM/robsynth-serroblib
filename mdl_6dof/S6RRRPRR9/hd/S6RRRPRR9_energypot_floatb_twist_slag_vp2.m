% Calculate potential energy for
% S6RRRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 19:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRPRR9_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR9_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRPRR9_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR9_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR9_energypot_floatb_twist_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR9_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR9_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:56:25
% EndTime: 2019-03-09 18:56:26
% DurationCPUTime: 1.04s
% Computational Cost: add. (619->133), mult. (1478->173), div. (0->0), fcn. (1852->16), ass. (0->68)
t96 = -m(1) - m(2);
t95 = -m(6) - m(7);
t94 = -mrSges(4,3) - mrSges(5,3);
t93 = -m(4) * pkin(10) + t94;
t92 = -m(7) * pkin(12) + mrSges(6,2) - mrSges(7,3);
t53 = sin(qJ(6));
t58 = cos(qJ(6));
t91 = mrSges(7,1) * t53 + mrSges(7,2) * t58 - mrSges(5,2) + mrSges(6,3);
t90 = -m(7) * pkin(5) - mrSges(7,1) * t58 + mrSges(7,2) * t53 - mrSges(6,1);
t55 = sin(qJ(3));
t89 = pkin(3) * t55;
t49 = sin(pkin(6));
t57 = sin(qJ(1));
t88 = t49 * t57;
t61 = cos(qJ(2));
t87 = t49 * t61;
t62 = cos(qJ(1));
t86 = t49 * t62;
t51 = cos(pkin(7));
t52 = cos(pkin(6));
t85 = t51 * t52;
t84 = t51 * t61;
t56 = sin(qJ(2));
t83 = t56 * t62;
t82 = t57 * t56;
t81 = t57 * t61;
t80 = t61 * t62;
t79 = pkin(10) + qJ(4);
t78 = pkin(8) + r_base(3);
t77 = pkin(1) * t57 + r_base(2);
t76 = pkin(9) * t52 + t78;
t75 = pkin(1) * t62 + pkin(9) * t88 + r_base(1);
t47 = sin(pkin(13));
t50 = cos(pkin(13));
t60 = cos(qJ(3));
t74 = t47 * t60 + t50 * t55;
t40 = -t47 * t55 + t50 * t60;
t35 = t52 * t80 - t82;
t48 = sin(pkin(7));
t20 = -t35 * t48 - t51 * t86;
t73 = t35 * t51 - t48 * t86;
t37 = -t52 * t81 - t83;
t21 = -t37 * t48 + t51 * t88;
t72 = t37 * t51 + t48 * t88;
t71 = -pkin(9) * t86 + t77;
t32 = t48 * t89 + t51 * t79;
t33 = -t48 * t79 + t51 * t89;
t43 = pkin(3) * t60 + pkin(2);
t70 = t43 * t49 * t56 + t32 * t52 + t33 * t87 + t76;
t38 = -t52 * t82 + t80;
t69 = t32 * t88 + t33 * t37 + t38 * t43 + t75;
t68 = t40 * t48;
t67 = t49 * t68;
t36 = t52 * t83 + t81;
t66 = t35 * t33 + t36 * t43 + (-pkin(9) - t32) * t86 + t77;
t59 = cos(qJ(5));
t54 = sin(qJ(5));
t34 = -t48 * t87 + t85;
t31 = t74 * t51;
t30 = t40 * t51;
t29 = t74 * t48;
t15 = t29 * t52 + (t31 * t61 + t40 * t56) * t49;
t14 = (t30 * t61 - t56 * t74) * t49 + t52 * t68;
t12 = t29 * t88 + t31 * t37 + t38 * t40;
t11 = t37 * t30 - t38 * t74 + t57 * t67;
t10 = -t29 * t86 + t31 * t35 + t36 * t40;
t9 = t30 * t35 - t36 * t74 - t62 * t67;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t78 - mrSges(2,3) - m(3) * t76 - m(4) * (pkin(10) * t85 + t76) - m(5) * t70 - t15 * mrSges(5,1) + t95 * (pkin(4) * t15 - pkin(11) * t14 + t70) + t92 * (t15 * t54 - t34 * t59) + (-mrSges(3,3) - (mrSges(4,1) * t55 + mrSges(4,2) * t60) * t48) * t52 + (-t56 * mrSges(3,1) - t61 * mrSges(3,2) - m(4) * (-pkin(10) * t48 * t61 + pkin(2) * t56) - (t55 * t84 + t56 * t60) * mrSges(4,1) - (-t55 * t56 + t60 * t84) * mrSges(4,2)) * t49 + t94 * t34 + t90 * (t15 * t59 + t34 * t54) + t91 * t14) * g(3) + (mrSges(3,3) * t86 - (t36 * t60 + t55 * t73) * mrSges(4,1) - (-t36 * t55 + t60 * t73) * mrSges(4,2) - m(4) * (pkin(2) * t36 + t71) - m(3) * t71 - m(5) * t66 - t62 * mrSges(2,2) - t57 * mrSges(2,1) - t35 * mrSges(3,2) - t36 * mrSges(3,1) - t10 * mrSges(5,1) - mrSges(1,2) + t96 * r_base(2) + t95 * (pkin(4) * t10 - t9 * pkin(11) + t66) + t90 * (t10 * t59 + t20 * t54) + t91 * t9 + t93 * t20 + t92 * (t10 * t54 - t20 * t59)) * g(2) + (-m(4) * (pkin(2) * t38 + t75) - m(3) * t75 - (t38 * t60 + t55 * t72) * mrSges(4,1) - (-t38 * t55 + t60 * t72) * mrSges(4,2) - m(5) * t69 - t62 * mrSges(2,1) + t57 * mrSges(2,2) - t37 * mrSges(3,2) - t38 * mrSges(3,1) - t12 * mrSges(5,1) - mrSges(1,1) - mrSges(3,3) * t88 + t96 * r_base(1) + t95 * (pkin(4) * t12 - pkin(11) * t11 + t69) + t92 * (t12 * t54 - t21 * t59) + t93 * t21 + t90 * (t12 * t59 + t21 * t54) + t91 * t11) * g(1);
U  = t1;
