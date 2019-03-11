% Calculate potential energy for
% S6PRRPRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
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
% Datum: 2019-03-08 22:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRPRR3_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR3_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PRRPRR3_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR3_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR3_energypot_floatb_twist_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR3_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR3_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:02:06
% EndTime: 2019-03-08 22:02:07
% DurationCPUTime: 1.04s
% Computational Cost: add. (619->133), mult. (1478->176), div. (0->0), fcn. (1852->16), ass. (0->67)
t95 = -m(1) - m(2);
t94 = -m(6) - m(7);
t93 = -mrSges(4,3) - mrSges(5,3);
t92 = -m(4) * pkin(9) + t93;
t91 = -m(7) * pkin(11) + mrSges(6,2) - mrSges(7,3);
t55 = sin(qJ(6));
t59 = cos(qJ(6));
t90 = t55 * mrSges(7,1) + t59 * mrSges(7,2) - mrSges(5,2) + mrSges(6,3);
t89 = -m(7) * pkin(5) - t59 * mrSges(7,1) + t55 * mrSges(7,2) - mrSges(6,1);
t57 = sin(qJ(3));
t88 = pkin(3) * t57;
t48 = sin(pkin(12));
t50 = sin(pkin(6));
t87 = t48 * t50;
t52 = cos(pkin(12));
t86 = t50 * t52;
t53 = cos(pkin(7));
t85 = t50 * t53;
t62 = cos(qJ(2));
t84 = t50 * t62;
t83 = t53 * t62;
t54 = cos(pkin(6));
t82 = t54 * t53;
t58 = sin(qJ(2));
t81 = t54 * t58;
t80 = t54 * t62;
t79 = pkin(9) + qJ(4);
t78 = t48 * pkin(1) + r_base(2);
t77 = qJ(1) + r_base(3);
t76 = t52 * pkin(1) + pkin(8) * t87 + r_base(1);
t75 = t54 * pkin(8) + t77;
t47 = sin(pkin(13));
t51 = cos(pkin(13));
t61 = cos(qJ(3));
t74 = t47 * t61 + t51 * t57;
t40 = -t47 * t57 + t51 * t61;
t35 = -t48 * t58 + t52 * t80;
t49 = sin(pkin(7));
t20 = -t35 * t49 - t52 * t85;
t73 = t35 * t53 - t49 * t86;
t37 = -t48 * t80 - t52 * t58;
t21 = -t37 * t49 + t48 * t85;
t72 = t37 * t53 + t49 * t87;
t71 = -pkin(8) * t86 + t78;
t32 = t49 * t88 + t53 * t79;
t33 = -t49 * t79 + t53 * t88;
t38 = -t48 * t81 + t52 * t62;
t43 = pkin(3) * t61 + pkin(2);
t70 = t32 * t87 + t37 * t33 + t38 * t43 + t76;
t69 = t40 * t49;
t68 = t50 * t58 * t43 + t54 * t32 + t33 * t84 + t75;
t67 = t50 * t69;
t36 = t48 * t62 + t52 * t81;
t66 = t35 * t33 + t36 * t43 + (-pkin(8) - t32) * t86 + t78;
t60 = cos(qJ(5));
t56 = sin(qJ(5));
t34 = -t49 * t84 + t82;
t31 = t74 * t53;
t30 = t40 * t53;
t29 = t74 * t49;
t15 = t54 * t29 + (t31 * t62 + t40 * t58) * t50;
t14 = (t30 * t62 - t58 * t74) * t50 + t54 * t69;
t10 = t29 * t87 + t31 * t37 + t38 * t40;
t9 = t37 * t30 - t38 * t74 + t48 * t67;
t8 = -t29 * t86 + t31 * t35 + t36 * t40;
t7 = t30 * t35 - t36 * t74 - t52 * t67;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t77 - mrSges(2,3) - m(3) * t75 - m(4) * (pkin(9) * t82 + t75) - m(5) * t68 - t15 * mrSges(5,1) + t94 * (t15 * pkin(4) - pkin(10) * t14 + t68) + (-mrSges(3,3) - (t57 * mrSges(4,1) + t61 * mrSges(4,2)) * t49) * t54 + (-t58 * mrSges(3,1) - t62 * mrSges(3,2) - m(4) * (-pkin(9) * t49 * t62 + pkin(2) * t58) - (t57 * t83 + t58 * t61) * mrSges(4,1) - (-t57 * t58 + t61 * t83) * mrSges(4,2)) * t50 + t93 * t34 + t89 * (t15 * t60 + t34 * t56) + t90 * t14 + t91 * (t15 * t56 - t34 * t60)) * g(3) + (-(t36 * t61 + t57 * t73) * mrSges(4,1) - (-t36 * t57 + t61 * t73) * mrSges(4,2) - m(4) * (pkin(2) * t36 + t71) - m(5) * t66 - m(3) * t71 + mrSges(3,3) * t86 - t48 * mrSges(2,1) - t52 * mrSges(2,2) - t35 * mrSges(3,2) - t36 * mrSges(3,1) - t8 * mrSges(5,1) - mrSges(1,2) + t95 * r_base(2) + t94 * (t8 * pkin(4) - pkin(10) * t7 + t66) + t89 * (t20 * t56 + t60 * t8) + t90 * t7 + t92 * t20 + t91 * (-t20 * t60 + t56 * t8)) * g(2) + (-m(4) * (pkin(2) * t38 + t76) - (t38 * t61 + t57 * t72) * mrSges(4,1) - (-t38 * t57 + t61 * t72) * mrSges(4,2) - m(5) * t70 - m(3) * t76 - mrSges(3,3) * t87 + t48 * mrSges(2,2) - t52 * mrSges(2,1) - t37 * mrSges(3,2) - t38 * mrSges(3,1) - t10 * mrSges(5,1) - mrSges(1,1) + t95 * r_base(1) + t94 * (t10 * pkin(4) - pkin(10) * t9 + t70) + t89 * (t10 * t60 + t21 * t56) + t90 * t9 + t91 * (t10 * t56 - t21 * t60) + t92 * t21) * g(1);
U  = t1;
