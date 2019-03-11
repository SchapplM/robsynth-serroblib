% Calculate potential energy for
% S6PPRPRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d5,d6,theta1,theta2,theta4]';
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
% Datum: 2019-03-08 18:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PPRPRR1_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRPRR1_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PPRPRR1_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRPRR1_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRPRR1_energypot_floatb_twist_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRPRR1_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRPRR1_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:41:54
% EndTime: 2019-03-08 18:41:55
% DurationCPUTime: 1.06s
% Computational Cost: add. (619->133), mult. (1478->178), div. (0->0), fcn. (1852->16), ass. (0->67)
t95 = -m(1) - m(2);
t94 = -m(6) - m(7);
t93 = -mrSges(4,3) - mrSges(5,3);
t92 = -m(4) * pkin(8) + t93;
t91 = -m(7) * pkin(10) + mrSges(6,2) - mrSges(7,3);
t57 = sin(qJ(6));
t60 = cos(qJ(6));
t90 = t57 * mrSges(7,1) + t60 * mrSges(7,2) - mrSges(5,2) + mrSges(6,3);
t89 = -m(7) * pkin(5) - t60 * mrSges(7,1) + t57 * mrSges(7,2) - mrSges(6,1);
t49 = sin(pkin(11));
t51 = sin(pkin(6));
t88 = t49 * t51;
t56 = cos(pkin(6));
t87 = t49 * t56;
t53 = cos(pkin(12));
t86 = t51 * t53;
t54 = cos(pkin(11));
t85 = t51 * t54;
t55 = cos(pkin(7));
t84 = t51 * t55;
t83 = t54 * t56;
t82 = t55 * t56;
t59 = sin(qJ(3));
t81 = t55 * t59;
t80 = pkin(8) + qJ(4);
t79 = qJ(2) * t51;
t78 = t49 * pkin(1) + r_base(2);
t77 = qJ(1) + r_base(3);
t76 = t54 * pkin(1) + t49 * t79 + r_base(1);
t75 = t56 * qJ(2) + t77;
t47 = sin(pkin(13));
t52 = cos(pkin(13));
t62 = cos(qJ(3));
t74 = t47 * t62 + t59 * t52;
t40 = -t59 * t47 + t52 * t62;
t48 = sin(pkin(12));
t35 = -t48 * t49 + t53 * t83;
t50 = sin(pkin(7));
t20 = -t35 * t50 - t54 * t84;
t73 = t35 * t55 - t50 * t85;
t37 = -t48 * t54 - t53 * t87;
t21 = -t37 * t50 + t49 * t84;
t72 = t37 * t55 + t50 * t88;
t71 = -t54 * t79 + t78;
t32 = pkin(3) * t50 * t59 + t55 * t80;
t33 = pkin(3) * t81 - t50 * t80;
t38 = -t48 * t87 + t53 * t54;
t43 = pkin(3) * t62 + pkin(2);
t70 = t32 * t88 + t37 * t33 + t38 * t43 + t76;
t69 = t40 * t50;
t68 = t51 * t48 * t43 + t56 * t32 + t33 * t86 + t75;
t67 = t51 * t69;
t36 = t48 * t83 + t49 * t53;
t65 = t35 * t33 + t36 * t43 + (-qJ(2) - t32) * t85 + t78;
t61 = cos(qJ(5));
t58 = sin(qJ(5));
t34 = -t50 * t86 + t82;
t31 = t74 * t55;
t30 = t40 * t55;
t29 = t74 * t50;
t15 = t29 * t56 + (t31 * t53 + t40 * t48) * t51;
t14 = (t30 * t53 - t48 * t74) * t51 + t56 * t69;
t10 = t29 * t88 + t31 * t37 + t38 * t40;
t9 = t37 * t30 - t38 * t74 + t49 * t67;
t8 = -t29 * t85 + t31 * t35 + t36 * t40;
t7 = t30 * t35 - t36 * t74 - t54 * t67;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t77 - mrSges(2,3) - m(3) * t75 - m(4) * (pkin(8) * t82 + t75) - m(5) * t68 - t15 * mrSges(5,1) + t94 * (t15 * pkin(4) - pkin(9) * t14 + t68) + (-mrSges(3,3) - (t59 * mrSges(4,1) + t62 * mrSges(4,2)) * t50) * t56 + (-t48 * mrSges(3,1) - t53 * mrSges(3,2) - m(4) * (-pkin(8) * t50 * t53 + pkin(2) * t48) - (t48 * t62 + t53 * t81) * mrSges(4,1) - (t53 * t55 * t62 - t48 * t59) * mrSges(4,2)) * t51 + t93 * t34 + t89 * (t15 * t61 + t34 * t58) + t90 * t14 + t91 * (t15 * t58 - t34 * t61)) * g(3) + (-m(4) * (pkin(2) * t36 + t71) - (t36 * t62 + t59 * t73) * mrSges(4,1) - (-t36 * t59 + t62 * t73) * mrSges(4,2) - m(5) * t65 - m(3) * t71 + mrSges(3,3) * t85 - t54 * mrSges(2,2) - t49 * mrSges(2,1) - t35 * mrSges(3,2) - t36 * mrSges(3,1) - t8 * mrSges(5,1) - mrSges(1,2) + t95 * r_base(2) + t94 * (t8 * pkin(4) - pkin(9) * t7 + t65) + t89 * (t20 * t58 + t61 * t8) + t90 * t7 + t92 * t20 + t91 * (-t20 * t61 + t58 * t8)) * g(2) + (-m(4) * (pkin(2) * t38 + t76) - (t38 * t62 + t59 * t72) * mrSges(4,1) - (-t38 * t59 + t62 * t72) * mrSges(4,2) - m(5) * t70 - m(3) * t76 - mrSges(3,3) * t88 - t54 * mrSges(2,1) + t49 * mrSges(2,2) - t37 * mrSges(3,2) - t38 * mrSges(3,1) - t10 * mrSges(5,1) - mrSges(1,1) + t95 * r_base(1) + t94 * (t10 * pkin(4) - pkin(9) * t9 + t70) + t89 * (t10 * t61 + t21 * t58) + t90 * t9 + t91 * (t10 * t58 - t21 * t61) + t92 * t21) * g(1);
U  = t1;
