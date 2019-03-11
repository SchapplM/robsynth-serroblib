% Calculate potential energy for
% S6RPRPRR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
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
% Datum: 2019-03-09 04:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRPRR9_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR9_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPRPRR9_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR9_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR9_energypot_floatb_twist_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR9_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR9_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:01:07
% EndTime: 2019-03-09 04:01:08
% DurationCPUTime: 1.06s
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
t50 = sin(pkin(6));
t52 = cos(pkin(12));
t88 = t50 * t52;
t58 = sin(qJ(1));
t87 = t50 * t58;
t62 = cos(qJ(1));
t86 = t50 * t62;
t53 = cos(pkin(7));
t54 = cos(pkin(6));
t85 = t53 * t54;
t57 = sin(qJ(3));
t84 = t53 * t57;
t83 = t54 * t62;
t48 = sin(pkin(12));
t82 = t58 * t48;
t81 = t58 * t52;
t80 = pkin(9) + qJ(4);
t79 = qJ(2) * t50;
t78 = pkin(8) + r_base(3);
t77 = t58 * pkin(1) + r_base(2);
t76 = t54 * qJ(2) + t78;
t75 = t62 * pkin(1) + t58 * t79 + r_base(1);
t47 = sin(pkin(13));
t51 = cos(pkin(13));
t61 = cos(qJ(3));
t74 = t47 * t61 + t51 * t57;
t40 = -t47 * t57 + t51 * t61;
t35 = t52 * t83 - t82;
t49 = sin(pkin(7));
t20 = -t35 * t49 - t53 * t86;
t73 = t35 * t53 - t49 * t86;
t37 = -t48 * t62 - t54 * t81;
t21 = -t37 * t49 + t53 * t87;
t72 = t37 * t53 + t49 * t87;
t32 = pkin(3) * t49 * t57 + t53 * t80;
t33 = pkin(3) * t84 - t49 * t80;
t43 = pkin(3) * t61 + pkin(2);
t71 = t50 * t48 * t43 + t54 * t32 + t33 * t88 + t76;
t70 = -t62 * t79 + t77;
t38 = t52 * t62 - t54 * t82;
t69 = t32 * t87 + t37 * t33 + t38 * t43 + t75;
t68 = t40 * t49;
t67 = t50 * t68;
t36 = t48 * t83 + t81;
t64 = t35 * t33 + t36 * t43 + (-qJ(2) - t32) * t86 + t77;
t60 = cos(qJ(5));
t56 = sin(qJ(5));
t34 = -t49 * t88 + t85;
t31 = t74 * t53;
t30 = t40 * t53;
t29 = t74 * t49;
t15 = t29 * t54 + (t31 * t52 + t40 * t48) * t50;
t14 = (t30 * t52 - t48 * t74) * t50 + t54 * t68;
t12 = t29 * t87 + t31 * t37 + t38 * t40;
t11 = t37 * t30 - t38 * t74 + t58 * t67;
t10 = -t29 * t86 + t35 * t31 + t36 * t40;
t9 = t30 * t35 - t36 * t74 - t62 * t67;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t78 - mrSges(2,3) - m(3) * t76 - m(4) * (pkin(9) * t85 + t76) - m(5) * t71 - t15 * mrSges(5,1) + t94 * (t15 * pkin(4) - pkin(10) * t14 + t71) + t91 * (t15 * t56 - t34 * t60) + (-mrSges(3,3) - (t57 * mrSges(4,1) + t61 * mrSges(4,2)) * t49) * t54 + (-t48 * mrSges(3,1) - t52 * mrSges(3,2) - m(4) * (-pkin(9) * t49 * t52 + pkin(2) * t48) - (t48 * t61 + t52 * t84) * mrSges(4,1) - (t52 * t53 * t61 - t48 * t57) * mrSges(4,2)) * t50 + t93 * t34 + t89 * (t15 * t60 + t34 * t56) + t90 * t14) * g(3) + (-(t36 * t61 + t57 * t73) * mrSges(4,1) - (-t36 * t57 + t61 * t73) * mrSges(4,2) - m(4) * (t36 * pkin(2) + t70) - m(5) * t64 - m(3) * t70 + mrSges(3,3) * t86 - t62 * mrSges(2,2) - t58 * mrSges(2,1) - t35 * mrSges(3,2) - t36 * mrSges(3,1) - t10 * mrSges(5,1) - mrSges(1,2) + t95 * r_base(2) + t94 * (t10 * pkin(4) - t9 * pkin(10) + t64) + t89 * (t10 * t60 + t20 * t56) + t90 * t9 + t92 * t20 + t91 * (t10 * t56 - t20 * t60)) * g(2) + (-m(4) * (pkin(2) * t38 + t75) - (t38 * t61 + t57 * t72) * mrSges(4,1) - (-t38 * t57 + t61 * t72) * mrSges(4,2) - m(5) * t69 - m(3) * t75 - mrSges(3,3) * t87 - t62 * mrSges(2,1) + t58 * mrSges(2,2) - t37 * mrSges(3,2) - t38 * mrSges(3,1) - t12 * mrSges(5,1) - mrSges(1,1) + t95 * r_base(1) + t94 * (t12 * pkin(4) - pkin(10) * t11 + t69) + t91 * (t12 * t56 - t21 * t60) + t92 * t21 + t89 * (t12 * t60 + t21 * t56) + t90 * t11) * g(1);
U  = t1;
