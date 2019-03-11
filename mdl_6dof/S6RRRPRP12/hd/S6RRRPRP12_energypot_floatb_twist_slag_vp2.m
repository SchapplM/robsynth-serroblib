% Calculate potential energy for
% S6RRRPRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5]';
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
% Datum: 2019-03-09 18:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRPRP12_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP12_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRPRP12_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP12_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP12_energypot_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP12_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP12_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:51:41
% EndTime: 2019-03-09 17:51:42
% DurationCPUTime: 0.62s
% Computational Cost: add. (322->97), mult. (664->110), div. (0->0), fcn. (780->10), ass. (0->47)
t73 = -m(1) - m(2);
t72 = -m(6) - m(7);
t71 = mrSges(4,2) - mrSges(5,3);
t70 = mrSges(4,3) + mrSges(5,1);
t69 = mrSges(3,2) - t70;
t68 = -mrSges(4,1) + mrSges(5,2) - mrSges(6,3) - mrSges(7,2);
t67 = -m(7) * pkin(5) - mrSges(6,1) - mrSges(7,1);
t66 = m(7) * qJ(6) - mrSges(6,2) + mrSges(7,3);
t65 = cos(qJ(3));
t35 = sin(pkin(6));
t38 = sin(qJ(2));
t64 = t35 * t38;
t39 = sin(qJ(1));
t63 = t35 * t39;
t41 = cos(qJ(2));
t62 = t35 * t41;
t42 = cos(qJ(1));
t61 = t35 * t42;
t60 = cos(pkin(6));
t59 = pkin(7) + r_base(3);
t58 = pkin(9) * t62;
t57 = t35 * t65;
t56 = t60 * pkin(8) + t59;
t55 = t39 * t60;
t54 = t42 * t60;
t53 = t42 * pkin(1) + pkin(8) * t63 + r_base(1);
t52 = pkin(2) * t64 + t56;
t51 = t39 * pkin(1) - pkin(8) * t61 + r_base(2);
t24 = t42 * t38 + t41 * t55;
t25 = -t38 * t55 + t42 * t41;
t50 = t25 * pkin(2) + pkin(9) * t24 + t53;
t37 = sin(qJ(3));
t20 = t37 * t64 - t60 * t65;
t21 = t60 * t37 + t38 * t57;
t49 = t21 * pkin(3) + qJ(4) * t20 + t52;
t22 = t38 * t39 - t41 * t54;
t23 = t38 * t54 + t39 * t41;
t48 = t23 * pkin(2) + t22 * pkin(9) + t51;
t13 = t25 * t37 - t39 * t57;
t14 = t25 * t65 + t37 * t63;
t47 = t14 * pkin(3) + qJ(4) * t13 + t50;
t11 = t23 * t37 + t42 * t57;
t12 = t23 * t65 - t37 * t61;
t46 = t12 * pkin(3) + t11 * qJ(4) + t48;
t40 = cos(qJ(5));
t36 = sin(qJ(5));
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t59 - mrSges(2,3) - m(3) * t56 - t60 * mrSges(3,3) - (t38 * mrSges(3,1) + t41 * mrSges(3,2)) * t35 - m(4) * (t52 - t58) - m(5) * (t49 - t58) + t72 * (pkin(10) * t21 + (-pkin(4) - pkin(9)) * t62 + t49) + t66 * (t20 * t40 + t36 * t62) + t70 * t62 + t71 * t20 + t67 * (t20 * t36 - t40 * t62) + t68 * t21) * g(3) + (-m(3) * t51 - m(4) * t48 - m(5) * t46 - t39 * mrSges(2,1) - t23 * mrSges(3,1) - t42 * mrSges(2,2) + mrSges(3,3) * t61 - mrSges(1,2) + t73 * r_base(2) + t72 * (t22 * pkin(4) + t12 * pkin(10) + t46) + t67 * (t11 * t36 + t22 * t40) + t71 * t11 - t66 * (-t11 * t40 + t22 * t36) + t69 * t22 + t68 * t12) * g(2) + (-m(3) * t53 - m(4) * t50 - m(5) * t47 - t42 * mrSges(2,1) - t25 * mrSges(3,1) + t39 * mrSges(2,2) - mrSges(3,3) * t63 - mrSges(1,1) + t73 * r_base(1) + t72 * (t24 * pkin(4) + pkin(10) * t14 + t47) + t67 * (t13 * t36 + t24 * t40) - t66 * (-t13 * t40 + t24 * t36) + t71 * t13 + t69 * t24 + t68 * t14) * g(1);
U  = t1;
