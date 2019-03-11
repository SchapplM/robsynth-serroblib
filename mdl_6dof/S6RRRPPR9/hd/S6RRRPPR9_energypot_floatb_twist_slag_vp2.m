% Calculate potential energy for
% S6RRRPPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
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
% Datum: 2019-03-09 16:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRPPR9_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR9_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRPPR9_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR9_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR9_energypot_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR9_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR9_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:11:31
% EndTime: 2019-03-09 16:11:31
% DurationCPUTime: 0.65s
% Computational Cost: add. (382->105), mult. (821->117), div. (0->0), fcn. (997->12), ass. (0->57)
t75 = -m(2) - m(1);
t74 = -mrSges(4,3) + mrSges(3,2);
t73 = mrSges(4,2) - mrSges(5,3) - mrSges(6,2) - m(7) * (-pkin(10) + qJ(4)) + mrSges(7,3);
t38 = sin(qJ(6));
t42 = cos(qJ(6));
t72 = -t38 * mrSges(7,1) - t42 * mrSges(7,2) + mrSges(5,2) - mrSges(6,3);
t71 = -m(7) * pkin(5) - t42 * mrSges(7,1) + t38 * mrSges(7,2) - mrSges(5,1) - mrSges(6,1);
t70 = cos(qJ(3));
t36 = sin(pkin(6));
t40 = sin(qJ(2));
t69 = t36 * t40;
t41 = sin(qJ(1));
t68 = t36 * t41;
t43 = cos(qJ(2));
t67 = t36 * t43;
t44 = cos(qJ(1));
t66 = t36 * t44;
t62 = cos(pkin(6));
t58 = t41 * t62;
t26 = -t40 * t58 + t44 * t43;
t39 = sin(qJ(3));
t60 = t36 * t70;
t14 = t26 * t39 - t41 * t60;
t65 = qJ(4) * t14;
t21 = t39 * t69 - t62 * t70;
t64 = qJ(4) * t21;
t57 = t44 * t62;
t24 = t40 * t57 + t41 * t43;
t12 = t24 * t39 + t44 * t60;
t63 = t12 * qJ(4);
t61 = pkin(7) + r_base(3);
t59 = t62 * pkin(8) + t61;
t56 = t44 * pkin(1) + pkin(8) * t68 + r_base(1);
t54 = t41 * pkin(1) - pkin(8) * t66 + r_base(2);
t25 = t44 * t40 + t43 * t58;
t53 = t26 * pkin(2) + pkin(9) * t25 + t56;
t52 = pkin(2) * t69 - pkin(9) * t67 + t59;
t15 = t26 * t70 + t39 * t68;
t51 = t15 * pkin(3) + t53;
t22 = t39 * t62 + t40 * t60;
t50 = t22 * pkin(3) + t52;
t23 = t40 * t41 - t43 * t57;
t49 = t24 * pkin(2) + t23 * pkin(9) + t54;
t13 = t24 * t70 - t39 * t66;
t48 = t13 * pkin(3) + t49;
t35 = sin(pkin(11));
t37 = cos(pkin(11));
t5 = t15 * t35 - t25 * t37;
t6 = t15 * t37 + t25 * t35;
t47 = t6 * pkin(4) + qJ(5) * t5 + t51;
t10 = t22 * t35 + t37 * t67;
t11 = t22 * t37 - t35 * t67;
t46 = t11 * pkin(4) + qJ(5) * t10 + t50;
t3 = t13 * t35 - t23 * t37;
t4 = t13 * t37 + t23 * t35;
t45 = t4 * pkin(4) + t3 * qJ(5) + t48;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t61 - mrSges(2,3) - m(3) * t59 - t62 * mrSges(3,3) - (mrSges(3,1) * t40 + mrSges(3,2) * t43) * t36 - m(4) * t52 - t22 * mrSges(4,1) + mrSges(4,3) * t67 - m(5) * (t50 + t64) - m(6) * (t46 + t64) - m(7) * t46 + t71 * t11 + t72 * t10 + t73 * t21) * g(3) + (-m(5) * (t48 + t63) - m(6) * (t45 + t63) - m(7) * t45 - m(4) * t49 - m(3) * t54 + mrSges(3,3) * t66 - mrSges(1,2) - t13 * mrSges(4,1) - t24 * mrSges(3,1) - t41 * mrSges(2,1) - t44 * mrSges(2,2) + t75 * r_base(2) + t71 * t4 + t72 * t3 + t74 * t23 + t73 * t12) * g(2) + (-m(5) * (t51 + t65) - m(6) * (t47 + t65) - m(7) * t47 - m(4) * t53 - m(3) * t56 - mrSges(3,3) * t68 - mrSges(1,1) - t15 * mrSges(4,1) - t26 * mrSges(3,1) + t41 * mrSges(2,2) - t44 * mrSges(2,1) + t75 * r_base(1) + t71 * t6 + t72 * t5 + t74 * t25 + t73 * t14) * g(1);
U  = t1;
