% Calculate potential energy for
% S6PRPRPR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1]';
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
% Datum: 2019-03-08 19:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRPRPR7_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR7_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PRPRPR7_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR7_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRPR7_energypot_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR7_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR7_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:51:05
% EndTime: 2019-03-08 19:51:05
% DurationCPUTime: 0.68s
% Computational Cost: add. (270->100), mult. (528->104), div. (0->0), fcn. (592->10), ass. (0->55)
t73 = -m(1) - m(2);
t72 = mrSges(4,2) - mrSges(3,1);
t71 = mrSges(3,3) + mrSges(4,1);
t70 = -mrSges(4,3) + mrSges(3,2);
t69 = m(7) * pkin(9) + mrSges(5,1) - mrSges(6,2) + mrSges(7,3);
t33 = sin(qJ(6));
t36 = cos(qJ(6));
t68 = -t36 * mrSges(7,1) + t33 * mrSges(7,2) - mrSges(6,1) - mrSges(5,3);
t67 = t33 * mrSges(7,1) + t36 * mrSges(7,2) - mrSges(5,2) + mrSges(6,3);
t66 = -m(7) * (pkin(5) + pkin(8)) + t68 + t72;
t29 = sin(pkin(10));
t31 = cos(pkin(10));
t38 = cos(qJ(2));
t32 = cos(pkin(6));
t35 = sin(qJ(2));
t56 = t32 * t35;
t13 = t29 * t38 + t31 * t56;
t64 = pkin(8) * t13;
t15 = -t29 * t56 + t31 * t38;
t63 = pkin(8) * t15;
t30 = sin(pkin(6));
t62 = t29 * t30;
t61 = t30 * t31;
t34 = sin(qJ(4));
t60 = t30 * t34;
t59 = t30 * t35;
t37 = cos(qJ(4));
t58 = t30 * t37;
t57 = t30 * t38;
t55 = t32 * t38;
t54 = qJ(3) * t38;
t53 = pkin(7) * t61;
t52 = t29 * pkin(1) + r_base(2);
t51 = qJ(1) + r_base(3);
t50 = t30 * t54;
t49 = t31 * pkin(1) + pkin(7) * t62 + r_base(1);
t48 = t32 * pkin(7) + t51;
t47 = pkin(2) * t59 + t48;
t12 = t29 * t35 - t31 * t55;
t46 = t13 * pkin(2) + qJ(3) * t12 + t52;
t45 = t32 * pkin(3) + pkin(8) * t59 + t47;
t14 = t29 * t55 + t31 * t35;
t44 = t15 * pkin(2) + qJ(3) * t14 + t49;
t43 = pkin(3) * t62 + t44;
t16 = t32 * t34 + t37 * t57;
t17 = t32 * t37 - t34 * t57;
t42 = t17 * pkin(4) + t16 * qJ(5) + t45;
t41 = (-pkin(3) - pkin(7)) * t61 + t46;
t3 = -t14 * t37 + t29 * t60;
t4 = t14 * t34 + t29 * t58;
t40 = t4 * pkin(4) + qJ(5) * t3 + t43;
t5 = t12 * t37 + t31 * t60;
t6 = -t12 * t34 + t31 * t58;
t39 = -t6 * pkin(4) - qJ(5) * t5 + t41;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t51 - mrSges(2,3) - m(3) * t48 - m(4) * t47 - m(5) * (t45 - t50) - m(6) * (t42 - t50) - m(7) * t42 - t71 * t32 + (m(7) * t54 + (m(4) * qJ(3) - t70) * t38 + (-m(7) * pkin(5) + t72) * t35) * t30 + t68 * t59 - t67 * t16 - t69 * t17) * g(3) + (-m(5) * (t41 + t64) - m(6) * (t39 + t64) - m(7) * t39 - m(3) * (t52 - t53) - m(4) * (t46 - t53) - mrSges(1,2) - t29 * mrSges(2,1) - t31 * mrSges(2,2) + t73 * r_base(2) + t71 * t61 + t70 * t12 + t69 * t6 + t67 * t5 + t66 * t13) * g(2) + (-m(7) * t40 - m(5) * (t43 + t63) - m(6) * (t40 + t63) - mrSges(1,1) - m(3) * t49 - m(4) * t44 + t29 * mrSges(2,2) - t31 * mrSges(2,1) + t73 * r_base(1) - t71 * t62 + t70 * t14 - t69 * t4 - t67 * t3 + t66 * t15) * g(1);
U  = t1;
