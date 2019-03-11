% Calculate potential energy for
% S6PRRPRP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-03-08 21:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRPRP1_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP1_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PRRPRP1_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP1_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP1_energypot_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP1_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRP1_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:22:32
% EndTime: 2019-03-08 21:22:33
% DurationCPUTime: 0.81s
% Computational Cost: add. (364->119), mult. (559->136), div. (0->0), fcn. (633->12), ass. (0->51)
t75 = -m(1) - m(2);
t74 = -mrSges(6,1) - mrSges(7,1);
t73 = -mrSges(6,2) - mrSges(7,2);
t72 = m(4) * pkin(8) - mrSges(3,2) + mrSges(4,3);
t71 = -mrSges(5,3) - t72;
t70 = -m(6) * pkin(9) + m(7) * (-qJ(6) - pkin(9)) + mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t44 = sin(qJ(3));
t69 = pkin(3) * t44;
t37 = sin(pkin(10));
t39 = cos(pkin(10));
t45 = sin(qJ(2));
t40 = cos(pkin(6));
t48 = cos(qJ(2));
t59 = t40 * t48;
t17 = t37 * t45 - t39 * t59;
t43 = sin(qJ(5));
t68 = t17 * t43;
t19 = t37 * t59 + t39 * t45;
t67 = t19 * t43;
t38 = sin(pkin(6));
t66 = t37 * t38;
t65 = t38 * t39;
t64 = t38 * t44;
t63 = t38 * t45;
t47 = cos(qJ(3));
t62 = t38 * t47;
t61 = t38 * t48;
t60 = t40 * t45;
t58 = pkin(1) * t37 + r_base(2);
t57 = t37 * t64;
t56 = t43 * t61;
t55 = qJ(1) + r_base(3);
t54 = pkin(1) * t39 + pkin(7) * t66 + r_base(1);
t53 = pkin(7) * t40 + t55;
t52 = -pkin(7) * t65 + t58;
t20 = -t37 * t60 + t39 * t48;
t30 = pkin(3) * t47 + pkin(2);
t42 = -qJ(4) - pkin(8);
t51 = pkin(3) * t57 - t19 * t42 + t20 * t30 + t54;
t50 = t30 * t63 + t40 * t69 + t42 * t61 + t53;
t18 = t37 * t48 + t39 * t60;
t49 = t18 * t30 - t17 * t42 + (-pkin(7) - t69) * t65 + t58;
t46 = cos(qJ(5));
t36 = qJ(3) + pkin(11);
t32 = cos(t36);
t31 = sin(t36);
t29 = pkin(5) * t46 + pkin(4);
t14 = t31 * t40 + t32 * t63;
t8 = t20 * t32 + t31 * t66;
t6 = t18 * t32 - t31 * t65;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t55 - mrSges(2,3) - m(5) * t50 - t14 * mrSges(5,1) + mrSges(5,3) * t61 - m(6) * (pkin(4) * t14 + t50) - m(7) * (-pkin(5) * t56 + t14 * t29 + t50) + t73 * (-t14 * t43 - t46 * t61) + (-m(3) - m(4)) * t53 + (-t44 * mrSges(4,1) - t47 * mrSges(4,2) - mrSges(3,3)) * t40 + (t72 * t48 + (-m(4) * pkin(2) - mrSges(4,1) * t47 + mrSges(4,2) * t44 - mrSges(3,1)) * t45) * t38 + t74 * (t14 * t46 - t56) + t70 * (t31 * t63 - t32 * t40)) * g(3) + (-m(7) * (pkin(5) * t68 + t29 * t6 + t49) - (-t18 * t44 - t39 * t62) * mrSges(4,2) - (t18 * t47 - t39 * t64) * mrSges(4,1) + mrSges(3,3) * t65 - m(4) * (pkin(2) * t18 + t52) - m(3) * t52 - m(6) * (pkin(4) * t6 + t49) - m(5) * t49 - mrSges(1,2) - t6 * mrSges(5,1) - t18 * mrSges(3,1) - t37 * mrSges(2,1) - t39 * mrSges(2,2) + t75 * r_base(2) + t74 * (t46 * t6 + t68) + t73 * (t17 * t46 - t43 * t6) + t70 * (t18 * t31 + t32 * t65) + t71 * t17) * g(2) + (-m(7) * (pkin(5) * t67 + t29 * t8 + t51) - (-t20 * t44 + t37 * t62) * mrSges(4,2) - (t20 * t47 + t57) * mrSges(4,1) - m(4) * (pkin(2) * t20 + t54) - m(3) * t54 - m(6) * (pkin(4) * t8 + t51) - m(5) * t51 - mrSges(3,3) * t66 - mrSges(1,1) - t8 * mrSges(5,1) - t20 * mrSges(3,1) + t37 * mrSges(2,2) - t39 * mrSges(2,1) + t75 * r_base(1) + t74 * (t46 * t8 + t67) + t73 * (t19 * t46 - t43 * t8) + t70 * (t20 * t31 - t32 * t66) + t71 * t19) * g(1);
U  = t1;
