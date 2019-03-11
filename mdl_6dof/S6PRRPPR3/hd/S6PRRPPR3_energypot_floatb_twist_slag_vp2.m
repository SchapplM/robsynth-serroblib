% Calculate potential energy for
% S6PRRPPR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1]';
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
% Datum: 2019-03-08 21:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRPPR3_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR3_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PRRPPR3_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR3_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPPR3_energypot_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR3_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPPR3_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:08:03
% EndTime: 2019-03-08 21:08:03
% DurationCPUTime: 0.65s
% Computational Cost: add. (299->98), mult. (604->99), div. (0->0), fcn. (697->10), ass. (0->52)
t32 = sin(qJ(6));
t35 = cos(qJ(6));
t71 = -t32 * mrSges(7,1) - t35 * mrSges(7,2) + mrSges(5,2) + mrSges(4,3) - mrSges(6,3);
t70 = -m(1) - m(2);
t33 = sin(qJ(3));
t34 = sin(qJ(2));
t30 = sin(pkin(6));
t62 = cos(qJ(3));
t52 = t30 * t62;
t54 = cos(pkin(6));
t18 = t33 * t54 + t34 * t52;
t36 = cos(qJ(2));
t57 = t30 * t36;
t69 = pkin(4) * t18 + qJ(5) * t57;
t68 = -m(7) * pkin(9) - mrSges(4,1) - mrSges(5,1) + mrSges(6,2) - mrSges(7,3);
t66 = -m(7) * (pkin(5) + qJ(4)) - t35 * mrSges(7,1) + t32 * mrSges(7,2) - mrSges(6,1) + mrSges(4,2) - mrSges(5,3);
t65 = mrSges(3,2) - t71 + (-m(6) - m(7)) * (pkin(8) - qJ(5));
t29 = sin(pkin(10));
t31 = cos(pkin(10));
t50 = t36 * t54;
t13 = t29 * t34 - t31 * t50;
t64 = pkin(8) * t13;
t15 = t29 * t50 + t31 * t34;
t63 = pkin(8) * t15;
t61 = t29 * t30;
t60 = t30 * t31;
t59 = t30 * t33;
t58 = t30 * t34;
t53 = qJ(1) + r_base(3);
t51 = t34 * t54;
t49 = pkin(1) * t31 + pkin(7) * t61 + r_base(1);
t48 = pkin(7) * t54 + t53;
t16 = -t29 * t51 + t31 * t36;
t47 = pkin(2) * t16 + t49;
t8 = t16 * t62 + t29 * t59;
t45 = pkin(3) * t8 + t47;
t44 = pkin(1) * t29 - pkin(7) * t60 + r_base(2);
t14 = t29 * t36 + t31 * t51;
t43 = pkin(2) * t14 + t44;
t6 = t14 * t62 - t31 * t59;
t42 = pkin(3) * t6 + t43;
t7 = t16 * t33 - t29 * t52;
t41 = qJ(4) * t7 + t45;
t40 = pkin(2) * t58 - pkin(8) * t57 + t48;
t39 = pkin(3) * t18 + t40;
t5 = t14 * t33 + t31 * t52;
t38 = qJ(4) * t5 + t42;
t17 = t33 * t58 - t54 * t62;
t37 = qJ(4) * t17 + t39;
t3 = t8 * pkin(4);
t1 = t6 * pkin(4);
t2 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t53 - mrSges(2,3) - m(3) * t48 - t54 * mrSges(3,3) - (t34 * mrSges(3,1) + t36 * mrSges(3,2)) * t30 - m(4) * t40 - m(5) * t37 - m(6) * (t37 + t69) - m(7) * (t39 + t69) + t71 * t57 + t66 * t17 + t68 * t18) * g(3) + (-m(4) * (t43 + t64) - m(5) * (t38 + t64) + mrSges(3,3) * t60 - m(7) * (t1 + t42) - m(3) * t44 - m(6) * (t1 + t38) - mrSges(1,2) - t14 * mrSges(3,1) - t29 * mrSges(2,1) - t31 * mrSges(2,2) + t70 * r_base(2) + t68 * t6 + t66 * t5 + t65 * t13) * g(2) + (-mrSges(3,3) * t61 - m(4) * (t47 + t63) - m(5) * (t41 + t63) - m(7) * (t3 + t45) - m(3) * t49 - m(6) * (t3 + t41) - mrSges(1,1) - t16 * mrSges(3,1) + t29 * mrSges(2,2) - t31 * mrSges(2,1) + t70 * r_base(1) + t68 * t8 + t66 * t7 + t65 * t15) * g(1);
U  = t2;
