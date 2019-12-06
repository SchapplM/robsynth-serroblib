% Calculate potential energy for
% S5PRRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRRRP8_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP8_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5PRRRP8_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP8_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP8_energypot_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP8_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP8_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:58:08
% EndTime: 2019-12-05 16:58:08
% DurationCPUTime: 0.52s
% Computational Cost: add. (240->84), mult. (493->102), div. (0->0), fcn. (575->10), ass. (0->42)
t65 = -m(1) - m(2);
t64 = -m(5) - m(6);
t63 = mrSges(3,2) - mrSges(4,3);
t62 = mrSges(4,2) - mrSges(5,3) - mrSges(6,2);
t61 = -m(6) * pkin(4) - mrSges(5,1) - mrSges(6,1);
t60 = -m(6) * qJ(5) + mrSges(5,2) - mrSges(6,3);
t32 = sin(pkin(9));
t33 = sin(pkin(5));
t59 = t32 * t33;
t34 = cos(pkin(9));
t58 = t33 * t34;
t37 = sin(qJ(3));
t57 = t33 * t37;
t38 = sin(qJ(2));
t56 = t33 * t38;
t40 = cos(qJ(3));
t55 = t33 * t40;
t41 = cos(qJ(2));
t54 = t33 * t41;
t35 = cos(pkin(5));
t53 = t35 * t38;
t52 = t35 * t41;
t51 = qJ(1) + r_base(3);
t50 = t34 * pkin(1) + pkin(6) * t59 + r_base(1);
t49 = t35 * pkin(6) + t51;
t48 = t32 * pkin(1) - pkin(6) * t58 + r_base(2);
t20 = t32 * t52 + t34 * t38;
t21 = -t32 * t53 + t34 * t41;
t47 = t21 * pkin(2) + t20 * pkin(7) + t50;
t46 = pkin(2) * t56 - pkin(7) * t54 + t49;
t18 = t32 * t38 - t34 * t52;
t19 = t32 * t41 + t34 * t53;
t45 = t19 * pkin(2) + t18 * pkin(7) + t48;
t39 = cos(qJ(4));
t36 = sin(qJ(4));
t23 = t35 * t37 + t38 * t55;
t22 = -t35 * t40 + t37 * t56;
t10 = t21 * t40 + t32 * t57;
t9 = t21 * t37 - t32 * t55;
t8 = t19 * t40 - t34 * t57;
t7 = t19 * t37 + t34 * t55;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t51 - mrSges(2,3) - m(3) * t49 - t35 * mrSges(3,3) - (t38 * mrSges(3,1) + t41 * mrSges(3,2)) * t33 - m(4) * t46 - t23 * mrSges(4,1) + mrSges(4,3) * t54 + t64 * (t23 * pkin(3) + t22 * pkin(8) + t46) + t61 * (t23 * t39 - t36 * t54) + t60 * (t23 * t36 + t39 * t54) + t62 * t22) * g(3) + (-m(3) * t48 - m(4) * t45 - t32 * mrSges(2,1) - t19 * mrSges(3,1) - t8 * mrSges(4,1) - t34 * mrSges(2,2) + mrSges(3,3) * t58 - mrSges(1,2) + t65 * r_base(2) + t64 * (t8 * pkin(3) + pkin(8) * t7 + t45) + t61 * (t18 * t36 + t8 * t39) + t63 * t18 + t60 * (-t18 * t39 + t8 * t36) + t62 * t7) * g(2) + (-m(3) * t50 - m(4) * t47 - t34 * mrSges(2,1) - t21 * mrSges(3,1) - t10 * mrSges(4,1) + t32 * mrSges(2,2) - mrSges(3,3) * t59 - mrSges(1,1) + t65 * r_base(1) + t64 * (t10 * pkin(3) + pkin(8) * t9 + t47) + t61 * (t10 * t39 + t20 * t36) + t60 * (t10 * t36 - t20 * t39) + t63 * t20 + t62 * t9) * g(1);
U = t1;
