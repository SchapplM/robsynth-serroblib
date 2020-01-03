% Calculate potential energy for
% S5RRPRR16
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
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPRR16_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR16_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRPRR16_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR16_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR16_energypot_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR16_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR16_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:44:37
% EndTime: 2019-12-31 20:44:37
% DurationCPUTime: 0.40s
% Computational Cost: add. (193->105), mult. (358->130), div. (0->0), fcn. (405->10), ass. (0->44)
t29 = sin(qJ(1));
t56 = g(1) * t29;
t33 = cos(qJ(1));
t55 = g(2) * t33;
t54 = pkin(9) + rSges(6,3);
t24 = sin(pkin(5));
t28 = sin(qJ(2));
t53 = t24 * t28;
t52 = t24 * t29;
t32 = cos(qJ(2));
t51 = t24 * t32;
t50 = t24 * t33;
t49 = t29 * t28;
t48 = t29 * t32;
t47 = t33 * t28;
t46 = t33 * t32;
t45 = t32 * qJ(3);
t44 = pkin(6) + r_base(3);
t43 = t29 * pkin(1) + r_base(2);
t42 = (-pkin(3) - pkin(7)) * t33;
t25 = cos(pkin(5));
t41 = t25 * pkin(7) + t44;
t40 = t33 * pkin(1) + pkin(7) * t52 + r_base(1);
t39 = pkin(2) * t53 + t41;
t10 = -t25 * t46 + t49;
t11 = t25 * t47 + t48;
t38 = t11 * pkin(2) + t10 * qJ(3) + t43;
t37 = t25 * pkin(3) + pkin(8) * t53 + t39;
t12 = t25 * t48 + t47;
t13 = -t25 * t49 + t46;
t36 = t13 * pkin(2) + t12 * qJ(3) + t40;
t35 = pkin(3) * t52 + t36;
t34 = t11 * pkin(8) + t38;
t31 = cos(qJ(4));
t30 = cos(qJ(5));
t27 = sin(qJ(4));
t26 = sin(qJ(5));
t9 = t25 * t31 - t27 * t51;
t8 = t25 * t27 + t31 * t51;
t4 = t10 * t27 - t31 * t50;
t3 = t10 * t31 + t27 * t50;
t2 = t12 * t27 + t31 * t52;
t1 = -t12 * t31 + t27 * t52;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t33 * rSges(2,1) - t29 * rSges(2,2) + r_base(1)) + g(2) * (t29 * rSges(2,1) + t33 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t44)) - m(3) * (g(1) * (t13 * rSges(3,1) - t12 * rSges(3,2) + t40) + g(2) * (t11 * rSges(3,1) - t10 * rSges(3,2) + t43) + g(3) * (t25 * rSges(3,3) + t41) + (rSges(3,3) * t56 + g(3) * (rSges(3,1) * t28 + rSges(3,2) * t32) + (-rSges(3,3) - pkin(7)) * t55) * t24) - m(4) * (g(1) * (-t13 * rSges(4,2) + t12 * rSges(4,3) + t36) + g(2) * (-t11 * rSges(4,2) + t10 * rSges(4,3) + t38) + g(3) * (t25 * rSges(4,1) + t39) + (rSges(4,1) * t56 + g(3) * (-rSges(4,2) * t28 - rSges(4,3) * t32 - t45) + (-rSges(4,1) - pkin(7)) * t55) * t24) - m(5) * (g(1) * (t2 * rSges(5,1) - t1 * rSges(5,2) + (rSges(5,3) + pkin(8)) * t13 + t35) + g(2) * (t4 * rSges(5,1) + t3 * rSges(5,2) + t11 * rSges(5,3) + t34) + g(3) * (t9 * rSges(5,1) - t8 * rSges(5,2) + t37) + (g(3) * (rSges(5,3) * t28 - t45) + g(2) * t42) * t24) - m(6) * (g(1) * (t2 * pkin(4) + t13 * pkin(8) + (t13 * t26 + t2 * t30) * rSges(6,1) + (t13 * t30 - t2 * t26) * rSges(6,2) + t54 * t1 + t35) + g(2) * (t4 * pkin(4) + (t11 * t26 + t4 * t30) * rSges(6,1) + (t11 * t30 - t4 * t26) * rSges(6,2) - t54 * t3 + t24 * t42 + t34) + g(3) * (t9 * pkin(4) - t24 * t45 + (t26 * t53 + t9 * t30) * rSges(6,1) + (-t9 * t26 + t30 * t53) * rSges(6,2) + t54 * t8 + t37));
U = t5;
