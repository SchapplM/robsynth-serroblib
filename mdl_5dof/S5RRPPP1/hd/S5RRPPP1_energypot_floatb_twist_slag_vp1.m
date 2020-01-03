% Calculate potential energy for
% S5RRPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,theta3]';
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
% Datum: 2019-12-31 19:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPPP1_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPP1_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRPPP1_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPP1_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPP1_energypot_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPP1_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPP1_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:23:16
% EndTime: 2019-12-31 19:23:17
% DurationCPUTime: 0.36s
% Computational Cost: add. (186->94), mult. (367->107), div. (0->0), fcn. (406->8), ass. (0->44)
t30 = sin(qJ(2));
t29 = cos(pkin(5));
t31 = sin(qJ(1));
t52 = t31 * t29;
t27 = sin(pkin(5));
t33 = cos(qJ(1));
t57 = t27 * t33;
t61 = t30 * t52 + t57;
t47 = rSges(6,2) + qJ(4);
t45 = rSges(6,3) + qJ(5);
t60 = rSges(6,1) + pkin(4);
t59 = t27 * t31;
t32 = cos(qJ(2));
t58 = t27 * t32;
t56 = t29 * t32;
t26 = sin(pkin(8));
t55 = t30 * t26;
t28 = cos(pkin(8));
t54 = t30 * t28;
t53 = t30 * t31;
t51 = t31 * t32;
t50 = t32 * t33;
t49 = t33 * t29;
t48 = qJ(3) * t29;
t46 = rSges(5,3) + qJ(4);
t44 = pkin(6) + r_base(3);
t43 = t31 * pkin(1) + r_base(2);
t41 = qJ(3) * t27 * t30;
t40 = t30 * pkin(2) + t44;
t39 = t33 * pkin(1) + t31 * pkin(7) + r_base(1);
t9 = t26 * t56 + t54;
t38 = t9 * pkin(3) + t40;
t37 = pkin(2) * t50 + t31 * t48 + t33 * t41 + t39;
t6 = t26 * t59 + t28 * t50 - t49 * t55;
t36 = t6 * pkin(3) + t37;
t35 = t31 * t41 + pkin(2) * t51 + (-pkin(7) - t48) * t33 + t43;
t4 = -t61 * t26 + t28 * t51;
t34 = t4 * pkin(3) + t35;
t11 = t30 * t57 + t52;
t10 = t27 * t53 - t49;
t8 = -t28 * t56 + t55;
t5 = -t28 * t59 + (t26 * t32 + t29 * t54) * t33;
t3 = t26 * t51 + t61 * t28;
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t33 * rSges(2,1) - t31 * rSges(2,2) + r_base(1)) + g(2) * (t31 * rSges(2,1) + t33 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t44)) - m(3) * (g(1) * (t31 * rSges(3,3) + t39) + g(2) * (rSges(3,1) * t51 - rSges(3,2) * t53 + t43) + g(3) * (t30 * rSges(3,1) + t32 * rSges(3,2) + t44) + (g(1) * (rSges(3,1) * t32 - rSges(3,2) * t30) + g(2) * (-rSges(3,3) - pkin(7))) * t33) - m(4) * (g(1) * (t6 * rSges(4,1) - t5 * rSges(4,2) + t11 * rSges(4,3) + t37) + g(2) * (t4 * rSges(4,1) - t3 * rSges(4,2) + t10 * rSges(4,3) + t35) + g(3) * (t9 * rSges(4,1) - t8 * rSges(4,2) + (-rSges(4,3) - qJ(3)) * t58 + t40)) - m(5) * (g(1) * (t11 * rSges(5,1) - t6 * rSges(5,2) + t46 * t5 + t36) + g(2) * (t10 * rSges(5,1) - t4 * rSges(5,2) + t46 * t3 + t34) + g(3) * (-t9 * rSges(5,2) + t46 * t8 + (-rSges(5,1) - qJ(3)) * t58 + t38)) - m(6) * (g(1) * (t60 * t11 + t45 * t6 + t47 * t5 + t36) + g(2) * (t60 * t10 + t47 * t3 + t45 * t4 + t34) + (t38 + (-qJ(3) - t60) * t58 + t45 * t9 + t47 * t8) * g(3));
U = t1;
