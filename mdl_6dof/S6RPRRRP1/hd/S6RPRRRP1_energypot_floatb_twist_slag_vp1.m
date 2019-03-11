% Calculate potential energy for
% S6RPRRRP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRRP1_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP1_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPRRRP1_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP1_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP1_energypot_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP1_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRP1_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:55:20
% EndTime: 2019-03-09 05:55:21
% DurationCPUTime: 0.32s
% Computational Cost: add. (259->95), mult. (185->105), div. (0->0), fcn. (169->10), ass. (0->38)
t49 = rSges(7,1) + pkin(5);
t48 = rSges(7,3) + qJ(6);
t47 = rSges(4,3) + pkin(7);
t22 = qJ(1) + pkin(10);
t15 = sin(t22);
t23 = qJ(3) + qJ(4);
t17 = sin(t23);
t46 = t15 * t17;
t16 = cos(t22);
t45 = t16 * t17;
t18 = cos(t23);
t44 = t16 * t18;
t24 = sin(qJ(5));
t43 = t18 * t24;
t27 = cos(qJ(5));
t42 = t18 * t27;
t41 = pkin(6) + r_base(3);
t26 = sin(qJ(1));
t40 = t26 * pkin(1) + r_base(2);
t29 = cos(qJ(1));
t39 = t29 * pkin(1) + r_base(1);
t28 = cos(qJ(3));
t14 = t28 * pkin(3) + pkin(2);
t38 = t16 * t14 + t39;
t37 = qJ(2) + t41;
t30 = -pkin(8) - pkin(7);
t36 = t15 * t14 + t16 * t30 + t40;
t25 = sin(qJ(3));
t35 = t25 * pkin(3) + t37;
t34 = t17 * pkin(4) + t35;
t33 = rSges(4,1) * t28 - rSges(4,2) * t25 + pkin(2);
t32 = t15 * t18 * pkin(4) + pkin(9) * t46 + t36;
t31 = pkin(4) * t44 + pkin(9) * t45 - t15 * t30 + t38;
t4 = t15 * t24 + t16 * t42;
t3 = -t15 * t27 + t16 * t43;
t2 = t15 * t42 - t16 * t24;
t1 = t15 * t43 + t16 * t27;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t29 * rSges(2,1) - t26 * rSges(2,2) + r_base(1)) + g(2) * (t26 * rSges(2,1) + t29 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t41)) - m(3) * (g(1) * (t16 * rSges(3,1) - t15 * rSges(3,2) + t39) + g(2) * (t15 * rSges(3,1) + t16 * rSges(3,2) + t40) + g(3) * (rSges(3,3) + t37)) - m(4) * (g(1) * t39 + g(2) * t40 + g(3) * (t25 * rSges(4,1) + t28 * rSges(4,2) + t37) + (g(1) * t33 - g(2) * t47) * t16 + (g(1) * t47 + g(2) * t33) * t15) - m(5) * (g(1) * (rSges(5,1) * t44 - rSges(5,2) * t45 + t38) + g(2) * (-t16 * rSges(5,3) + t36) + g(3) * (t17 * rSges(5,1) + t18 * rSges(5,2) + t35) + (g(1) * (rSges(5,3) - t30) + g(2) * (rSges(5,1) * t18 - rSges(5,2) * t17)) * t15) - m(6) * (g(1) * (t4 * rSges(6,1) - t3 * rSges(6,2) + rSges(6,3) * t45 + t31) + g(2) * (t2 * rSges(6,1) - t1 * rSges(6,2) + rSges(6,3) * t46 + t32) + g(3) * ((-rSges(6,3) - pkin(9)) * t18 + (rSges(6,1) * t27 - rSges(6,2) * t24) * t17 + t34)) - m(7) * (g(1) * (t48 * t3 + t49 * t4 + t31) + g(2) * (t48 * t1 + t49 * t2 + t32) + g(3) * (t34 + (-rSges(7,2) - pkin(9)) * t18) + (g(3) * (t48 * t24 + t49 * t27) + (g(1) * t16 + g(2) * t15) * rSges(7,2)) * t17);
U  = t5;
