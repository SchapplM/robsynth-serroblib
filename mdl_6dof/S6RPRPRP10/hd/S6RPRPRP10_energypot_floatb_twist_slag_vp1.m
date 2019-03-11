% Calculate potential energy for
% S6RPRPRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5]';
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
% Datum: 2019-03-09 03:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRPRP10_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP10_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPRPRP10_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP10_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPRP10_energypot_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP10_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRP10_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:30:52
% EndTime: 2019-03-09 03:30:52
% DurationCPUTime: 0.38s
% Computational Cost: add. (151->97), mult. (202->107), div. (0->0), fcn. (186->6), ass. (0->34)
t47 = rSges(7,1) + pkin(5);
t46 = rSges(7,3) + qJ(6);
t45 = -pkin(3) - pkin(8);
t21 = sin(qJ(1));
t44 = g(1) * t21;
t24 = cos(qJ(1));
t43 = g(2) * t24;
t20 = sin(qJ(3));
t42 = t20 * t21;
t23 = cos(qJ(3));
t41 = t21 * t23;
t19 = sin(qJ(5));
t40 = t24 * t19;
t22 = cos(qJ(5));
t39 = t24 * t22;
t38 = qJ(4) * t23;
t37 = pkin(6) + r_base(3);
t36 = t21 * pkin(1) + r_base(2);
t35 = pkin(2) + t37;
t34 = t24 * pkin(1) + t21 * qJ(2) + r_base(1);
t33 = t21 * pkin(7) + t36;
t32 = t24 * t38 + t33;
t31 = t24 * pkin(7) + t34;
t30 = t23 * pkin(3) + t20 * qJ(4) + t35;
t29 = -rSges(5,2) * t20 - rSges(5,3) * t23;
t28 = pkin(3) * t42 + t31;
t27 = t23 * pkin(8) + t30;
t26 = t21 * pkin(4) - t24 * qJ(2) + t32;
t25 = t24 * pkin(4) + pkin(8) * t42 - t21 * t38 + t28;
t4 = -t19 * t41 + t39;
t3 = t22 * t41 + t40;
t2 = t21 * t22 + t23 * t40;
t1 = t21 * t19 - t23 * t39;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t24 * rSges(2,1) - t21 * rSges(2,2) + r_base(1)) + g(2) * (t21 * rSges(2,1) + t24 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t37)) - m(3) * (g(1) * (-t24 * rSges(3,2) + t21 * rSges(3,3) + t34) + g(2) * (-t21 * rSges(3,2) + (-rSges(3,3) - qJ(2)) * t24 + t36) + g(3) * (rSges(3,1) + t37)) - m(4) * (g(1) * (rSges(4,1) * t42 + rSges(4,2) * t41 + t31) + g(2) * (t21 * rSges(4,3) + t33) + g(3) * (t23 * rSges(4,1) - t20 * rSges(4,2) + t35) + (g(1) * rSges(4,3) + g(2) * (-rSges(4,1) * t20 - rSges(4,2) * t23 - qJ(2))) * t24) - m(5) * (g(1) * t28 + g(2) * t32 + g(3) * (-t23 * rSges(5,2) + t20 * rSges(5,3) + t30) + (g(1) * (t29 - t38) + g(2) * rSges(5,1)) * t21 + (g(1) * rSges(5,1) + g(2) * (-pkin(3) * t20 - qJ(2) - t29)) * t24) - m(6) * (g(1) * (t4 * rSges(6,1) - t3 * rSges(6,2) + t25) + g(2) * (t2 * rSges(6,1) - t1 * rSges(6,2) + t26) + g(3) * (t23 * rSges(6,3) + t27) + (rSges(6,3) * t44 + g(3) * (rSges(6,1) * t19 + rSges(6,2) * t22) + (-rSges(6,3) + t45) * t43) * t20) - m(7) * (g(1) * (t46 * t3 + t47 * t4 + t25) + g(2) * (t46 * t1 + t47 * t2 + t26) + g(3) * (t23 * rSges(7,2) + t27) + (rSges(7,2) * t44 + g(3) * (t47 * t19 - t46 * t22) + (-rSges(7,2) + t45) * t43) * t20);
U  = t5;
