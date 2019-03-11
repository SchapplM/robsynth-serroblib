% Calculate potential energy for
% S6RPRRPP8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4]';
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
% Datum: 2019-03-09 04:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRPP8_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP8_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPRRPP8_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP8_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP8_energypot_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP8_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPP8_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:53:48
% EndTime: 2019-03-09 04:53:48
% DurationCPUTime: 0.39s
% Computational Cost: add. (164->98), mult. (236->108), div. (0->0), fcn. (230->6), ass. (0->33)
t44 = rSges(7,1) + pkin(5);
t47 = rSges(7,3) + qJ(6);
t22 = sin(qJ(1));
t46 = g(1) * t22;
t25 = cos(qJ(1));
t45 = g(2) * t25;
t24 = cos(qJ(3));
t43 = rSges(4,2) * t24;
t21 = sin(qJ(3));
t42 = t21 * t22;
t20 = sin(qJ(4));
t41 = t22 * t20;
t23 = cos(qJ(4));
t40 = t22 * t23;
t39 = t25 * t20;
t38 = t25 * t23;
t37 = pkin(6) + r_base(3);
t36 = t22 * pkin(1) + r_base(2);
t35 = pkin(2) + t37;
t34 = t25 * pkin(1) + t22 * qJ(2) + r_base(1);
t33 = t22 * pkin(7) + t36;
t32 = t25 * pkin(7) + t34;
t31 = t24 * pkin(3) + t21 * pkin(8) + t35;
t30 = pkin(3) * t42 + t32;
t29 = t31 + (pkin(4) * t23 + qJ(5) * t20) * t24;
t3 = t21 * t41 - t38;
t4 = t21 * t40 + t39;
t28 = t4 * pkin(4) + t3 * qJ(5) + t30;
t27 = t33 + (-pkin(3) * t21 + t24 * pkin(8) - qJ(2)) * t25;
t5 = t21 * t39 + t40;
t6 = t21 * t38 - t41;
t26 = -t6 * pkin(4) - t5 * qJ(5) + t27;
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t25 * rSges(2,1) - t22 * rSges(2,2) + r_base(1)) + g(2) * (t22 * rSges(2,1) + t25 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t37)) - m(3) * (g(1) * (-t25 * rSges(3,2) + t22 * rSges(3,3) + t34) + g(2) * (-t22 * rSges(3,2) + (-rSges(3,3) - qJ(2)) * t25 + t36) + g(3) * (rSges(3,1) + t37)) - m(4) * (g(1) * (rSges(4,1) * t42 + t22 * t43 + t32) + g(2) * (t22 * rSges(4,3) + t33) + g(3) * (t24 * rSges(4,1) - t21 * rSges(4,2) + t35) + (g(1) * rSges(4,3) + g(2) * (-rSges(4,1) * t21 - qJ(2) - t43)) * t25) - m(5) * (g(1) * (t4 * rSges(5,1) - t3 * rSges(5,2) + t30) + g(2) * (-t6 * rSges(5,1) + t5 * rSges(5,2) + t27) + g(3) * (t21 * rSges(5,3) + t31) + (rSges(5,3) * t45 + g(3) * (rSges(5,1) * t23 - rSges(5,2) * t20) + (-rSges(5,3) - pkin(8)) * t46) * t24) - m(6) * (g(1) * (-t4 * rSges(6,2) + t3 * rSges(6,3) + t28) + g(2) * (t6 * rSges(6,2) - t5 * rSges(6,3) + t26) + g(3) * (t21 * rSges(6,1) + t29) + (rSges(6,1) * t45 + g(3) * (-rSges(6,2) * t23 + rSges(6,3) * t20) + (-rSges(6,1) - pkin(8)) * t46) * t24) - m(7) * (g(1) * (t3 * rSges(7,2) + t47 * t4 + t28) + g(2) * (-t5 * rSges(7,2) - t47 * t6 + t26) + g(3) * (t44 * t21 + t29) + (g(3) * (rSges(7,2) * t20 + t47 * t23) + t44 * t45 + (-pkin(8) - t44) * t46) * t24);
U  = t1;
