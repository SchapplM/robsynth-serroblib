% Calculate potential energy for
% S4RPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,theta2]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:46
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function U = S4RPPP1_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPP1_energypot_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S4RPPP1_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPP1_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPP1_energypot_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPP1_energypot_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPPP1_energypot_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:45:22
% EndTime: 2018-11-14 13:45:22
% DurationCPUTime: 0.28s
% Computational Cost: add. (207->74), mult. (208->77), div. (0->0), fcn. (177->10), ass. (0->34)
t46 = rSges(5,2) + qJ(3);
t45 = rSges(5,3) + qJ(4);
t44 = rSges(5,1) + pkin(3);
t21 = sin(pkin(4));
t24 = sin(qJ(1));
t43 = t21 * t24;
t25 = cos(qJ(1));
t42 = t21 * t25;
t41 = rSges(4,3) + qJ(3);
t40 = pkin(5) + r_base(3);
t39 = t24 * pkin(1) + r_base(2);
t38 = pkin(4) - pkin(6);
t37 = pkin(4) + pkin(6);
t23 = cos(pkin(4));
t36 = t23 * qJ(2) + t40;
t22 = cos(pkin(6));
t28 = sin(t37) / 0.2e1;
t32 = sin(t38);
t26 = t28 - t32 / 0.2e1;
t4 = t24 * t22 + t25 * t26;
t35 = t4 * pkin(2) + t39;
t34 = t25 * pkin(1) + qJ(2) * t43 + r_base(1);
t33 = cos(t37);
t29 = cos(t38) / 0.2e1;
t11 = t29 - t33 / 0.2e1;
t31 = t11 * pkin(2) + t36;
t6 = t25 * t22 - t24 * t26;
t30 = t6 * pkin(2) + t34;
t27 = t29 + t33 / 0.2e1;
t20 = sin(pkin(6));
t10 = t28 + t32 / 0.2e1;
t5 = t25 * t20 + t24 * t27;
t3 = t24 * t20 - t25 * t27;
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t25 * rSges(2,1) - t24 * rSges(2,2) + r_base(1)) + g(2) * (t24 * rSges(2,1) + t25 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t40)) - m(3) * (g(1) * (t6 * rSges(3,1) - t5 * rSges(3,2) + rSges(3,3) * t43 + t34) + g(2) * (t4 * rSges(3,1) - t3 * rSges(3,2) + (-rSges(3,3) - qJ(2)) * t42 + t39) + g(3) * (t11 * rSges(3,1) + t10 * rSges(3,2) + t23 * rSges(3,3) + t36)) - m(4) * (g(1) * (rSges(4,1) * t43 - t6 * rSges(4,2) + t41 * t5 + t30) + g(2) * (-t4 * rSges(4,2) + t41 * t3 + (-rSges(4,1) - qJ(2)) * t42 + t35) + g(3) * (t23 * rSges(4,1) - t11 * rSges(4,2) - t41 * t10 + t31)) - m(5) * (g(1) * (t45 * t6 + t46 * t5 + t30) + g(2) * (t46 * t3 + t45 * t4 + t35) + g(3) * (-t46 * t10 + t45 * t11 + t44 * t23 + t31) + (g(1) * t44 * t24 + g(2) * (-qJ(2) - t44) * t25) * t21);
U  = t1;
