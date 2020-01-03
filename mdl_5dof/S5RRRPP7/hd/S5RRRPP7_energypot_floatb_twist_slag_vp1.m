% Calculate potential energy for
% S5RRRPP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
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
% Datum: 2019-12-31 21:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRPP7_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP7_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRRPP7_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP7_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP7_energypot_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP7_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPP7_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:03:52
% EndTime: 2019-12-31 21:03:53
% DurationCPUTime: 0.39s
% Computational Cost: add. (139->86), mult. (216->96), div. (0->0), fcn. (214->6), ass. (0->29)
t41 = rSges(6,1) + pkin(4);
t40 = rSges(6,2) + qJ(4);
t19 = sin(qJ(2));
t20 = sin(qJ(1));
t39 = t19 * t20;
t23 = cos(qJ(1));
t38 = t19 * t23;
t22 = cos(qJ(2));
t37 = t20 * t22;
t18 = sin(qJ(3));
t36 = t23 * t18;
t21 = cos(qJ(3));
t35 = t23 * t21;
t34 = rSges(5,3) + qJ(4);
t33 = -rSges(6,3) - qJ(5);
t32 = pkin(5) + r_base(3);
t31 = t20 * pkin(1) + r_base(2);
t30 = t19 * pkin(2) + t32;
t29 = t23 * pkin(1) + t20 * pkin(6) + r_base(1);
t28 = t30 + (pkin(3) * t21 + qJ(4) * t18) * t19;
t27 = t23 * t22 * pkin(2) + pkin(7) * t38 + t29;
t6 = t20 * t18 + t22 * t35;
t26 = t6 * pkin(3) + t27;
t25 = pkin(2) * t37 - t23 * pkin(6) + pkin(7) * t39 + t31;
t4 = t21 * t37 - t36;
t24 = t4 * pkin(3) + t25;
t5 = -t20 * t21 + t22 * t36;
t3 = t18 * t37 + t35;
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t23 * rSges(2,1) - t20 * rSges(2,2) + r_base(1)) + g(2) * (t20 * rSges(2,1) + t23 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t32)) - m(3) * (g(1) * (t20 * rSges(3,3) + t29) + g(2) * (rSges(3,1) * t37 - rSges(3,2) * t39 + t31) + g(3) * (t19 * rSges(3,1) + t22 * rSges(3,2) + t32) + (g(1) * (rSges(3,1) * t22 - rSges(3,2) * t19) + g(2) * (-rSges(3,3) - pkin(6))) * t23) - m(4) * (g(1) * (t6 * rSges(4,1) - t5 * rSges(4,2) + rSges(4,3) * t38 + t27) + g(2) * (t4 * rSges(4,1) - t3 * rSges(4,2) + rSges(4,3) * t39 + t25) + g(3) * ((-rSges(4,3) - pkin(7)) * t22 + (rSges(4,1) * t21 - rSges(4,2) * t18) * t19 + t30)) - m(5) * (g(1) * (t6 * rSges(5,1) + rSges(5,2) * t38 + t34 * t5 + t26) + g(2) * (t4 * rSges(5,1) + rSges(5,2) * t39 + t3 * t34 + t24) + g(3) * ((-rSges(5,2) - pkin(7)) * t22 + (rSges(5,1) * t21 + rSges(5,3) * t18) * t19 + t28)) - m(6) * (g(1) * (t40 * t5 + t41 * t6 + t26) + g(2) * (t40 * t3 + t41 * t4 + t24) + (g(1) * t23 + g(2) * t20) * t19 * t33 + (t28 + (-pkin(7) - t33) * t22 + (rSges(6,2) * t18 + t41 * t21) * t19) * g(3));
U = t1;
