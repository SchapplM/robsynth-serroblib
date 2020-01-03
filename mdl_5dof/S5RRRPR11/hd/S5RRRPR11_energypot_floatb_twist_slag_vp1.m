% Calculate potential energy for
% S5RRRPR11
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
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
% Datum: 2019-12-31 21:36
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRPR11_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR11_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRRPR11_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR11_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR11_energypot_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR11_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR11_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:32:20
% EndTime: 2019-12-31 21:32:20
% DurationCPUTime: 0.49s
% Computational Cost: add. (149->95), mult. (242->113), div. (0->0), fcn. (250->8), ass. (0->29)
t41 = -rSges(6,3) - pkin(8);
t20 = sin(qJ(2));
t21 = sin(qJ(1));
t40 = t21 * t20;
t24 = cos(qJ(2));
t39 = t21 * t24;
t19 = sin(qJ(3));
t25 = cos(qJ(1));
t38 = t25 * t19;
t37 = t25 * t20;
t23 = cos(qJ(3));
t36 = t25 * t23;
t35 = rSges(5,3) + qJ(4);
t34 = pkin(5) + r_base(3);
t33 = t21 * pkin(1) + r_base(2);
t32 = t20 * pkin(2) + t34;
t31 = t25 * pkin(1) + t21 * pkin(6) + r_base(1);
t30 = t32 + (pkin(3) * t23 + qJ(4) * t19) * t20;
t29 = t25 * t24 * pkin(2) + pkin(7) * t37 + t31;
t6 = t21 * t19 + t24 * t36;
t28 = t6 * pkin(3) + t29;
t27 = pkin(2) * t39 - t25 * pkin(6) + pkin(7) * t40 + t33;
t4 = t23 * t39 - t38;
t26 = t4 * pkin(3) + t27;
t22 = cos(qJ(5));
t18 = sin(qJ(5));
t5 = -t21 * t23 + t24 * t38;
t3 = t19 * t39 + t36;
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t25 * rSges(2,1) - t21 * rSges(2,2) + r_base(1)) + g(2) * (t21 * rSges(2,1) + t25 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t34)) - m(3) * (g(1) * (t21 * rSges(3,3) + t31) + g(2) * (rSges(3,1) * t39 - rSges(3,2) * t40 + t33) + g(3) * (t20 * rSges(3,1) + t24 * rSges(3,2) + t34) + (g(1) * (rSges(3,1) * t24 - rSges(3,2) * t20) + g(2) * (-rSges(3,3) - pkin(6))) * t25) - m(4) * (g(1) * (t6 * rSges(4,1) - t5 * rSges(4,2) + rSges(4,3) * t37 + t29) + g(2) * (t4 * rSges(4,1) - t3 * rSges(4,2) + rSges(4,3) * t40 + t27) + g(3) * ((-rSges(4,3) - pkin(7)) * t24 + (rSges(4,1) * t23 - rSges(4,2) * t19) * t20 + t32)) - m(5) * (g(1) * (t6 * rSges(5,1) + rSges(5,2) * t37 + t35 * t5 + t28) + g(2) * (t4 * rSges(5,1) + rSges(5,2) * t40 + t3 * t35 + t26) + g(3) * ((-rSges(5,2) - pkin(7)) * t24 + (rSges(5,1) * t23 + rSges(5,3) * t19) * t20 + t30)) - m(6) * (g(1) * (t6 * pkin(4) + t5 * qJ(4) + (t5 * t18 + t6 * t22) * rSges(6,1) + (-t6 * t18 + t5 * t22) * rSges(6,2) + t28) + g(2) * (t4 * pkin(4) + t3 * qJ(4) + (t3 * t18 + t4 * t22) * rSges(6,1) + (-t4 * t18 + t3 * t22) * rSges(6,2) + t26) + (g(1) * t25 + g(2) * t21) * t20 * t41 + (t30 + (-pkin(7) - t41) * t24 + (t23 * pkin(4) + (t18 * t19 + t22 * t23) * rSges(6,1) + (-t18 * t23 + t19 * t22) * rSges(6,2)) * t20) * g(3));
U = t1;
