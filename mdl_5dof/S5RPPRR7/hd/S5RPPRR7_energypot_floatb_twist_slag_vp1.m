% Calculate potential energy for
% S5RPPRR7
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
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
% Datum: 2019-12-31 18:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPPRR7_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR7_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RPPRR7_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR7_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR7_energypot_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR7_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR7_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:59:27
% EndTime: 2019-12-31 17:59:28
% DurationCPUTime: 0.31s
% Computational Cost: add. (141->73), mult. (107->76), div. (0->0), fcn. (87->8), ass. (0->24)
t33 = rSges(6,3) + pkin(7);
t12 = sin(qJ(4));
t15 = cos(qJ(4));
t32 = rSges(5,1) * t12 + rSges(5,2) * t15;
t31 = t12 * pkin(4);
t11 = sin(qJ(5));
t27 = t11 * t12;
t14 = cos(qJ(5));
t26 = t12 * t14;
t25 = pkin(5) + r_base(3);
t13 = sin(qJ(1));
t24 = t13 * pkin(1) + r_base(2);
t16 = cos(qJ(1));
t23 = t16 * pkin(1) + r_base(1);
t10 = qJ(1) + pkin(8);
t6 = sin(t10);
t22 = t6 * pkin(2) + t24;
t21 = qJ(2) + t25;
t7 = cos(t10);
t20 = t7 * pkin(2) + t6 * qJ(3) + t23;
t19 = t6 * pkin(6) + t22;
t18 = pkin(3) + t21;
t17 = t7 * pkin(6) + t20;
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t16 * rSges(2,1) - t13 * rSges(2,2) + r_base(1)) + g(2) * (t13 * rSges(2,1) + t16 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t25)) - m(3) * (g(1) * (t7 * rSges(3,1) - t6 * rSges(3,2) + t23) + g(2) * (t6 * rSges(3,1) + t7 * rSges(3,2) + t24) + g(3) * (rSges(3,3) + t21)) - m(4) * (g(1) * (-t7 * rSges(4,2) + t6 * rSges(4,3) + t20) + g(2) * (-t6 * rSges(4,2) + (-rSges(4,3) - qJ(3)) * t7 + t22) + g(3) * (rSges(4,1) + t21)) - m(5) * (g(1) * (t32 * t6 + t17) + g(2) * (t6 * rSges(5,3) + t19) + g(3) * (t15 * rSges(5,1) - t12 * rSges(5,2) + t18) + (g(1) * rSges(5,3) + g(2) * (-qJ(3) - t32)) * t7) - m(6) * (g(1) * (t6 * t31 + (t7 * t11 + t6 * t26) * rSges(6,1) + (t7 * t14 - t6 * t27) * rSges(6,2) + t17) + g(2) * (t19 + (t11 * rSges(6,1) + t14 * rSges(6,2)) * t6 + (-t26 * rSges(6,1) + t27 * rSges(6,2) - qJ(3) - t31) * t7) + g(3) * (t33 * t12 + t18) + (g(3) * (rSges(6,1) * t14 - rSges(6,2) * t11 + pkin(4)) + (-g(1) * t6 + g(2) * t7) * t33) * t15);
U = t1;
