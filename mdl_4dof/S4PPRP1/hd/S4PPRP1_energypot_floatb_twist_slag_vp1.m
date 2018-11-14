% Calculate potential energy for
% S4PPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,theta1]';
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
% Datum: 2018-11-14 13:39
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function U = S4PPRP1_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRP1_energypot_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S4PPRP1_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRP1_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPRP1_energypot_floatb_twist_slag_vp1: pkin has to be [5x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRP1_energypot_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PPRP1_energypot_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:39:02
% EndTime: 2018-11-14 13:39:02
% DurationCPUTime: 0.13s
% Computational Cost: add. (75->50), mult. (80->48), div. (0->0), fcn. (72->4), ass. (0->15)
t22 = -rSges(5,1) - pkin(3);
t21 = cos(qJ(3));
t20 = sin(qJ(3));
t19 = rSges(5,3) + qJ(4);
t11 = sin(pkin(5));
t18 = t11 * pkin(1) + r_base(2);
t17 = qJ(1) + r_base(3);
t12 = cos(pkin(5));
t16 = t12 * pkin(1) + t11 * qJ(2) + r_base(1);
t15 = -pkin(4) + t17;
t14 = t12 * pkin(2) + t16;
t13 = t11 * pkin(2) - t12 * qJ(2) + t18;
t2 = -t11 * t21 + t12 * t20;
t1 = -t11 * t20 - t12 * t21;
t3 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t12 - t11 * rSges(2,2) + r_base(1)) + g(2) * (t11 * rSges(2,1) + rSges(2,2) * t12 + r_base(2)) + g(3) * (rSges(2,3) + t17)) - m(3) * (g(1) * (rSges(3,1) * t12 + t11 * rSges(3,3) + t16) + g(2) * (t11 * rSges(3,1) + (-rSges(3,3) - qJ(2)) * t12 + t18) + g(3) * (rSges(3,2) + t17)) - m(4) * (g(1) * (-rSges(4,1) * t1 - rSges(4,2) * t2 + t14) + g(2) * (-t2 * rSges(4,1) + t1 * rSges(4,2) + t13) + g(3) * (-rSges(4,3) + t15)) - m(5) * (g(1) * t14 + g(2) * t13 + g(3) * (-rSges(5,2) + t15) + (g(1) * t19 + g(2) * t22) * t2 + (g(1) * t22 - g(2) * t19) * t1);
U  = t3;
