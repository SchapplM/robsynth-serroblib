% Calculate potential energy for
% S4PPRR5
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
%   pkin=[a2,a3,a4,d3,d4,theta1]';
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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4PPRR5_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR5_energypot_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S4PPRR5_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR5_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR5_energypot_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR5_energypot_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PPRR5_energypot_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:19:38
% EndTime: 2019-12-31 16:19:38
% DurationCPUTime: 0.24s
% Computational Cost: add. (76->61), mult. (91->66), div. (0->0), fcn. (75->6), ass. (0->18)
t25 = rSges(5,3) + pkin(5);
t11 = cos(qJ(3));
t9 = sin(qJ(3));
t24 = rSges(4,1) * t9 + rSges(4,2) * t11;
t8 = sin(qJ(4));
t23 = t8 * t9;
t22 = t9 * pkin(3);
t10 = cos(qJ(4));
t20 = t10 * t9;
t6 = sin(pkin(6));
t17 = t6 * pkin(1) + r_base(2);
t16 = qJ(1) + r_base(3);
t7 = cos(pkin(6));
t15 = t7 * pkin(1) + t6 * qJ(2) + r_base(1);
t14 = t6 * pkin(4) + t17;
t13 = pkin(2) + t16;
t12 = t7 * pkin(4) + t15;
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t7 - rSges(2,2) * t6 + r_base(1)) + g(2) * (rSges(2,1) * t6 + rSges(2,2) * t7 + r_base(2)) + g(3) * (rSges(2,3) + t16)) - m(3) * (g(1) * (-rSges(3,2) * t7 + rSges(3,3) * t6 + t15) + g(2) * (-rSges(3,2) * t6 + (-rSges(3,3) - qJ(2)) * t7 + t17) + g(3) * (rSges(3,1) + t16)) - m(4) * (g(1) * (t24 * t6 + t12) + g(2) * (t6 * rSges(4,3) + t14) + g(3) * (rSges(4,1) * t11 - t9 * rSges(4,2) + t13) + (g(1) * rSges(4,3) + g(2) * (-qJ(2) - t24)) * t7) - m(5) * (g(1) * (t6 * t22 + (t6 * t20 + t7 * t8) * rSges(5,1) + (t10 * t7 - t6 * t23) * rSges(5,2) + t12) + g(2) * (t14 + (t8 * rSges(5,1) + t10 * rSges(5,2)) * t6 + (-t20 * rSges(5,1) + t23 * rSges(5,2) - qJ(2) - t22) * t7) + g(3) * (t25 * t9 + t13) + (g(3) * (rSges(5,1) * t10 - rSges(5,2) * t8 + pkin(3)) + (-g(1) * t6 + g(2) * t7) * t25) * t11);
U = t1;
