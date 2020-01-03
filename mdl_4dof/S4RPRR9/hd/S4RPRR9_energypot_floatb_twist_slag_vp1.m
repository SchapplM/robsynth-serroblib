% Calculate potential energy for
% S4RPRR9
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
%   pkin=[a2,a3,a4,d1,d3,d4]';
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
% Datum: 2019-12-31 16:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4RPRR9_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR9_energypot_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S4RPRR9_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR9_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR9_energypot_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR9_energypot_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRR9_energypot_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:56:05
% EndTime: 2019-12-31 16:56:06
% DurationCPUTime: 0.24s
% Computational Cost: add. (76->61), mult. (91->67), div. (0->0), fcn. (75->6), ass. (0->19)
t26 = rSges(5,3) + pkin(6);
t10 = cos(qJ(3));
t7 = sin(qJ(3));
t25 = rSges(4,1) * t7 + rSges(4,2) * t10;
t6 = sin(qJ(4));
t8 = sin(qJ(1));
t24 = t8 * t6;
t9 = cos(qJ(4));
t23 = t8 * t9;
t11 = cos(qJ(1));
t21 = t11 * t7;
t20 = t11 * t9;
t17 = pkin(4) + r_base(3);
t16 = t8 * pkin(1) + r_base(2);
t15 = pkin(2) + t17;
t14 = t11 * pkin(1) + t8 * qJ(2) + r_base(1);
t13 = t8 * pkin(5) + t16;
t12 = t11 * pkin(5) + t14;
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t11 - t8 * rSges(2,2) + r_base(1)) + g(2) * (t8 * rSges(2,1) + rSges(2,2) * t11 + r_base(2)) + g(3) * (rSges(2,3) + t17)) - m(3) * (g(1) * (-rSges(3,2) * t11 + t8 * rSges(3,3) + t14) + g(2) * (-t8 * rSges(3,2) + (-rSges(3,3) - qJ(2)) * t11 + t16) + g(3) * (rSges(3,1) + t17)) - m(4) * (g(1) * (t25 * t8 + t12) + g(2) * (t8 * rSges(4,3) + t13) + g(3) * (rSges(4,1) * t10 - rSges(4,2) * t7 + t15) + (g(1) * rSges(4,3) + g(2) * (-qJ(2) - t25)) * t11) - m(5) * (g(1) * (t8 * t7 * pkin(3) + (t11 * t6 + t7 * t23) * rSges(5,1) + (-t7 * t24 + t20) * rSges(5,2) + t12) + g(2) * (-pkin(3) * t21 - t11 * qJ(2) + (-t7 * t20 + t24) * rSges(5,1) + (t6 * t21 + t23) * rSges(5,2) + t13) + g(3) * (t26 * t7 + t15) + (g(3) * (rSges(5,1) * t9 - rSges(5,2) * t6 + pkin(3)) + (-g(1) * t8 + g(2) * t11) * t26) * t10);
U = t1;
