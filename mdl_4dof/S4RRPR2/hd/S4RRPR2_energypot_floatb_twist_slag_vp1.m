% Calculate potential energy for
% S4RRPR2
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
%   pkin=[a2,a3,a4,d1,d2]';
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
% Datum: 2019-07-18 18:16
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4RRPR2_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR2_energypot_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S4RRPR2_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR2_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPR2_energypot_floatb_twist_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR2_energypot_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPR2_energypot_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 18:16:28
% EndTime: 2019-07-18 18:16:28
% DurationCPUTime: 0.12s
% Computational Cost: add. (90->51), mult. (60->48), div. (0->0), fcn. (44->6), ass. (0->16)
t20 = pkin(4) + r_base(3);
t12 = sin(qJ(1));
t19 = t12 * pkin(1) + r_base(2);
t14 = cos(qJ(1));
t18 = t14 * pkin(1) + r_base(1);
t17 = pkin(5) + t20;
t10 = qJ(1) + qJ(2);
t6 = sin(t10);
t16 = t6 * pkin(2) + t19;
t7 = cos(t10);
t15 = t7 * pkin(2) + t6 * qJ(3) + t18;
t13 = cos(qJ(4));
t11 = sin(qJ(4));
t2 = -t11 * t7 + t13 * t6;
t1 = -t11 * t6 - t13 * t7;
t3 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t14 - t12 * rSges(2,2) + r_base(1)) + g(2) * (t12 * rSges(2,1) + rSges(2,2) * t14 + r_base(2)) + g(3) * (rSges(2,3) + t20)) - m(3) * (g(1) * (rSges(3,1) * t7 - rSges(3,2) * t6 + t18) + g(2) * (rSges(3,1) * t6 + rSges(3,2) * t7 + t19) + g(3) * (rSges(3,3) + t17)) - m(4) * (g(1) * (rSges(4,1) * t7 + rSges(4,3) * t6 + t15) + g(2) * (rSges(4,1) * t6 + (-rSges(4,3) - qJ(3)) * t7 + t16) + g(3) * (rSges(4,2) + t17)) - m(5) * (g(1) * (-t1 * rSges(5,1) + t2 * rSges(5,2) + t7 * pkin(3) + t15) + g(2) * (t2 * rSges(5,1) + t1 * rSges(5,2) + t6 * pkin(3) - t7 * qJ(3) + t16) + g(3) * (-rSges(5,3) + t17));
U  = t3;
