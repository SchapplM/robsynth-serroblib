% Calculate potential energy for
% S4PRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4]';
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
% Datum: 2019-07-18 13:27
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4PRRR2_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(2,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR2_energypot_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S4PRRR2_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR2_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S4PRRR2_energypot_floatb_twist_slag_vp1: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR2_energypot_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRR2_energypot_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:27:20
% EndTime: 2019-07-18 13:27:20
% DurationCPUTime: 0.11s
% Computational Cost: add. (59->42), mult. (38->36), div. (0->0), fcn. (18->6), ass. (0->12)
t7 = qJ(2) + qJ(3);
t9 = cos(qJ(2));
t12 = t9 * pkin(1) + r_base(3);
t11 = r_base(2) - qJ(1);
t8 = sin(qJ(2));
t10 = -t8 * pkin(1) + r_base(1);
t5 = qJ(4) + t7;
t4 = cos(t7);
t3 = sin(t7);
t2 = cos(t5);
t1 = sin(t5);
t6 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (r_base(1) - rSges(2,2)) + g(2) * (-rSges(2,3) + t11) + g(3) * (r_base(3) + rSges(2,1))) - m(3) * (g(1) * (-t8 * rSges(3,1) - t9 * rSges(3,2) + r_base(1)) + g(2) * (-rSges(3,3) + t11) + g(3) * (t9 * rSges(3,1) - t8 * rSges(3,2) + r_base(3))) - m(4) * (g(1) * (-t3 * rSges(4,1) - t4 * rSges(4,2) + t10) + g(2) * (-rSges(4,3) + t11) + g(3) * (t4 * rSges(4,1) - t3 * rSges(4,2) + t12)) - m(5) * (g(1) * (-t1 * rSges(5,1) - t2 * rSges(5,2) - pkin(2) * t3 + t10) + g(2) * (-rSges(5,3) + t11) + g(3) * (t2 * rSges(5,1) - t1 * rSges(5,2) + pkin(2) * t4 + t12));
U  = t6;
