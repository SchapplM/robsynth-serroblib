% Calculate potential energy for
% S4PPPR4
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
%   pkin=[a2,a3,a4,d4,theta3]';
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
% Datum: 2018-11-14 14:04
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function U = S4PPPR4_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPPR4_energypot_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S4PPPR4_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPPR4_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPPR4_energypot_floatb_twist_slag_vp1: pkin has to be [5x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PPPR4_energypot_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PPPR4_energypot_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:04:10
% EndTime: 2018-11-14 14:04:10
% DurationCPUTime: 0.10s
% Computational Cost: add. (54->42), mult. (30->30), div. (0->0), fcn. (10->4), ass. (0->11)
t10 = r_base(3) - pkin(1);
t9 = qJ(1) + r_base(1);
t8 = r_base(2) - qJ(2);
t7 = pkin(2) + t9;
t6 = -qJ(3) + t10;
t5 = cos(pkin(5));
t4 = sin(pkin(5));
t3 = pkin(5) + qJ(4);
t2 = cos(t3);
t1 = sin(t3);
t11 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,3) + t9) + g(2) * (r_base(2) + rSges(2,2)) + g(3) * (r_base(3) - rSges(2,1))) - m(3) * (g(1) * (rSges(3,1) + t9) + g(2) * (-rSges(3,3) + t8) + g(3) * (rSges(3,2) + t10)) - m(4) * (g(1) * (rSges(4,1) * t5 - t4 * rSges(4,2) + t7) + g(2) * (-t4 * rSges(4,1) - rSges(4,2) * t5 + t8) + g(3) * (-rSges(4,3) + t6)) - m(5) * (g(1) * (t2 * rSges(5,1) - t1 * rSges(5,2) + pkin(3) * t5 + t7) + g(2) * (-rSges(5,1) * t1 - rSges(5,2) * t2 - pkin(3) * t4 + t8) + g(3) * (-pkin(4) - rSges(5,3) + t6));
U  = t11;
