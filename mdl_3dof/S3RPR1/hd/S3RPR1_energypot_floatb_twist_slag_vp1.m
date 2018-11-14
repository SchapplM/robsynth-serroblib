% Calculate potential energy for
% S3RPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d3]';
% m_mdh [4x1]
%   mass of all robot links (including the base)
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 10:14
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function U = S3RPR1_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,1),zeros(4,3)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RPR1_energypot_floatb_twist_slag_vp1: qJ has to be [3x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S3RPR1_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3RPR1_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3RPR1_energypot_floatb_twist_slag_vp1: pkin has to be [4x1] (double)');
assert( isreal(m) && all(size(m) == [4 1]), ...
  'S3RPR1_energypot_floatb_twist_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'S3RPR1_energypot_floatb_twist_slag_vp1: rSges has to be [4x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:14:21
% EndTime: 2018-11-14 10:14:21
% DurationCPUTime: 0.10s
% Computational Cost: add. (47->40), mult. (46->38), div. (0->0), fcn. (34->4), ass. (0->10)
t12 = pkin(3) + r_base(3);
t7 = sin(qJ(1));
t11 = t7 * pkin(1) + r_base(2);
t9 = cos(qJ(1));
t10 = t9 * pkin(1) + t7 * qJ(2) + r_base(1);
t8 = cos(qJ(3));
t6 = sin(qJ(3));
t2 = -t6 * t9 + t7 * t8;
t1 = -t7 * t6 - t8 * t9;
t3 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t9 * rSges(2,1) - t7 * rSges(2,2) + r_base(1)) + g(2) * (t7 * rSges(2,1) + t9 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t12)) - m(3) * (g(1) * (t9 * rSges(3,1) + t7 * rSges(3,3) + t10) + g(2) * (t7 * rSges(3,1) + (-rSges(3,3) - qJ(2)) * t9 + t11) + g(3) * (rSges(3,2) + t12)) - m(4) * (g(1) * (-t1 * rSges(4,1) + t2 * rSges(4,2) + t9 * pkin(2) + t10) + g(2) * (t2 * rSges(4,1) + t1 * rSges(4,2) + t7 * pkin(2) - t9 * qJ(2) + t11) + g(3) * (-pkin(4) - rSges(4,3) + t12));
U  = t3;
