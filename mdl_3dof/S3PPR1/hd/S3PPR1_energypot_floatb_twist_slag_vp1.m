% Calculate potential energy for
% S3PPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d3]';
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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S3PPR1_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,3)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PPR1_energypot_floatb_twist_slag_vp1: qJ has to be [3x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S3PPR1_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3PPR1_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S3PPR1_energypot_floatb_twist_slag_vp1: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'S3PPR1_energypot_floatb_twist_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'S3PPR1_energypot_floatb_twist_slag_vp1: rSges has to be [4x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:01:49
% EndTime: 2019-03-08 18:01:49
% DurationCPUTime: 0.08s
% Computational Cost: add. (34->30), mult. (20->20), div. (0->0), fcn. (4->2), ass. (0->6)
t5 = r_base(3) - pkin(1);
t4 = qJ(1) + r_base(2);
t3 = qJ(2) + r_base(1);
t2 = cos(qJ(3));
t1 = sin(qJ(3));
t6 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (r_base(1) - rSges(2,2)) + g(2) * (rSges(2,3) + t4) + g(3) * (r_base(3) - rSges(2,1))) - m(3) * (g(1) * (rSges(3,3) + t3) + g(2) * (rSges(3,1) + t4) + g(3) * (rSges(3,2) + t5)) - m(4) * (g(1) * (t1 * rSges(4,1) + rSges(4,2) * t2 + t3) + g(2) * (rSges(4,1) * t2 - t1 * rSges(4,2) + pkin(2) + t4) + g(3) * (-pkin(3) - rSges(4,3) + t5));
U  = t6;
