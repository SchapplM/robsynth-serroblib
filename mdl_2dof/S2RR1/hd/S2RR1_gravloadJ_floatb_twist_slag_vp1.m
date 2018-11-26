% Calculate Gravitation load on the joints for
% S2RR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d2]';
% m_mdh [3x1]
%   mass of all robot links (including the base)
% rSges [3x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% taug [2x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-16 16:44
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S2RR1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(3,1),zeros(1,1),zeros(3,1),zeros(3,3)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S2RR1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S2RR1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [1x1] (double)');
assert( isreal(m) && all(size(m) == [3 1]), ...
  'S2RR1_gravloadJ_floatb_twist_slag_vp1: m has to be [3x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [3,3]), ...
  'S2RR1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [3x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-16 16:44:19
% EndTime: 2018-11-16 16:44:20
% DurationCPUTime: 0.08s
% Computational Cost: add. (17->13), mult. (36->23), div. (0->0), fcn. (26->4), ass. (0->7)
t7 = rSges(3,3) + pkin(1);
t1 = sin(qJ(2));
t3 = cos(qJ(2));
t6 = rSges(3,1) * t3 - rSges(3,2) * t1;
t4 = cos(qJ(1));
t2 = sin(qJ(1));
t5 = [-m(2) * (g(1) * (-rSges(2,1) * t4 + t2 * rSges(2,2)) + g(3) * (t2 * rSges(2,1) + rSges(2,2) * t4)) - m(3) * ((-g(1) * t6 + g(3) * t7) * t4 + (g(1) * t7 + g(3) * t6) * t2) -m(3) * (-g(2) * t6 + (g(1) * t2 + g(3) * t4) * (rSges(3,1) * t1 + rSges(3,2) * t3))];
taug  = t5(:);
