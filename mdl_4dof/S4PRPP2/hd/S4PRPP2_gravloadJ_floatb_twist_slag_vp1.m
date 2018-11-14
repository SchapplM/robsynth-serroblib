% Calculate Gravitation load on the joints for
% S4PRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,theta3]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 14:00
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S4PRPP2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPP2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRPP2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [5x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPP2_gravloadJ_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRPP2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:00:25
% EndTime: 2018-11-14 14:00:25
% DurationCPUTime: 0.09s
% Computational Cost: add. (39->23), mult. (38->30), div. (0->0), fcn. (22->4), ass. (0->11)
t5 = sin(qJ(2));
t10 = pkin(2) * t5;
t9 = -m(4) - m(5);
t8 = rSges(5,1) + pkin(3);
t7 = rSges(5,3) + qJ(4);
t6 = cos(qJ(2));
t4 = qJ(2) + pkin(5);
t3 = t6 * pkin(2);
t2 = cos(t4);
t1 = sin(t4);
t11 = [(-m(2) - m(3) + t9) * g(2), -m(3) * (g(1) * (-t5 * rSges(3,1) - rSges(3,2) * t6) + g(2) * (rSges(3,1) * t6 - t5 * rSges(3,2))) - m(4) * (g(1) * (-rSges(4,1) * t1 - rSges(4,2) * t2 - t10) + g(2) * (rSges(4,1) * t2 - rSges(4,2) * t1 + t3)) - m(5) * (-g(1) * t10 + g(2) * t3 + (g(1) * t7 + g(2) * t8) * t2 + (-g(1) * t8 + g(2) * t7) * t1) t9 * g(3), -m(5) * (g(1) * t1 - g(2) * t2)];
taug  = t11(:);
