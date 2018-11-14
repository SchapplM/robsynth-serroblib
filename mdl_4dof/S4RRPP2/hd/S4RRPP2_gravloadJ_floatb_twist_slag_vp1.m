% Calculate Gravitation load on the joints for
% S4RRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
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
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:53
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S4RRPP2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [5x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP2_gravloadJ_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPP2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:52:29
% EndTime: 2018-11-14 13:52:29
% DurationCPUTime: 0.11s
% Computational Cost: add. (103->38), mult. (82->46), div. (0->0), fcn. (58->4), ass. (0->19)
t29 = rSges(5,1) + pkin(3);
t16 = qJ(1) + qJ(2);
t14 = cos(t16);
t3 = t14 * qJ(3);
t28 = t14 * rSges(5,2) + t3;
t13 = sin(t16);
t27 = g(1) * t13;
t17 = sin(qJ(1));
t26 = t17 * pkin(1);
t25 = t14 * pkin(2) + t13 * qJ(3);
t24 = t14 * rSges(4,1) + t13 * rSges(4,3) + t25;
t23 = t14 * rSges(3,1) - t13 * rSges(3,2);
t22 = t13 * rSges(5,2) + t29 * t14 + t25;
t21 = (-pkin(2) - t29) * t27;
t20 = -t13 * rSges(3,1) - t14 * rSges(3,2);
t19 = t3 + t14 * rSges(4,3) + (-rSges(4,1) - pkin(2)) * t13;
t18 = cos(qJ(1));
t15 = t18 * pkin(1);
t1 = [-m(2) * (g(1) * (-t17 * rSges(2,1) - t18 * rSges(2,2)) + g(2) * (t18 * rSges(2,1) - t17 * rSges(2,2))) - m(3) * (g(1) * (t20 - t26) + g(2) * (t15 + t23)) - m(4) * (g(1) * (t19 - t26) + g(2) * (t15 + t24)) - m(5) * (g(1) * (-t26 + t28) + g(2) * (t15 + t22) + t21) -m(3) * (g(1) * t20 + g(2) * t23) - m(4) * (g(1) * t19 + g(2) * t24) - m(5) * (g(1) * t28 + g(2) * t22 + t21) (-m(4) - m(5)) * (-g(2) * t14 + t27) m(5) * g(3)];
taug  = t1(:);
