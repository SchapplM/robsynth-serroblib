% Calculate Gravitation load on the joints for
% S4PRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
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
% Datum: 2018-11-14 13:45
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S4PRRR1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR1_gravloadJ_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRR1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:44:18
% EndTime: 2018-11-14 13:44:18
% DurationCPUTime: 0.10s
% Computational Cost: add. (101->28), mult. (54->35), div. (0->0), fcn. (32->6), ass. (0->18)
t13 = pkin(7) + qJ(2);
t10 = sin(t13);
t20 = pkin(2) * t10;
t12 = qJ(3) + t13;
t7 = sin(t12);
t8 = cos(t12);
t19 = t8 * rSges(4,1) - rSges(4,2) * t7;
t9 = qJ(4) + t12;
t4 = sin(t9);
t5 = cos(t9);
t18 = t5 * rSges(5,1) - rSges(5,2) * t4;
t17 = pkin(3) * t8 + t18;
t16 = -rSges(4,1) * t7 - rSges(4,2) * t8;
t15 = -rSges(5,1) * t4 - rSges(5,2) * t5;
t14 = -pkin(3) * t7 + t15;
t11 = cos(t13);
t6 = pkin(2) * t11;
t1 = [(-m(2) - m(3) - m(4) - m(5)) * g(3), -m(3) * (g(1) * (-rSges(3,1) * t10 - rSges(3,2) * t11) + g(2) * (rSges(3,1) * t11 - rSges(3,2) * t10)) - m(4) * (g(1) * (t16 - t20) + g(2) * (t19 + t6)) - m(5) * (g(1) * (t14 - t20) + g(2) * (t17 + t6)) -m(4) * (g(1) * t16 + g(2) * t19) - m(5) * (g(1) * t14 + g(2) * t17) -m(5) * (g(1) * t15 + g(2) * t18)];
taug  = t1(:);
