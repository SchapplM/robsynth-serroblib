% Calculate Gravitation load on the joints for
% S4PRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:24
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4PRPR6_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR6_gravloadJ_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR6_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR6_gravloadJ_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR6_gravloadJ_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRPR6_gravloadJ_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:24:20
% EndTime: 2019-12-31 16:24:20
% DurationCPUTime: 0.16s
% Computational Cost: add. (75->31), mult. (127->44), div. (0->0), fcn. (109->8), ass. (0->15)
t7 = sin(pkin(6));
t9 = cos(pkin(6));
t24 = g(1) * t9 + g(2) * t7;
t5 = pkin(7) + qJ(4);
t3 = sin(t5);
t4 = cos(t5);
t8 = cos(pkin(7));
t23 = m(4) * (rSges(4,1) * t8 - rSges(4,2) * sin(pkin(7)) + pkin(2)) + m(5) * (rSges(5,1) * t4 - rSges(5,2) * t3 + pkin(3) * t8 + pkin(2)) + m(3) * rSges(3,1);
t20 = -m(4) - m(5);
t12 = cos(qJ(2));
t18 = t12 * t7;
t17 = t12 * t9;
t14 = m(3) * rSges(3,2) - m(4) * (rSges(4,3) + qJ(3)) - m(5) * (rSges(5,3) + pkin(5) + qJ(3));
t11 = sin(qJ(2));
t1 = [(-m(2) - m(3) + t20) * g(3), (t14 * t11 - t23 * t12) * g(3) + t24 * (t23 * t11 + t14 * t12), t20 * (-g(3) * t12 + t24 * t11), -m(5) * (g(1) * ((-t3 * t17 + t7 * t4) * rSges(5,1) + (-t4 * t17 - t7 * t3) * rSges(5,2)) + g(2) * ((-t3 * t18 - t9 * t4) * rSges(5,1) + (-t4 * t18 + t9 * t3) * rSges(5,2)) + g(3) * (-rSges(5,1) * t3 - rSges(5,2) * t4) * t11)];
taug = t1(:);
