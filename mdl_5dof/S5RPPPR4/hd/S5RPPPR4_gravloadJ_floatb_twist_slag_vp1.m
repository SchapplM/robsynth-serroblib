% Calculate Gravitation load on the joints for
% S5RPPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPPR4_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR4_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR4_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR4_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR4_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPPR4_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:45:06
% EndTime: 2019-12-31 17:45:07
% DurationCPUTime: 0.20s
% Computational Cost: add. (130->47), mult. (108->59), div. (0->0), fcn. (80->8), ass. (0->23)
t12 = qJ(1) + pkin(7);
t7 = sin(t12);
t9 = cos(t12);
t32 = g(1) * t7 - g(2) * t9;
t29 = -m(5) - m(6);
t16 = sin(qJ(1));
t28 = pkin(1) * t16;
t27 = rSges(6,3) + pkin(6) + qJ(4);
t26 = rSges(5,3) + qJ(4);
t25 = -m(4) + t29;
t17 = cos(qJ(1));
t10 = t17 * pkin(1);
t24 = t9 * pkin(2) + t7 * qJ(3) + t10;
t23 = t9 * qJ(3) - t28;
t11 = pkin(8) + qJ(5);
t6 = sin(t11);
t8 = cos(t11);
t21 = -rSges(6,1) * t6 - rSges(6,2) * t8;
t13 = sin(pkin(8));
t20 = rSges(5,1) * t13 + rSges(5,2) * cos(pkin(8));
t19 = pkin(4) * t13 - t21;
t18 = g(1) * t23 + g(2) * t24;
t1 = [-m(2) * (g(1) * (-t16 * rSges(2,1) - rSges(2,2) * t17) + g(2) * (rSges(2,1) * t17 - t16 * rSges(2,2))) - m(3) * (g(1) * (-rSges(3,1) * t7 - rSges(3,2) * t9 - t28) + g(2) * (rSges(3,1) * t9 - rSges(3,2) * t7 + t10)) - m(4) * (g(1) * (rSges(4,3) * t9 + (rSges(4,2) - pkin(2)) * t7 + t23) + g(2) * (-rSges(4,2) * t9 + rSges(4,3) * t7 + t24)) - m(5) * ((g(1) * t20 + g(2) * t26) * t9 + (g(1) * (-pkin(2) - t26) + g(2) * t20) * t7 + t18) - m(6) * ((g(1) * t19 + g(2) * t27) * t9 + (g(1) * (-pkin(2) - t27) + g(2) * t19) * t7 + t18), (-m(3) + t25) * g(3), t25 * t32, t29 * (g(1) * t9 + g(2) * t7), -m(6) * (g(3) * t21 + t32 * (rSges(6,1) * t8 - rSges(6,2) * t6))];
taug = t1(:);
