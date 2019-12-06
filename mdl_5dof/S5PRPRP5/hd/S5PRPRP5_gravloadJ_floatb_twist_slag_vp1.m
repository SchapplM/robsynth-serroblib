% Calculate Gravitation load on the joints for
% S5PRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-05 15:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRPRP5_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP5_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP5_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP5_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP5_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRP5_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:37:35
% EndTime: 2019-12-05 15:37:36
% DurationCPUTime: 0.32s
% Computational Cost: add. (167->54), mult. (243->79), div. (0->0), fcn. (224->8), ass. (0->29)
t33 = rSges(6,1) + pkin(4);
t14 = sin(pkin(7));
t16 = cos(pkin(7));
t38 = g(1) * t16 + g(2) * t14;
t27 = rSges(6,3) + qJ(5);
t12 = pkin(8) + qJ(4);
t10 = sin(t12);
t11 = cos(t12);
t37 = t27 * t10 + t33 * t11;
t18 = sin(qJ(2));
t34 = g(3) * t18;
t19 = cos(qJ(2));
t32 = t14 * t19;
t31 = t16 * t19;
t17 = -pkin(6) - qJ(3);
t30 = rSges(6,2) - t17;
t29 = rSges(5,3) - t17;
t28 = rSges(4,3) + qJ(3);
t26 = -m(4) - m(5) - m(6);
t24 = rSges(5,1) * t11 - rSges(5,2) * t10;
t15 = cos(pkin(8));
t23 = rSges(4,1) * t15 - rSges(4,2) * sin(pkin(8)) + pkin(2);
t9 = t15 * pkin(3) + pkin(2);
t6 = t19 * t9;
t4 = t14 * t10 + t11 * t31;
t3 = t10 * t31 - t14 * t11;
t2 = -t16 * t10 + t11 * t32;
t1 = t10 * t32 + t16 * t11;
t5 = [(-m(2) - m(3) + t26) * g(3), -m(3) * (g(3) * (t19 * rSges(3,1) - t18 * rSges(3,2)) + t38 * (-rSges(3,1) * t18 - rSges(3,2) * t19)) - m(4) * (g(3) * (t28 * t18 + t23 * t19) + t38 * (-t23 * t18 + t28 * t19)) - m(5) * (g(3) * (t29 * t18 + t24 * t19 + t6) + t38 * (t29 * t19 + (-t24 - t9) * t18)) - m(6) * (g(3) * t6 + (g(3) * t37 + t38 * t30) * t19 + (g(3) * t30 + t38 * (-t37 - t9)) * t18), t26 * (-g(3) * t19 + t38 * t18), -m(5) * (g(1) * (-rSges(5,1) * t3 - rSges(5,2) * t4) + g(2) * (-rSges(5,1) * t1 - rSges(5,2) * t2)) - m(6) * (g(1) * (t27 * t4 - t33 * t3) + g(2) * (-t33 * t1 + t27 * t2)) + (-m(5) * (-rSges(5,1) * t10 - rSges(5,2) * t11) - m(6) * (-t33 * t10 + t27 * t11)) * t34, -m(6) * (g(1) * t3 + g(2) * t1 + t10 * t34)];
taug = t5(:);
