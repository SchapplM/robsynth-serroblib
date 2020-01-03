% Calculate Gravitation load on the joints for
% S5RPRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3]';
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
% Datum: 2019-12-31 18:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRPP5_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(6,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP5_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP5_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RPRPP5_gravloadJ_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP5_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPP5_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:16:08
% EndTime: 2019-12-31 18:16:09
% DurationCPUTime: 0.30s
% Computational Cost: add. (105->73), mult. (191->94), div. (0->0), fcn. (157->4), ass. (0->28)
t38 = rSges(6,1) + pkin(4);
t14 = cos(qJ(3));
t8 = t14 * qJ(4);
t37 = -rSges(5,3) * t14 - t8;
t36 = -rSges(6,2) * t14 - t8;
t35 = -m(5) - m(6);
t34 = -pkin(1) - pkin(6);
t13 = sin(qJ(1));
t33 = g(1) * t13;
t15 = cos(qJ(1));
t32 = g(2) * t15;
t31 = -rSges(5,1) - pkin(3);
t30 = t15 * pkin(1) + t13 * qJ(2);
t29 = rSges(4,2) * t14;
t12 = sin(qJ(3));
t26 = t12 * t13;
t25 = t12 * t15;
t24 = -rSges(6,3) - qJ(5);
t23 = -pkin(3) - t38;
t22 = t15 * pkin(6) + t30;
t21 = g(1) * (t13 * t14 * pkin(3) + qJ(4) * t26);
t20 = -t32 + t33;
t19 = rSges(4,1) * t12 + t29;
t18 = rSges(5,1) * t12 + t37;
t9 = t15 * qJ(2);
t17 = g(1) * (pkin(3) * t25 + t9) + g(2) * (pkin(3) * t26 + t22);
t16 = t38 * t12 + t36;
t1 = [-m(2) * (g(1) * (-t13 * rSges(2,1) - rSges(2,2) * t15) + g(2) * (rSges(2,1) * t15 - t13 * rSges(2,2))) - m(3) * (g(1) * (rSges(3,3) * t15 + t9 + (rSges(3,2) - pkin(1)) * t13) + g(2) * (-rSges(3,2) * t15 + t13 * rSges(3,3) + t30)) - m(4) * (g(1) * (rSges(4,1) * t25 + t15 * t29 + t9) + g(2) * (rSges(4,3) * t15 + t22) + (g(1) * (-rSges(4,3) + t34) + g(2) * t19) * t13) - m(5) * ((g(2) * rSges(5,2) + g(1) * t18) * t15 + (g(1) * (-rSges(5,2) + t34) + g(2) * t18) * t13 + t17) - m(6) * ((g(1) * t16 + g(2) * t24) * t15 + (g(1) * (-t24 + t34) + g(2) * t16) * t13 + t17), (-m(3) - m(4) + t35) * t20, m(4) * g(3) * t19 - m(5) * (t21 + g(3) * (t31 * t12 - t37)) - m(6) * (t21 + g(3) * (t23 * t12 - t36)) + (-m(4) * (rSges(4,1) * t14 - rSges(4,2) * t12) - m(5) * (rSges(5,1) * t14 + rSges(5,3) * t12) - m(6) * (rSges(6,2) * t12 + t38 * t14)) * t33 + ((-m(4) * rSges(4,2) - m(5) * (-rSges(5,3) - qJ(4)) - m(6) * (-rSges(6,2) - qJ(4))) * t12 + (m(4) * rSges(4,1) - m(5) * t31 - m(6) * t23) * t14) * t32, t35 * (g(3) * t12 - t20 * t14), -m(6) * (-g(1) * t15 - g(2) * t13)];
taug = t1(:);
