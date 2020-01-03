% Calculate Gravitation load on the joints for
% S4RPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
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
% Datum: 2019-12-31 16:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RPRP5_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP5_gravloadJ_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP5_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP5_gravloadJ_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP5_gravloadJ_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRP5_gravloadJ_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:44:46
% EndTime: 2019-12-31 16:44:47
% DurationCPUTime: 0.24s
% Computational Cost: add. (103->42), mult. (123->54), div. (0->0), fcn. (99->6), ass. (0->20)
t22 = rSges(5,3) + qJ(4);
t24 = rSges(5,1) + pkin(3);
t6 = pkin(6) + qJ(3);
t4 = sin(t6);
t5 = cos(t6);
t13 = t22 * t4 + t24 * t5;
t11 = cos(qJ(1));
t25 = g(2) * t11;
t10 = sin(qJ(1));
t23 = g(1) * t11 + g(2) * t10;
t8 = cos(pkin(6));
t3 = pkin(2) * t8 + pkin(1);
t21 = t3 * t25;
t9 = -pkin(5) - qJ(2);
t18 = rSges(5,2) - t9;
t17 = rSges(4,3) - t9;
t16 = rSges(3,3) + qJ(2);
t15 = rSges(4,1) * t5 - rSges(4,2) * t4;
t14 = rSges(3,1) * t8 - rSges(3,2) * sin(pkin(6)) + pkin(1);
t1 = [-m(2) * (g(1) * (-t10 * rSges(2,1) - rSges(2,2) * t11) + g(2) * (rSges(2,1) * t11 - t10 * rSges(2,2))) - m(3) * ((g(1) * t16 + g(2) * t14) * t11 + (-g(1) * t14 + g(2) * t16) * t10) - m(4) * (t21 + (g(1) * t17 + g(2) * t15) * t11 + (g(1) * (-t15 - t3) + g(2) * t17) * t10) - m(5) * (t21 + (g(1) * t18 + g(2) * t13) * t11 + (g(1) * (-t13 - t3) + g(2) * t18) * t10), (-m(3) - m(4) - m(5)) * (g(1) * t10 - t25), (-m(4) * t15 - m(5) * t13) * g(3) + t23 * (-m(4) * (-rSges(4,1) * t4 - rSges(4,2) * t5) - m(5) * (t22 * t5 - t24 * t4)), -m(5) * (-g(3) * t5 + t23 * t4)];
taug = t1(:);
