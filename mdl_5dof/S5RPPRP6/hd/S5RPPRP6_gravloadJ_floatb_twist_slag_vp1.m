% Calculate Gravitation load on the joints for
% S5RPPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
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
% Datum: 2019-12-31 17:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPRP6_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP6_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP6_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP6_gravloadJ_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP6_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRP6_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:54:58
% EndTime: 2019-12-31 17:54:59
% DurationCPUTime: 0.28s
% Computational Cost: add. (124->56), mult. (156->71), div. (0->0), fcn. (125->6), ass. (0->20)
t15 = sin(qJ(1));
t16 = cos(qJ(1));
t34 = -g(1) * t15 + g(2) * t16;
t31 = rSges(6,3) + qJ(5);
t32 = rSges(6,1) + pkin(4);
t11 = pkin(7) + qJ(4);
t6 = sin(t11);
t7 = cos(t11);
t33 = t31 * t7 - t32 * t6;
t12 = sin(pkin(7));
t29 = pkin(3) * t12;
t26 = t16 * pkin(1) + t15 * qJ(2);
t25 = rSges(4,3) + qJ(3);
t24 = -m(4) - m(5) - m(6);
t22 = rSges(5,1) * t6 + rSges(5,2) * t7;
t20 = rSges(4,1) * t12 + rSges(4,2) * cos(pkin(7));
t14 = -pkin(6) - qJ(3);
t9 = t16 * qJ(2);
t19 = g(1) * (t15 * t14 + t16 * t29 + t9) + g(2) * (t15 * t29 + t26);
t1 = [-m(2) * (g(1) * (-t15 * rSges(2,1) - rSges(2,2) * t16) + g(2) * (rSges(2,1) * t16 - t15 * rSges(2,2))) - m(3) * (g(1) * (rSges(3,3) * t16 + t9 + (rSges(3,2) - pkin(1)) * t15) + g(2) * (-rSges(3,2) * t16 + t15 * rSges(3,3) + t26)) - m(4) * (g(1) * t9 + g(2) * t26 + (g(1) * t20 + g(2) * t25) * t16 + (g(1) * (-pkin(1) - t25) + g(2) * t20) * t15) - m(5) * ((g(1) * t22 + g(2) * (rSges(5,3) - t14)) * t16 + (g(1) * (-rSges(5,3) - pkin(1)) + g(2) * t22) * t15 + t19) - m(6) * ((-g(1) * t33 + g(2) * (rSges(6,2) - t14)) * t16 + (g(1) * (-rSges(6,2) - pkin(1)) - g(2) * t33) * t15 + t19), -(-m(3) + t24) * t34, t24 * (g(1) * t16 + g(2) * t15), (m(5) * t22 - m(6) * t33) * g(3) + t34 * (m(5) * (rSges(5,1) * t7 - rSges(5,2) * t6) + m(6) * (t31 * t6 + t32 * t7)), -m(6) * (g(3) * t6 + t34 * t7)];
taug = t1(:);
