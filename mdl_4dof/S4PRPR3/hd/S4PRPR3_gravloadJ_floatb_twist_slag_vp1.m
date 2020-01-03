% Calculate Gravitation load on the joints for
% S4PRPR3
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
% Datum: 2019-12-31 16:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4PRPR3_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR3_gravloadJ_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR3_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR3_gravloadJ_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR3_gravloadJ_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRPR3_gravloadJ_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:20:47
% EndTime: 2019-12-31 16:20:48
% DurationCPUTime: 0.11s
% Computational Cost: add. (81->28), mult. (67->37), div. (0->0), fcn. (48->6), ass. (0->14)
t18 = -m(4) - m(5);
t17 = rSges(5,3) + pkin(5) + qJ(3);
t16 = rSges(4,3) + qJ(3);
t7 = pkin(7) + qJ(4);
t3 = sin(t7);
t5 = cos(t7);
t15 = rSges(5,1) * t5 - rSges(5,2) * t3;
t10 = cos(pkin(7));
t13 = pkin(3) * t10 + pkin(2) + t15;
t12 = rSges(4,1) * t10 - rSges(4,2) * sin(pkin(7)) + pkin(2);
t8 = pkin(6) + qJ(2);
t6 = cos(t8);
t4 = sin(t8);
t1 = [(-m(2) - m(3) + t18) * g(3), -m(3) * (g(1) * (-rSges(3,1) * t4 - rSges(3,2) * t6) + g(2) * (rSges(3,1) * t6 - rSges(3,2) * t4)) - m(4) * ((g(1) * t16 + g(2) * t12) * t6 + (-g(1) * t12 + g(2) * t16) * t4) - m(5) * ((g(1) * t17 + g(2) * t13) * t6 + (-g(1) * t13 + g(2) * t17) * t4), t18 * (g(1) * t4 - g(2) * t6), -m(5) * (g(3) * t15 + (g(1) * t6 + g(2) * t4) * (-rSges(5,1) * t3 - rSges(5,2) * t5))];
taug = t1(:);
