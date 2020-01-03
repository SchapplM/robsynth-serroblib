% Calculate Gravitation load on the joints for
% S4PRRR4
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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4PRRR4_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR4_gravloadJ_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR4_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR4_gravloadJ_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR4_gravloadJ_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRR4_gravloadJ_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:32:30
% EndTime: 2019-12-31 16:32:31
% DurationCPUTime: 0.15s
% Computational Cost: add. (104->32), mult. (94->43), div. (0->0), fcn. (69->6), ass. (0->18)
t11 = cos(qJ(3));
t9 = qJ(3) + qJ(4);
t5 = sin(t9);
t6 = cos(t9);
t18 = t6 * rSges(5,1) - rSges(5,2) * t5;
t25 = t11 * pkin(3) + t18;
t8 = pkin(7) + qJ(2);
t3 = sin(t8);
t4 = cos(t8);
t24 = g(1) * t4 + g(2) * t3;
t20 = rSges(4,3) + pkin(5);
t19 = rSges(5,3) + pkin(6) + pkin(5);
t17 = -rSges(5,1) * t5 - rSges(5,2) * t6;
t10 = sin(qJ(3));
t16 = rSges(4,1) * t11 - rSges(4,2) * t10;
t15 = pkin(2) + t25;
t14 = pkin(2) + t16;
t1 = [(-m(2) - m(3) - m(4) - m(5)) * g(3), -m(3) * (g(1) * (-rSges(3,1) * t3 - rSges(3,2) * t4) + g(2) * (rSges(3,1) * t4 - rSges(3,2) * t3)) - m(4) * ((g(1) * t20 + g(2) * t14) * t4 + (-g(1) * t14 + g(2) * t20) * t3) - m(5) * ((g(1) * t19 + g(2) * t15) * t4 + (-g(1) * t15 + g(2) * t19) * t3), (-m(4) * t16 - m(5) * t25) * g(3) + t24 * (-m(4) * (-rSges(4,1) * t10 - rSges(4,2) * t11) - m(5) * (-pkin(3) * t10 + t17)), -m(5) * (g(3) * t18 + t24 * t17)];
taug = t1(:);
