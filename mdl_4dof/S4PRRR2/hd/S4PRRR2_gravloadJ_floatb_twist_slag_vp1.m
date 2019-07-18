% Calculate Gravitation load on the joints for
% S4PRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4]';
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
% Datum: 2019-07-18 13:27
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4PRRR2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(2,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S4PRRR2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR2_gravloadJ_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRR2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:27:22
% EndTime: 2019-07-18 13:27:22
% DurationCPUTime: 0.12s
% Computational Cost: add. (69->27), mult. (54->35), div. (0->0), fcn. (32->6), ass. (0->17)
t9 = sin(qJ(2));
t18 = pkin(1) * t9;
t10 = cos(qJ(2));
t17 = pkin(1) * t10;
t8 = qJ(2) + qJ(3);
t5 = sin(t8);
t6 = cos(t8);
t16 = -t6 * rSges(4,1) + t5 * rSges(4,2);
t7 = qJ(4) + t8;
t3 = sin(t7);
t4 = cos(t7);
t15 = -t4 * rSges(5,1) + t3 * rSges(5,2);
t14 = -rSges(4,1) * t5 - rSges(4,2) * t6;
t13 = -rSges(5,1) * t3 - rSges(5,2) * t4;
t12 = -pkin(2) * t6 + t15;
t11 = -pkin(2) * t5 + t13;
t1 = [(m(2) + m(3) + m(4) + m(5)) * g(2), -m(3) * (g(1) * (-rSges(3,1) * t10 + t9 * rSges(3,2)) + g(3) * (-t9 * rSges(3,1) - rSges(3,2) * t10)) - m(4) * (g(1) * (t16 - t17) + g(3) * (t14 - t18)) - m(5) * (g(1) * (t12 - t17) + g(3) * (t11 - t18)), -m(4) * (g(1) * t16 + g(3) * t14) - m(5) * (g(1) * t12 + g(3) * t11), -m(5) * (g(1) * t15 + g(3) * t13)];
taug  = t1(:);
