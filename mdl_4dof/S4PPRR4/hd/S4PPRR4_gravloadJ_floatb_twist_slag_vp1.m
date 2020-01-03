% Calculate Gravitation load on the joints for
% S4PPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1,theta2]';
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
% Datum: 2019-12-31 16:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4PPRR4_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR4_gravloadJ_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR4_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PPRR4_gravloadJ_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR4_gravloadJ_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PPRR4_gravloadJ_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:18:26
% EndTime: 2019-12-31 16:18:27
% DurationCPUTime: 0.14s
% Computational Cost: add. (65->23), mult. (86->37), div. (0->0), fcn. (71->6), ass. (0->15)
t7 = sin(qJ(4));
t8 = cos(qJ(4));
t18 = m(4) * rSges(4,1) + m(5) * (rSges(5,1) * t8 - rSges(5,2) * t7 + pkin(3));
t5 = sin(pkin(6));
t17 = t5 * t7;
t16 = t5 * t8;
t6 = cos(pkin(6));
t15 = t6 * t7;
t14 = t6 * t8;
t12 = -m(3) - m(4) - m(5);
t10 = m(4) * rSges(4,2) - m(5) * (rSges(5,3) + pkin(5));
t4 = pkin(7) + qJ(3);
t3 = cos(t4);
t2 = sin(t4);
t1 = [(-m(2) + t12) * g(3), t12 * (g(1) * t5 - g(2) * t6), (t10 * t2 - t18 * t3) * g(3) + (g(1) * t6 + g(2) * t5) * (t10 * t3 + t18 * t2), -m(5) * (g(1) * ((-t3 * t15 + t16) * rSges(5,1) + (-t3 * t14 - t17) * rSges(5,2)) + g(2) * ((-t3 * t17 - t14) * rSges(5,1) + (-t3 * t16 + t15) * rSges(5,2)) + g(3) * (-rSges(5,1) * t7 - rSges(5,2) * t8) * t2)];
taug = t1(:);
