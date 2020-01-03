% Calculate Gravitation load on the joints for
% S4PRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
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
% Datum: 2019-12-31 16:29
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4PRRP5_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP5_gravloadJ_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP5_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP5_gravloadJ_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP5_gravloadJ_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRP5_gravloadJ_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:28:51
% EndTime: 2019-12-31 16:28:52
% DurationCPUTime: 0.20s
% Computational Cost: add. (71->36), mult. (157->54), div. (0->0), fcn. (140->6), ass. (0->17)
t24 = rSges(5,1) + pkin(3);
t6 = sin(pkin(6));
t7 = cos(pkin(6));
t23 = g(1) * t7 + g(2) * t6;
t11 = cos(qJ(3));
t9 = sin(qJ(3));
t22 = m(4) * (rSges(4,1) * t11 - rSges(4,2) * t9 + pkin(2)) + m(5) * (-rSges(5,2) * t9 + t24 * t11 + pkin(2)) + m(3) * rSges(3,1);
t12 = cos(qJ(2));
t18 = t12 * t9;
t17 = t11 * t12;
t3 = t6 * t11 - t7 * t18;
t1 = -t7 * t11 - t6 * t18;
t14 = m(3) * rSges(3,2) - m(4) * (rSges(4,3) + pkin(5)) - m(5) * (rSges(5,3) + qJ(4) + pkin(5));
t10 = sin(qJ(2));
t4 = -t7 * t17 - t6 * t9;
t2 = -t6 * t17 + t7 * t9;
t5 = [(-m(2) - m(3) - m(4) - m(5)) * g(3), (t14 * t10 - t22 * t12) * g(3) + t23 * (t22 * t10 + t14 * t12), -m(4) * (g(1) * (rSges(4,1) * t3 + rSges(4,2) * t4) + g(2) * (rSges(4,1) * t1 + rSges(4,2) * t2)) - m(5) * (g(1) * (t4 * rSges(5,2) + t24 * t3) + g(2) * (t2 * rSges(5,2) + t24 * t1)) + (-m(4) * (-rSges(4,1) * t9 - rSges(4,2) * t11) - m(5) * (-rSges(5,2) * t11 - t24 * t9)) * g(3) * t10, -m(5) * (-g(3) * t12 + t23 * t10)];
taug = t5(:);
