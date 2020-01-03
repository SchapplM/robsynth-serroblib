% Calculate Gravitation load on the joints for
% S5PRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1]';
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
% Datum: 2019-12-31 17:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRPP4_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP4_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPP4_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRRPP4_gravloadJ_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPP4_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPP4_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:40:46
% EndTime: 2019-12-31 17:40:47
% DurationCPUTime: 0.28s
% Computational Cost: add. (162->61), mult. (167->74), div. (0->0), fcn. (135->4), ass. (0->26)
t13 = pkin(7) + qJ(2);
t10 = cos(t13);
t9 = sin(t13);
t35 = g(1) * t10 + g(2) * t9;
t15 = cos(qJ(3));
t14 = sin(qJ(3));
t27 = t14 * rSges(6,2);
t30 = rSges(6,1) + pkin(4);
t36 = t30 * t15 + t27;
t33 = -m(5) - m(6);
t32 = t10 * pkin(2) + t9 * pkin(6);
t12 = t15 * pkin(3);
t29 = t10 * t15;
t28 = t14 * rSges(4,2);
t26 = t14 * rSges(5,3);
t11 = t14 * qJ(4);
t25 = t11 + t12;
t23 = -rSges(6,3) - qJ(5);
t22 = -pkin(3) - t30;
t21 = pkin(3) * t29 + t10 * t11 + t32;
t20 = -pkin(2) - t11;
t19 = t35 * qJ(4) * t15;
t18 = rSges(4,1) * t15 - t28;
t17 = rSges(5,1) * t15 + t26;
t7 = t10 * pkin(6);
t1 = [(-m(2) - m(3) - m(4) + t33) * g(3), -m(3) * (g(1) * (-rSges(3,1) * t9 - rSges(3,2) * t10) + g(2) * (rSges(3,1) * t10 - rSges(3,2) * t9)) - m(4) * (g(1) * (t10 * rSges(4,3) + t7) + g(2) * (rSges(4,1) * t29 - t10 * t28 + t32) + (g(1) * (-pkin(2) - t18) + g(2) * rSges(4,3)) * t9) - m(5) * (g(1) * (t10 * rSges(5,2) + t7) + g(2) * (rSges(5,1) * t29 + t10 * t26 + t21) + (g(1) * (-t17 + t20 - t12) + g(2) * rSges(5,2)) * t9) - m(6) * (g(1) * t7 + g(2) * t21 + (g(1) * t23 + g(2) * t36) * t10 + (g(2) * t23 + (t22 * t15 + t20 - t27) * g(1)) * t9), -m(4) * g(3) * t18 - m(5) * (g(3) * (t17 + t25) + t19) - m(6) * (g(3) * (t25 + t36) + t19) + t35 * ((m(4) * rSges(4,2) - m(5) * rSges(5,3) - m(6) * rSges(6,2)) * t15 + (m(4) * rSges(4,1) - m(5) * (-rSges(5,1) - pkin(3)) - m(6) * t22) * t14), t33 * (-g(3) * t15 + t35 * t14), -m(6) * (-g(1) * t9 + g(2) * t10)];
taug = t1(:);
