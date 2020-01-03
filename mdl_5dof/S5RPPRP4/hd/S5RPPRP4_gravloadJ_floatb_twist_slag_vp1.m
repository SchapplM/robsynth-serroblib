% Calculate Gravitation load on the joints for
% S5RPPRP4
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
% Datum: 2019-12-31 17:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPRP4_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP4_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP4_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP4_gravloadJ_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP4_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRP4_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:52:03
% EndTime: 2019-12-31 17:52:04
% DurationCPUTime: 0.24s
% Computational Cost: add. (118->50), mult. (202->70), div. (0->0), fcn. (207->6), ass. (0->22)
t13 = sin(qJ(4));
t14 = cos(qJ(4));
t32 = rSges(6,1) + pkin(4);
t33 = t13 * rSges(6,2) - t32 * t14;
t31 = rSges(5,3) + pkin(6);
t28 = sin(qJ(1));
t29 = cos(qJ(1));
t30 = t29 * pkin(1) + t28 * qJ(2);
t26 = rSges(6,3) + qJ(5) + pkin(6);
t25 = cos(pkin(7));
t24 = sin(pkin(7));
t23 = m(4) + m(5) + m(6);
t22 = t29 * pkin(2) + t30;
t21 = -t28 * pkin(1) + t29 * qJ(2);
t20 = -rSges(5,1) * t14 + t13 * rSges(5,2);
t19 = pkin(3) - t20;
t18 = pkin(3) - t33;
t17 = -t28 * pkin(2) + t21;
t16 = g(1) * t17 + g(2) * t22;
t2 = t29 * t24 - t28 * t25;
t1 = -t28 * t24 - t29 * t25;
t3 = [-m(2) * (g(1) * (-t28 * rSges(2,1) - t29 * rSges(2,2)) + g(2) * (t29 * rSges(2,1) - t28 * rSges(2,2))) - m(3) * (g(1) * (-t28 * rSges(3,1) + t29 * rSges(3,3) + t21) + g(2) * (t29 * rSges(3,1) + t28 * rSges(3,3) + t30)) - m(4) * (g(1) * (t2 * rSges(4,1) - t1 * rSges(4,2) + t17) + g(2) * (-rSges(4,1) * t1 - rSges(4,2) * t2 + t22)) - m(5) * ((g(1) * t19 + g(2) * t31) * t2 + (g(1) * t31 - g(2) * t19) * t1 + t16) - m(6) * ((g(1) * t18 + g(2) * t26) * t2 + (g(1) * t26 - g(2) * t18) * t1 + t16), (-m(3) - t23) * (g(1) * t28 - g(2) * t29), t23 * g(3), (-m(5) * t20 - m(6) * t33) * g(3) + (g(1) * t1 + g(2) * t2) * (-m(5) * (rSges(5,1) * t13 + rSges(5,2) * t14) - m(6) * (rSges(6,2) * t14 + t32 * t13)), -m(6) * (g(1) * t2 - g(2) * t1)];
taug = t3(:);
