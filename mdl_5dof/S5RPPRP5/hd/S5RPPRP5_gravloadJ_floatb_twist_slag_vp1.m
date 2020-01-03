% Calculate Gravitation load on the joints for
% S5RPPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
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
% Datum: 2019-12-31 17:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPRP5_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP5_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP5_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP5_gravloadJ_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP5_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRP5_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:53:20
% EndTime: 2019-12-31 17:53:22
% DurationCPUTime: 0.41s
% Computational Cost: add. (135->73), mult. (285->103), div. (0->0), fcn. (287->6), ass. (0->27)
t38 = rSges(6,1) + pkin(4);
t19 = cos(pkin(7));
t20 = sin(qJ(4));
t18 = sin(pkin(7));
t36 = cos(qJ(4));
t28 = t18 * t36;
t7 = -t19 * t20 + t28;
t30 = rSges(6,3) + qJ(5);
t22 = cos(qJ(1));
t35 = t18 * t22;
t33 = t19 * t22;
t21 = sin(qJ(1));
t32 = t22 * pkin(1) + t21 * qJ(2);
t31 = qJ(3) * t18;
t29 = -m(4) - m(5) - m(6);
t16 = t22 * qJ(2);
t27 = -t22 * pkin(6) + t16;
t26 = pkin(2) * t33 + t22 * t31 + t32;
t25 = pkin(3) * t33 + t26;
t24 = -pkin(2) * t19 - pkin(1) - t31;
t6 = t18 * t20 + t19 * t36;
t23 = g(1) * (-pkin(3) * t19 + t24);
t5 = t6 * t22;
t4 = t20 * t33 - t22 * t28;
t3 = t6 * t21;
t2 = t7 * t21;
t1 = [-m(2) * (g(1) * (-t21 * rSges(2,1) - t22 * rSges(2,2)) + g(2) * (t22 * rSges(2,1) - t21 * rSges(2,2))) - m(3) * (g(1) * (t22 * rSges(3,3) + t16) + g(2) * (rSges(3,1) * t33 - rSges(3,2) * t35 + t32) + (g(1) * (-rSges(3,1) * t19 + rSges(3,2) * t18 - pkin(1)) + g(2) * rSges(3,3)) * t21) - m(4) * (g(1) * (t22 * rSges(4,2) + t16) + g(2) * (rSges(4,1) * t33 + rSges(4,3) * t35 + t26) + (g(1) * (-rSges(4,1) * t19 - rSges(4,3) * t18 + t24) + g(2) * rSges(4,2)) * t21) - m(5) * (g(1) * (-t3 * rSges(5,1) - t2 * rSges(5,2) - t22 * rSges(5,3) + t27) + g(2) * (t5 * rSges(5,1) - t4 * rSges(5,2) + t25) + (t23 + g(2) * (-rSges(5,3) - pkin(6))) * t21) - m(6) * (g(1) * (-t22 * rSges(6,2) + t30 * t2 - t38 * t3 + t27) + g(2) * (t30 * t4 + t38 * t5 + t25) + (t23 + g(2) * (-rSges(6,2) - pkin(6))) * t21), (-m(3) + t29) * (g(1) * t21 - g(2) * t22), t29 * (-g(3) * t19 + (g(1) * t22 + g(2) * t21) * t18), -m(5) * (g(1) * (-t4 * rSges(5,1) - t5 * rSges(5,2)) + g(2) * (rSges(5,1) * t2 - rSges(5,2) * t3) + g(3) * (-rSges(5,1) * t6 - rSges(5,2) * t7)) - m(6) * (g(1) * (t30 * t5 - t38 * t4) + g(2) * (t2 * t38 + t30 * t3) + g(3) * (t30 * t7 - t38 * t6)), -m(6) * (g(1) * t4 - g(2) * t2 + g(3) * t6)];
taug = t1(:);
