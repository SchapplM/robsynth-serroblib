% Calculate Gravitation load on the joints for
% S5RPRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
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
% Datum: 2019-12-31 18:57
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRP12_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP12_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP12_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP12_gravloadJ_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP12_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP12_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:56:16
% EndTime: 2019-12-31 18:56:18
% DurationCPUTime: 0.54s
% Computational Cost: add. (144->81), mult. (282->109), div. (0->0), fcn. (258->6), ass. (0->30)
t14 = sin(qJ(1));
t17 = cos(qJ(1));
t47 = -g(1) * t14 + g(2) * t17;
t13 = sin(qJ(3));
t16 = cos(qJ(3));
t33 = rSges(5,3) + pkin(7);
t45 = t13 * pkin(3) - t33 * t16;
t44 = rSges(6,1) + pkin(4);
t25 = rSges(6,3) + qJ(5) + pkin(7);
t42 = t13 * rSges(4,1) + t16 * rSges(4,2);
t12 = sin(qJ(4));
t15 = cos(qJ(4));
t6 = t15 * pkin(4) + pkin(3);
t18 = m(5) * (rSges(5,1) * t15 - rSges(5,2) * t12 + pkin(3)) + m(6) * (rSges(6,1) * t15 - rSges(6,2) * t12 + t6) + m(4) * rSges(4,1);
t41 = -m(4) * rSges(4,2) + m(5) * t33 + m(6) * t25;
t40 = -pkin(1) - pkin(6);
t37 = pkin(4) * t12;
t32 = t17 * pkin(1) + t14 * qJ(2);
t30 = t14 * t12;
t29 = t14 * t15;
t27 = t17 * t12;
t26 = t17 * t15;
t24 = t17 * pkin(6) + t32;
t3 = t13 * t27 + t29;
t1 = -t13 * t30 + t26;
t20 = t13 * t6 - t25 * t16;
t8 = t17 * qJ(2);
t4 = t13 * t26 - t30;
t2 = t13 * t29 + t27;
t5 = [-m(2) * (g(1) * (-t14 * rSges(2,1) - t17 * rSges(2,2)) + g(2) * (t17 * rSges(2,1) - t14 * rSges(2,2))) - m(3) * (g(1) * (t17 * rSges(3,3) + t8 + (rSges(3,2) - pkin(1)) * t14) + g(2) * (-t17 * rSges(3,2) + t14 * rSges(3,3) + t32)) - m(4) * (g(1) * (t42 * t17 + t8) + g(2) * (t17 * rSges(4,3) + t24) + (g(1) * (-rSges(4,3) + t40) + g(2) * t42) * t14) - m(5) * ((t2 * rSges(5,1) + t1 * rSges(5,2) + t45 * t14 + t24) * g(2) + (t4 * rSges(5,1) - t3 * rSges(5,2) + t40 * t14 + t45 * t17 + t8) * g(1)) - m(6) * (g(1) * (t4 * rSges(6,1) - t3 * rSges(6,2) + t8) + g(2) * (t2 * rSges(6,1) + t1 * rSges(6,2) + t24) + (g(1) * t20 + g(2) * t37) * t17 + (g(1) * (-t37 + t40) + g(2) * t20) * t14), -(-m(3) - m(4) - m(5) - m(6)) * t47, (t18 * t13 - t16 * t41) * g(3) + t47 * (t41 * t13 + t18 * t16), -m(5) * (g(1) * (t1 * rSges(5,1) - t2 * rSges(5,2)) + g(2) * (t3 * rSges(5,1) + t4 * rSges(5,2))) - m(6) * (g(1) * (-t2 * rSges(6,2) + t44 * t1) + g(2) * (t4 * rSges(6,2) + t44 * t3)) + (-m(5) * (-rSges(5,1) * t12 - rSges(5,2) * t15) - m(6) * (-rSges(6,1) * t12 - rSges(6,2) * t15 - t37)) * g(3) * t16, -m(6) * (g(3) * t13 + t47 * t16)];
taug = t5(:);
