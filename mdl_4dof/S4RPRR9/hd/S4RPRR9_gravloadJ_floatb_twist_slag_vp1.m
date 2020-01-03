% Calculate Gravitation load on the joints for
% S4RPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
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
% Datum: 2019-12-31 16:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RPRR9_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR9_gravloadJ_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR9_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR9_gravloadJ_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR9_gravloadJ_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRR9_gravloadJ_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:56:08
% EndTime: 2019-12-31 16:56:09
% DurationCPUTime: 0.28s
% Computational Cost: add. (78->54), mult. (156->76), div. (0->0), fcn. (137->6), ass. (0->25)
t12 = sin(qJ(1));
t15 = cos(qJ(1));
t37 = -g(1) * t12 + g(2) * t15;
t11 = sin(qJ(3));
t14 = cos(qJ(3));
t25 = rSges(5,3) + pkin(6);
t35 = pkin(3) * t11 - t25 * t14;
t10 = sin(qJ(4));
t13 = cos(qJ(4));
t16 = m(4) * rSges(4,1) + m(5) * (rSges(5,1) * t13 - rSges(5,2) * t10 + pkin(3));
t33 = -m(4) * rSges(4,2) + m(5) * t25;
t32 = -pkin(1) - pkin(5);
t31 = t15 * pkin(1) + t12 * qJ(2);
t24 = rSges(4,2) * t14;
t23 = t11 * t15;
t22 = t12 * t10;
t21 = t12 * t13;
t20 = t13 * t15;
t19 = t15 * pkin(5) + t31;
t7 = t15 * qJ(2);
t4 = t11 * t20 - t22;
t3 = t10 * t23 + t21;
t2 = t10 * t15 + t11 * t21;
t1 = -t11 * t22 + t20;
t5 = [-m(2) * (g(1) * (-t12 * rSges(2,1) - rSges(2,2) * t15) + g(2) * (rSges(2,1) * t15 - t12 * rSges(2,2))) - m(3) * (g(1) * (rSges(3,3) * t15 + t7 + (rSges(3,2) - pkin(1)) * t12) + g(2) * (-rSges(3,2) * t15 + t12 * rSges(3,3) + t31)) - m(4) * (g(1) * (rSges(4,1) * t23 + t15 * t24 + t7) + g(2) * (rSges(4,3) * t15 + t19) + (g(1) * (-rSges(4,3) + t32) + g(2) * (rSges(4,1) * t11 + t24)) * t12) - m(5) * ((rSges(5,1) * t2 + rSges(5,2) * t1 + t35 * t12 + t19) * g(2) + (t4 * rSges(5,1) - t3 * rSges(5,2) + t32 * t12 + t35 * t15 + t7) * g(1)), -(-m(3) - m(4) - m(5)) * t37, (t16 * t11 - t14 * t33) * g(3) + t37 * (t33 * t11 + t16 * t14), -m(5) * (g(1) * (rSges(5,1) * t1 - rSges(5,2) * t2) + g(2) * (rSges(5,1) * t3 + rSges(5,2) * t4) + g(3) * (-rSges(5,1) * t10 - rSges(5,2) * t13) * t14)];
taug = t5(:);
