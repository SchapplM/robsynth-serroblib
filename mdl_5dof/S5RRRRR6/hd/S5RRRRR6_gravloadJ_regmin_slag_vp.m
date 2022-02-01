% Calculate minimal parameter regressor of gravitation load for
% S5RRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% 
% Output:
% taug_reg [5x27]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 12:09
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRRR6_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR6_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR6_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR6_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 12:08:22
% EndTime: 2022-01-20 12:08:23
% DurationCPUTime: 0.10s
% Computational Cost: add. (170->19), mult. (120->26), div. (0->0), fcn. (120->10), ass. (0->26)
t20 = qJ(3) + qJ(4);
t21 = qJ(1) + qJ(2);
t16 = sin(t21);
t18 = cos(t21);
t12 = g(1) * t18 + g(2) * t16;
t11 = g(1) * t16 - g(2) * t18;
t25 = cos(qJ(1));
t24 = cos(qJ(3));
t23 = sin(qJ(1));
t22 = sin(qJ(3));
t19 = qJ(5) + t20;
t17 = cos(t20);
t15 = sin(t20);
t14 = cos(t19);
t13 = sin(t19);
t10 = t11 * t24;
t9 = t11 * t22;
t8 = t11 * t17;
t7 = t11 * t15;
t6 = t11 * t14;
t5 = t11 * t13;
t4 = g(3) * t15 + t12 * t17;
t3 = -g(3) * t17 + t12 * t15;
t2 = g(3) * t13 + t12 * t14;
t1 = -g(3) * t14 + t12 * t13;
t26 = [0, g(1) * t23 - g(2) * t25, g(1) * t25 + g(2) * t23, 0, t11, t12, 0, 0, 0, 0, 0, t10, -t9, 0, 0, 0, 0, 0, t8, -t7, 0, 0, 0, 0, 0, t6, -t5; 0, 0, 0, 0, t11, t12, 0, 0, 0, 0, 0, t10, -t9, 0, 0, 0, 0, 0, t8, -t7, 0, 0, 0, 0, 0, t6, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t24 + t12 * t22, g(3) * t22 + t12 * t24, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg = t26;
