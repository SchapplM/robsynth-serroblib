% Calculate inertial parameters regressor of gravitation load for
% S5RPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:23
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRPR4_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR4_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR4_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR4_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:23:10
% EndTime: 2022-01-23 09:23:11
% DurationCPUTime: 0.17s
% Computational Cost: add. (175->47), mult. (135->54), div. (0->0), fcn. (121->10), ass. (0->28)
t24 = sin(qJ(1));
t30 = t24 * pkin(1);
t20 = qJ(3) + pkin(9);
t14 = cos(t20);
t25 = cos(qJ(3));
t17 = t25 * pkin(3);
t29 = pkin(4) * t14 + t17;
t22 = -qJ(4) - pkin(6);
t21 = qJ(1) + pkin(8);
t13 = sin(t21);
t15 = cos(t21);
t4 = g(1) * t15 + g(2) * t13;
t3 = g(1) * t13 - g(2) * t15;
t26 = cos(qJ(1));
t28 = g(1) * t24 - g(2) * t26;
t23 = sin(qJ(3));
t27 = -g(3) * t25 + t4 * t23;
t19 = pkin(7) - t22;
t18 = t26 * pkin(1);
t16 = qJ(5) + t20;
t12 = sin(t20);
t11 = t17 + pkin(2);
t10 = cos(t16);
t9 = sin(t16);
t5 = pkin(2) + t29;
t2 = g(3) * t9 + t4 * t10;
t1 = -g(3) * t10 + t4 * t9;
t6 = [0, 0, 0, 0, 0, 0, t28, g(1) * t26 + g(2) * t24, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, t28 * pkin(1), 0, 0, 0, 0, 0, 0, t3 * t25, -t3 * t23, -t4, -g(1) * (-t13 * pkin(2) + t15 * pkin(6) - t30) - g(2) * (t15 * pkin(2) + t13 * pkin(6) + t18), 0, 0, 0, 0, 0, 0, t3 * t14, -t3 * t12, -t4, -g(1) * (-t13 * t11 - t15 * t22 - t30) - g(2) * (t15 * t11 - t13 * t22 + t18), 0, 0, 0, 0, 0, 0, t3 * t10, -t3 * t9, -t4, -g(1) * (-t13 * t5 + t19 * t15 - t30) - g(2) * (t13 * t19 + t15 * t5 + t18); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, g(3) * t23 + t4 * t25, 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t14 + t4 * t12, g(3) * t12 + t4 * t14, 0, t27 * pkin(3), 0, 0, 0, 0, 0, 0, t1, t2, 0, -g(3) * t29 - t4 * (-t23 * pkin(3) - pkin(4) * t12); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg = t6;
