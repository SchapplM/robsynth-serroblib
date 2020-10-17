% Calculate inertial parameters regressor of gravitation load for
% S5PRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRRPP4_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP4_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPP4_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRRPP4_gravloadJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:41:12
% EndTime: 2019-12-31 17:41:13
% DurationCPUTime: 0.15s
% Computational Cost: add. (141->41), mult. (147->45), div. (0->0), fcn. (135->4), ass. (0->26)
t20 = sin(qJ(3));
t19 = pkin(7) + qJ(2);
t14 = sin(t19);
t15 = cos(t19);
t6 = g(1) * t15 + g(2) * t14;
t34 = t6 * t20;
t16 = t20 * qJ(4);
t21 = cos(qJ(3));
t25 = t21 * pkin(3) + t16;
t32 = pkin(3) * t20;
t31 = g(1) * t14;
t28 = t21 * pkin(4);
t27 = t15 * t21;
t26 = t15 * pkin(2) + t14 * pkin(6);
t24 = qJ(4) * t21;
t23 = pkin(3) * t27 + t15 * t16 + t26;
t5 = -g(2) * t15 + t31;
t22 = -pkin(2) - t25;
t12 = t15 * pkin(6);
t9 = t15 * t24;
t7 = t14 * t24;
t4 = t5 * t21;
t3 = t5 * t20;
t2 = g(3) * t20 + t6 * t21;
t1 = -g(3) * t21 + t34;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, 0, 0, 0, t4, -t3, -t6, -g(1) * (-t14 * pkin(2) + t12) - g(2) * t26, 0, 0, 0, 0, 0, 0, t4, -t6, t3, -g(1) * t12 - g(2) * t23 - t22 * t31, 0, 0, 0, 0, 0, 0, t4, t3, t6, -g(1) * (-t15 * qJ(5) + t12) - g(2) * (pkin(4) * t27 + t23) + (-g(1) * (t22 - t28) + g(2) * qJ(5)) * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t2, -g(1) * (-t15 * t32 + t9) - g(2) * (-t14 * t32 + t7) - g(3) * t25, 0, 0, 0, 0, 0, 0, t1, -t2, 0, -g(1) * t9 - g(2) * t7 - g(3) * (t25 + t28) + (pkin(3) + pkin(4)) * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5;];
taug_reg = t8;
