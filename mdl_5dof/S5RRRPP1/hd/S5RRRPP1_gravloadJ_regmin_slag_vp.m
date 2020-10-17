% Calculate minimal parameter regressor of gravitation load for
% S5RRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% 
% Output:
% taug_reg [5x19]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRPP1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP1_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP1_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:49:57
% EndTime: 2019-12-31 20:49:58
% DurationCPUTime: 0.15s
% Computational Cost: add. (189->41), mult. (159->47), div. (0->0), fcn. (145->8), ass. (0->31)
t19 = qJ(3) + pkin(8);
t13 = sin(t19);
t14 = cos(t19);
t39 = t14 * pkin(4) + t13 * qJ(5);
t20 = qJ(1) + qJ(2);
t15 = sin(t20);
t16 = cos(t20);
t7 = g(1) * t16 + g(2) * t15;
t23 = sin(qJ(1));
t35 = t23 * pkin(1);
t21 = -qJ(4) - pkin(7);
t34 = t16 * t21;
t24 = cos(qJ(3));
t17 = t24 * pkin(3);
t12 = t17 + pkin(2);
t10 = t16 * t12;
t32 = t39 * t16 + t10;
t31 = -t15 * t21 + t10;
t6 = g(1) * t15 - g(2) * t16;
t29 = -t15 * t12 - t34;
t22 = sin(qJ(3));
t27 = -g(3) * t24 + t7 * t22;
t26 = (-g(1) * (-t12 - t39) + g(2) * t21) * t15;
t25 = cos(qJ(1));
t18 = t25 * pkin(1);
t5 = t6 * t24;
t4 = t6 * t22;
t3 = t6 * t14;
t2 = t6 * t13;
t1 = -g(3) * t14 + t7 * t13;
t8 = [0, g(1) * t23 - g(2) * t25, g(1) * t25 + g(2) * t23, 0, t6, t7, 0, 0, 0, 0, 0, t5, -t4, -t7, -g(1) * (t29 - t35) - g(2) * (t18 + t31), t3, -t7, t2, -g(1) * (-t34 - t35) - g(2) * (t18 + t32) + t26; 0, 0, 0, 0, t6, t7, 0, 0, 0, 0, 0, t5, -t4, -t7, -g(1) * t29 - g(2) * t31, t3, -t7, t2, g(1) * t34 - g(2) * t32 + t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, g(3) * t22 + t7 * t24, 0, t27 * pkin(3), t1, 0, -g(3) * t13 - t7 * t14, -g(3) * (t17 + t39) + t7 * (pkin(3) * t22 + pkin(4) * t13 - qJ(5) * t14); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, 0, 0, 0, -t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t8;
