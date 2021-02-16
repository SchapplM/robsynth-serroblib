% Calculate minimal parameter regressor of gravitation load for
% S5RRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% taug_reg [5x25]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 20:30
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPRP6_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP6_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP6_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP6_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 20:29:58
% EndTime: 2021-01-15 20:29:59
% DurationCPUTime: 0.20s
% Computational Cost: add. (165->45), mult. (238->61), div. (0->0), fcn. (241->8), ass. (0->39)
t30 = cos(qJ(4));
t18 = t30 * pkin(4) + pkin(3);
t24 = qJ(2) + pkin(8);
t20 = sin(t24);
t21 = cos(t24);
t25 = -qJ(5) - pkin(7);
t49 = t21 * t18 - t20 * t25;
t29 = sin(qJ(1));
t32 = cos(qJ(1));
t15 = g(1) * t32 + g(2) * t29;
t37 = t32 * t30;
t27 = sin(qJ(4));
t40 = t29 * t27;
t10 = t21 * t40 + t37;
t38 = t32 * t27;
t39 = t29 * t30;
t12 = -t21 * t38 + t39;
t44 = g(3) * t20;
t1 = -g(1) * t12 + g(2) * t10 + t27 * t44;
t7 = -g(3) * t21 + t15 * t20;
t14 = g(1) * t29 - g(2) * t32;
t28 = sin(qJ(2));
t31 = cos(qJ(2));
t33 = -g(3) * t31 + t15 * t28;
t26 = -qJ(3) - pkin(6);
t22 = t31 * pkin(2);
t19 = t22 + pkin(1);
t17 = t26 * t32;
t16 = t32 * t19;
t13 = t21 * t37 + t40;
t11 = -t21 * t39 + t38;
t9 = t14 * t20;
t8 = t15 * t21 + t44;
t6 = t7 * t30;
t5 = t7 * t27;
t4 = -g(1) * t11 - g(2) * t13;
t3 = -g(1) * t10 - g(2) * t12;
t2 = g(1) * t13 - g(2) * t11 + t30 * t44;
t23 = [0, t14, t15, 0, 0, 0, 0, 0, t14 * t31, -t14 * t28, t14 * t21, -t9, -t15, -g(1) * (-t29 * t19 - t17) - g(2) * (-t29 * t26 + t16), 0, 0, 0, 0, 0, t4, t3, t4, t3, t9, -g(1) * (pkin(4) * t38 - t17) - g(2) * (t49 * t32 + t16) + (-g(1) * (-t19 - t49) - g(2) * (pkin(4) * t27 - t26)) * t29; 0, 0, 0, 0, 0, 0, 0, 0, t33, g(3) * t28 + t15 * t31, t7, t8, 0, t33 * pkin(2), 0, 0, 0, 0, 0, t6, -t5, t6, -t5, -t8, -g(3) * (t22 + t49) + t15 * (pkin(2) * t28 + t18 * t20 + t21 * t25); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, t1, t2, 0, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7;];
taug_reg = t23;
