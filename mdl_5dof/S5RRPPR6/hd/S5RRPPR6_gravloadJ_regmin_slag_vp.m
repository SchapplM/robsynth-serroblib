% Calculate minimal parameter regressor of gravitation load for
% S5RRPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% 
% Output:
% taug_reg [5x25]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 19:48
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPPR6_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR6_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR6_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR6_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 19:47:31
% EndTime: 2021-01-15 19:47:33
% DurationCPUTime: 0.20s
% Computational Cost: add. (156->51), mult. (202->77), div. (0->0), fcn. (206->12), ass. (0->43)
t30 = sin(qJ(1));
t32 = cos(qJ(1));
t13 = g(1) * t32 + g(2) * t30;
t23 = qJ(2) + pkin(8);
t17 = sin(t23);
t19 = cos(t23);
t45 = g(3) * t19;
t35 = t13 * t17 - t45;
t46 = g(3) * t17;
t22 = pkin(9) + qJ(5);
t16 = sin(t22);
t44 = t30 * t16;
t18 = cos(t22);
t43 = t30 * t18;
t24 = sin(pkin(9));
t42 = t30 * t24;
t26 = cos(pkin(9));
t41 = t30 * t26;
t28 = -qJ(3) - pkin(6);
t40 = t30 * t28;
t39 = t32 * t16;
t38 = t32 * t18;
t37 = t32 * t24;
t36 = t32 * t26;
t12 = g(1) * t30 - g(2) * t32;
t29 = sin(qJ(2));
t31 = cos(qJ(2));
t33 = -g(3) * t31 + t13 * t29;
t27 = cos(pkin(8));
t25 = sin(pkin(8));
t20 = t31 * pkin(2);
t15 = t20 + pkin(1);
t14 = t28 * t32;
t11 = -t25 * pkin(3) + qJ(4) * t27;
t10 = pkin(3) * t27 + qJ(4) * t25 + pkin(2);
t8 = t12 * t17;
t7 = t19 * t38 + t44;
t6 = -t19 * t39 + t43;
t5 = -t19 * t43 + t39;
t4 = t19 * t44 + t38;
t3 = t13 * t19 + t46;
t1 = t10 * t31 + t11 * t29 + pkin(1);
t2 = [0, t12, t13, 0, 0, 0, 0, 0, t12 * t31, -t12 * t29, t12 * t19, -t8, -t13, -g(1) * (-t30 * t15 - t14) - g(2) * (t32 * t15 - t40), -g(1) * (-t19 * t41 + t37) - g(2) * (t19 * t36 + t42), -g(1) * (t19 * t42 + t36) - g(2) * (-t19 * t37 + t41), t8, -g(1) * (-t1 * t30 - t14) - g(2) * (t1 * t32 - t40), 0, 0, 0, 0, 0, -g(1) * t5 - g(2) * t7, -g(1) * t4 - g(2) * t6; 0, 0, 0, 0, 0, 0, 0, 0, t33, g(3) * t29 + t13 * t31, t35, t3, 0, t33 * pkin(2), t35 * t26, -t35 * t24, -t3, -g(3) * (t19 * pkin(3) + t17 * qJ(4) + t20) - t13 * (-t10 * t29 + t11 * t31), 0, 0, 0, 0, 0, t35 * t18, -t35 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, 0, 0, 0, -t12, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45 - t13 * (t25 * t31 + t27 * t29), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t6 + g(2) * t4 + t16 * t46, g(1) * t7 - g(2) * t5 + t18 * t46;];
taug_reg = t2;
