% Calculate minimal parameter regressor of gravitation load for
% S5PRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% 
% Output:
% taug_reg [5x16]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 15:14
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRPRP3_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP3_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP3_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP3_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 15:14:00
% EndTime: 2021-01-15 15:14:01
% DurationCPUTime: 0.11s
% Computational Cost: add. (96->28), mult. (140->39), div. (0->0), fcn. (139->8), ass. (0->25)
t10 = sin(pkin(7));
t11 = cos(pkin(7));
t20 = g(1) * t11 + g(2) * t10;
t9 = qJ(2) + pkin(8);
t7 = sin(t9);
t8 = cos(t9);
t29 = -g(3) * t8 + t20 * t7;
t28 = g(3) * t7;
t13 = sin(qJ(4));
t24 = t10 * t13;
t15 = cos(qJ(4));
t23 = t10 * t15;
t22 = t11 * t13;
t21 = t11 * t15;
t14 = sin(qJ(2));
t16 = cos(qJ(2));
t17 = -g(3) * t16 + t20 * t14;
t1 = -g(1) * (-t8 * t22 + t23) - g(2) * (-t8 * t24 - t21) + t13 * t28;
t12 = -qJ(5) - pkin(6);
t6 = t15 * pkin(4) + pkin(3);
t5 = -g(1) * t10 + g(2) * t11;
t4 = t29 * t15;
t3 = t29 * t13;
t2 = -g(1) * (-t8 * t21 - t24) - g(2) * (-t8 * t23 + t22) + t15 * t28;
t18 = [-g(3), 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, t17, g(3) * t14 + t20 * t16, t17 * pkin(2), 0, 0, 0, 0, 0, t4, -t3, t4, -t3, -t20 * t8 - t28, -g(3) * (t16 * pkin(2) - t7 * t12 + t8 * t6) + t20 * (pkin(2) * t14 + t12 * t8 + t6 * t7); 0, 0, 0, 0, t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, t1, t2, 0, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29;];
taug_reg = t18;
