% Calculate minimal parameter regressor of gravitation load for
% S6RPRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% 
% Output:
% taug_reg [6x29]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRRP1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP1_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP1_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 01:09:27
% EndTime: 2019-05-06 01:09:27
% DurationCPUTime: 0.25s
% Computational Cost: add. (358->62), mult. (334->87), div. (0->0), fcn. (347->10), ass. (0->44)
t28 = qJ(3) + qJ(4);
t24 = sin(t28);
t27 = qJ(1) + pkin(10);
t22 = sin(t27);
t23 = cos(t27);
t43 = g(1) * t23 + g(2) * t22;
t38 = t43 * t24;
t25 = cos(t28);
t54 = t25 * pkin(4) + t24 * pkin(9);
t30 = sin(qJ(3));
t53 = pkin(3) * t30;
t52 = pkin(9) * t25;
t49 = g(3) * t24;
t29 = sin(qJ(5));
t48 = t25 * t29;
t32 = cos(qJ(5));
t47 = t25 * t32;
t46 = qJ(6) * t29;
t45 = pkin(5) * t47 + t25 * t46 + t54;
t10 = -t22 * t32 + t23 * t48;
t8 = t22 * t48 + t23 * t32;
t44 = g(1) * t8 - g(2) * t10;
t42 = g(1) * t22 - g(2) * t23;
t31 = sin(qJ(1));
t34 = cos(qJ(1));
t41 = g(1) * t31 - g(2) * t34;
t33 = cos(qJ(3));
t26 = t33 * pkin(3);
t40 = t26 + pkin(2) + t54;
t1 = g(1) * t10 + g(2) * t8 + t29 * t49;
t11 = t22 * t29 + t23 * t47;
t9 = t22 * t47 - t23 * t29;
t37 = g(1) * t11 + g(2) * t9 + t32 * t49;
t5 = -g(3) * t25 + t38;
t36 = (pkin(5) * t32 + pkin(4) + t46) * t38;
t35 = -pkin(8) - pkin(7);
t13 = t23 * t52;
t12 = t22 * t52;
t7 = t42 * t24;
t6 = t43 * t25 + t49;
t4 = t5 * t32;
t3 = -g(3) * t48 + t29 * t38;
t2 = g(1) * t9 - g(2) * t11;
t14 = [0, t41, g(1) * t34 + g(2) * t31, t41 * pkin(1), 0, 0, 0, 0, 0, t42 * t33, -t42 * t30, 0, 0, 0, 0, 0, t42 * t25, -t7, 0, 0, 0, 0, 0, t2, -t44, t2, t7, t44, -g(1) * (-t31 * pkin(1) - t9 * pkin(5) - t8 * qJ(6)) - g(2) * (t34 * pkin(1) + t11 * pkin(5) + t10 * qJ(6)) + (g(1) * t35 - g(2) * t40) * t23 + (g(1) * t40 + g(2) * t35) * t22; 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t33 + t43 * t30, g(3) * t30 + t43 * t33, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, t4, -t3, t4, -t6, t3, -g(1) * (-t23 * t53 + t13) - g(2) * (-t22 * t53 + t12) - g(3) * (t26 + t45) + t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, t4, -t3, t4, -t6, t3, -g(1) * t13 - g(2) * t12 - g(3) * t45 + t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t37, t1, 0, -t37, -g(1) * (-t10 * pkin(5) + t11 * qJ(6)) - g(2) * (-t8 * pkin(5) + t9 * qJ(6)) - (-pkin(5) * t29 + qJ(6) * t32) * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t14;
