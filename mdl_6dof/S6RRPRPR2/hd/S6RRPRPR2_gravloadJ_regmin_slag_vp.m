% Calculate minimal parameter regressor of gravitation load for
% S6RRPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
% 
% Output:
% taug_reg [6x30]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRPR2_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR2_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR2_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR2_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t34 = cos(qJ(2));
t25 = t34 * pkin(2);
t28 = qJ(2) + pkin(10);
t24 = qJ(4) + t28;
t21 = sin(t24);
t22 = cos(t24);
t42 = t22 * pkin(4) + t21 * qJ(5);
t50 = pkin(3) * cos(t28) + t25 + t42;
t32 = sin(qJ(1));
t35 = cos(qJ(1));
t16 = g(1) * t35 + g(2) * t32;
t4 = g(3) * t21 + t16 * t22;
t49 = pkin(4) * t21;
t47 = g(3) * t22;
t30 = sin(qJ(6));
t46 = t32 * t30;
t33 = cos(qJ(6));
t45 = t32 * t33;
t44 = t35 * t30;
t43 = t35 * t33;
t29 = -qJ(3) - pkin(7);
t40 = qJ(5) * t22;
t31 = sin(qJ(2));
t39 = -pkin(3) * sin(t28) - t31 * pkin(2) - t49;
t15 = g(1) * t32 - g(2) * t35;
t38 = pkin(1) + t50;
t36 = -g(3) * t34 + t16 * t31;
t27 = -pkin(8) + t29;
t23 = t25 + pkin(1);
t14 = t35 * t40;
t13 = t32 * t40;
t10 = -t21 * t46 + t43;
t9 = t21 * t45 + t44;
t8 = t21 * t44 + t45;
t7 = t21 * t43 - t46;
t6 = t15 * t22;
t5 = t15 * t21;
t3 = t16 * t21 - t47;
t2 = t4 * t33;
t1 = t4 * t30;
t11 = [0, t15, t16, 0, 0, 0, 0, 0, t15 * t34, -t15 * t31, -t16, -g(1) * (-t32 * t23 - t35 * t29) - g(2) * (t35 * t23 - t32 * t29) 0, 0, 0, 0, 0, t6, -t5, -t16, -t6, t5 (g(1) * t27 - g(2) * t38) * t35 + (g(1) * t38 + g(2) * t27) * t32, 0, 0, 0, 0, 0, -g(1) * t10 - g(2) * t8, g(1) * t9 - g(2) * t7; 0, 0, 0, 0, 0, 0, 0, 0, t36, g(3) * t31 + t16 * t34, 0, t36 * pkin(2), 0, 0, 0, 0, 0, t3, t4, 0, -t3, -t4, -g(1) * (t39 * t35 + t14) - g(2) * (t39 * t32 + t13) - g(3) * t50, 0, 0, 0, 0, 0, -t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, -t3, -t4, -g(1) * (-t35 * t49 + t14) - g(2) * (-t32 * t49 + t13) - g(3) * t42, 0, 0, 0, 0, 0, -t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t7 - g(2) * t9 + t33 * t47, g(1) * t8 - g(2) * t10 - t30 * t47;];
taug_reg  = t11;
