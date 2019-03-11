% Calculate minimal parameter regressor of gravitation load for
% S6RPRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
% 
% Output:
% taug_reg [6x32]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRPRR5_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR5_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR5_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR5_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t27 = sin(qJ(6));
t23 = pkin(10) + qJ(3);
t20 = sin(t23);
t21 = cos(t23);
t28 = sin(qJ(5));
t42 = cos(qJ(5));
t34 = t20 * t28 + t21 * t42;
t29 = sin(qJ(1));
t46 = t20 * t42 - t21 * t28;
t5 = t46 * t29;
t31 = cos(qJ(1));
t7 = t46 * t31;
t33 = -g(1) * t7 - g(2) * t5 + g(3) * t34;
t48 = t33 * t27;
t30 = cos(qJ(6));
t47 = t33 * t30;
t16 = g(1) * t31 + g(2) * t29;
t43 = g(3) * t46;
t15 = g(1) * t29 - g(2) * t31;
t6 = t34 * t29;
t39 = t6 * t27 - t31 * t30;
t38 = t31 * t27 + t6 * t30;
t37 = t21 * pkin(3) + t20 * qJ(4);
t25 = cos(pkin(10));
t35 = t25 * pkin(2) + pkin(1) + t37;
t8 = t34 * t31;
t32 = g(1) * t8 + g(2) * t6 + t43;
t26 = -pkin(7) - qJ(2);
t10 = t15 * t21;
t9 = t15 * t20;
t4 = g(3) * t20 + t16 * t21;
t3 = -g(3) * t21 + t16 * t20;
t2 = -t29 * t27 + t8 * t30;
t1 = -t8 * t27 - t29 * t30;
t11 = [0, t15, t16, t15 * t25, -t15 * sin(pkin(10)) -t16, -g(1) * (-t29 * pkin(1) + t31 * qJ(2)) - g(2) * (t31 * pkin(1) + t29 * qJ(2)) 0, 0, 0, 0, 0, t10, -t9, t10, -t16, t9 (g(1) * t26 - g(2) * t35) * t31 + (g(1) * t35 + g(2) * t26) * t29, 0, 0, 0, 0, 0, g(1) * t6 - g(2) * t8, g(1) * t5 - g(2) * t7, 0, 0, 0, 0, 0, g(1) * t38 - g(2) * t2, -g(1) * t39 - g(2) * t1; 0, 0, 0, 0, 0, 0, -t15, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, t3, 0, -t4, -g(3) * t37 + t16 * (pkin(3) * t20 - qJ(4) * t21) 0, 0, 0, 0, 0, -t33, -t32, 0, 0, 0, 0, 0, -t47, t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, t32, 0, 0, 0, 0, 0, t47, -t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t1 + g(2) * t39 + t27 * t43, g(1) * t2 + g(2) * t38 + t30 * t43;];
taug_reg  = t11;
