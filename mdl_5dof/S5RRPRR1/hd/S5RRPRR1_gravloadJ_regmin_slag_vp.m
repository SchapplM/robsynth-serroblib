% Calculate minimal parameter regressor of gravitation load for
% S5RRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d4,d5]';
% 
% Output:
% taug_reg [5x26]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 17:22
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPRR1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR1_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S5RRPRR1_gravloadJ_regmin_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t16 = sin(qJ(1));
t19 = cos(qJ(1));
t10 = g(1) * t19 + g(2) * t16;
t13 = qJ(2) + qJ(4);
t11 = sin(t13);
t12 = cos(t13);
t3 = -g(3) * t12 + t10 * t11;
t18 = cos(qJ(2));
t27 = pkin(1) * t18;
t26 = g(3) * t11;
t14 = sin(qJ(5));
t24 = t16 * t14;
t17 = cos(qJ(5));
t23 = t16 * t17;
t22 = t19 * t14;
t21 = t19 * t17;
t9 = g(1) * t16 - g(2) * t19;
t15 = sin(qJ(2));
t20 = -g(3) * t18 + t10 * t15;
t8 = t12 * t21 + t24;
t7 = -t12 * t22 + t23;
t6 = -t12 * t23 + t22;
t5 = t12 * t24 + t21;
t4 = t10 * t12 + t26;
t2 = t3 * t17;
t1 = t3 * t14;
t25 = [0, t9, t10, 0, 0, 0, 0, 0, t9 * t18, -t9 * t15, -t10, -g(1) * (t19 * qJ(3) - t16 * t27) - g(2) * (t16 * qJ(3) + t19 * t27), 0, 0, 0, 0, 0, t9 * t12, -t9 * t11, 0, 0, 0, 0, 0, -g(1) * t6 - g(2) * t8, -g(1) * t5 - g(2) * t7; 0, 0, 0, 0, 0, 0, 0, 0, t20, g(3) * t15 + t10 * t18, 0, t20 * pkin(1), 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t7 + g(2) * t5 + t14 * t26, g(1) * t8 - g(2) * t6 + t17 * t26;];
taug_reg  = t25;
