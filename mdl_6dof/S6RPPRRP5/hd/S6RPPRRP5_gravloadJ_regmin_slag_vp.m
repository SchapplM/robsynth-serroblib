% Calculate minimal parameter regressor of gravitation load for
% S6RPPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5]';
% 
% Output:
% taug_reg [6x25]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPPRRP5_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP5_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP5_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP5_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:58:26
% EndTime: 2019-05-05 14:58:26
% DurationCPUTime: 0.17s
% Computational Cost: add. (80->45), mult. (175->56), div. (0->0), fcn. (171->6), ass. (0->30)
t19 = sin(qJ(1));
t22 = cos(qJ(1));
t9 = g(1) * t22 + g(2) * t19;
t17 = sin(qJ(5));
t18 = sin(qJ(4));
t20 = cos(qJ(5));
t30 = t22 * t20;
t33 = t19 * t17;
t3 = t18 * t33 - t30;
t21 = cos(qJ(4));
t34 = g(3) * t21;
t31 = t22 * t17;
t32 = t19 * t20;
t5 = -t18 * t31 - t32;
t40 = -g(1) * t5 + g(2) * t3 + t17 * t34;
t2 = -g(3) * t18 + t9 * t21;
t29 = -pkin(1) - qJ(3);
t28 = t22 * pkin(1) + t19 * qJ(2);
t26 = pkin(5) * t17 + pkin(7);
t25 = g(2) * (t22 * qJ(3) + t28);
t8 = g(1) * t19 - g(2) * t22;
t10 = t20 * pkin(5) + pkin(4);
t16 = -qJ(6) - pkin(8);
t23 = t18 * t10 + t21 * t16;
t13 = t22 * qJ(2);
t7 = t8 * t21;
t6 = t18 * t30 - t33;
t4 = -t18 * t32 - t31;
t1 = t9 * t18 + t34;
t11 = [0, t8, t9, -t8, -t9, -g(1) * (-t19 * pkin(1) + t13) - g(2) * t28, -t9, t8, -g(1) * (t29 * t19 + t13) - t25, 0, 0, 0, 0, 0, t8 * t18, t7, 0, 0, 0, 0, 0, -g(1) * t4 - g(2) * t6, -g(1) * t3 - g(2) * t5, -t7, -g(1) * t13 - t25 + (g(1) * t26 - g(2) * t23) * t22 + (-g(1) * (-t23 + t29) + g(2) * t26) * t19; 0, 0, 0, 0, 0, -t8, 0, 0, -t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8; 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, t1, 0, 0, 0, 0, 0, -t2 * t20, t2 * t17, -t1, g(3) * t23 - t9 * (t10 * t21 - t16 * t18); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, g(1) * t6 - g(2) * t4 + t20 * t34, 0, t40 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2;];
taug_reg  = t11;
