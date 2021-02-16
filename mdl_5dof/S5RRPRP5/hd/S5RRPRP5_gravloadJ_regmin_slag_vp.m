% Calculate minimal parameter regressor of gravitation load for
% S5RRPRP5
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
% Datum: 2021-01-15 20:19
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPRP5_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP5_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP5_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP5_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 20:18:18
% EndTime: 2021-01-15 20:18:19
% DurationCPUTime: 0.15s
% Computational Cost: add. (197->43), mult. (177->55), div. (0->0), fcn. (164->8), ass. (0->30)
t23 = qJ(2) + pkin(8);
t18 = cos(t23);
t27 = cos(qJ(2));
t20 = t27 * pkin(2);
t19 = qJ(4) + t23;
t14 = sin(t19);
t15 = cos(t19);
t34 = t15 * pkin(4) + t14 * qJ(5);
t36 = pkin(3) * t18 + t20 + t34;
t35 = pkin(4) * t14;
t24 = -qJ(3) - pkin(6);
t32 = qJ(5) * t15;
t17 = sin(t23);
t25 = sin(qJ(2));
t31 = -t25 * pkin(2) - pkin(3) * t17 - t35;
t26 = sin(qJ(1));
t28 = cos(qJ(1));
t10 = g(1) * t28 + g(2) * t26;
t9 = g(1) * t26 - g(2) * t28;
t30 = pkin(1) + t36;
t29 = -g(3) * t27 + t10 * t25;
t22 = -pkin(7) + t24;
t16 = t20 + pkin(1);
t8 = t28 * t32;
t7 = t26 * t32;
t4 = t9 * t15;
t3 = t9 * t14;
t2 = g(3) * t14 + t10 * t15;
t1 = -g(3) * t15 + t10 * t14;
t5 = [0, t9, t10, 0, 0, 0, 0, 0, t9 * t27, -t9 * t25, t9 * t18, -t9 * t17, -t10, -g(1) * (-t26 * t16 - t24 * t28) - g(2) * (t28 * t16 - t26 * t24), 0, 0, 0, 0, 0, t4, -t3, t4, -t10, t3, (g(1) * t22 - g(2) * t30) * t28 + (g(1) * t30 + g(2) * t22) * t26; 0, 0, 0, 0, 0, 0, 0, 0, t29, g(3) * t25 + t10 * t27, -g(3) * t18 + t10 * t17, g(3) * t17 + t10 * t18, 0, t29 * pkin(2), 0, 0, 0, 0, 0, t1, t2, t1, 0, -t2, -g(1) * (t31 * t28 + t8) - g(2) * (t31 * t26 + t7) - g(3) * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, t1, 0, -t2, -g(1) * (-t28 * t35 + t8) - g(2) * (-t26 * t35 + t7) - g(3) * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t5;
