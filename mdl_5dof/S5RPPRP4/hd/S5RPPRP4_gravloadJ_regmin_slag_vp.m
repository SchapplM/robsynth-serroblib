% Calculate minimal parameter regressor of gravitation load for
% S5RPPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
% 
% Output:
% taug_reg [5x18]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 17:13
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPPRP4_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP4_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP4_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP4_gravloadJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 17:13:15
% EndTime: 2021-01-15 17:13:16
% DurationCPUTime: 0.11s
% Computational Cost: add. (74->31), mult. (144->40), div. (0->0), fcn. (155->6), ass. (0->24)
t16 = sin(pkin(7));
t17 = cos(pkin(7));
t18 = qJ(5) + pkin(6);
t21 = cos(qJ(4));
t27 = pkin(4) * t21 + pkin(3);
t29 = t27 * t16 - t18 * t17 + qJ(2);
t19 = sin(qJ(4));
t20 = sin(qJ(1));
t22 = cos(qJ(1));
t6 = -t20 * t16 - t22 * t17;
t7 = t22 * t16 - t20 * t17;
t24 = -g(1) * t6 - g(2) * t7;
t1 = g(3) * t21 + t24 * t19;
t23 = pkin(1) + pkin(2);
t25 = g(1) * t7 - g(2) * t6;
t14 = t22 * qJ(2);
t13 = t20 * qJ(2);
t9 = g(1) * t22 + g(2) * t20;
t8 = g(1) * t20 - g(2) * t22;
t5 = t18 * t16 + t27 * t17 + t23;
t4 = t25 * t21;
t3 = t25 * t19;
t2 = -g(3) * t19 + t24 * t21;
t10 = [0, t8, t9, t8, -t9, -g(1) * (-t20 * pkin(1) + t14) - g(2) * (t22 * pkin(1) + t13), -g(1) * (-t23 * t20 + t14) - g(2) * (t23 * t22 + t13), 0, 0, 0, 0, 0, -t4, t3, -t4, t3, t24, -g(1) * (-t5 * t20 + t29 * t22) - g(2) * (t29 * t20 + t5 * t22); 0, 0, 0, 0, 0, -t8, -t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8; 0, 0, 0, 0, 0, 0, g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, t1, t2, 0, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t25;];
taug_reg = t10;
