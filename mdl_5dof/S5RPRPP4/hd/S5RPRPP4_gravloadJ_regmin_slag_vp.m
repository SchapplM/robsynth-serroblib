% Calculate minimal parameter regressor of gravitation load for
% S5RPRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta4]';
% 
% Output:
% taug_reg [5x19]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRPP4_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP4_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP4_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP4_gravloadJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:14:50
% EndTime: 2019-12-31 18:14:50
% DurationCPUTime: 0.11s
% Computational Cost: add. (83->34), mult. (127->41), div. (0->0), fcn. (115->6), ass. (0->19)
t16 = sin(qJ(1));
t18 = cos(qJ(1));
t29 = -g(1) * t16 + g(2) * t18;
t15 = sin(qJ(3));
t26 = t15 * pkin(3);
t25 = t18 * pkin(1) + t16 * qJ(2);
t24 = -t16 * pkin(1) + t18 * qJ(2);
t13 = qJ(3) + pkin(7);
t7 = sin(t13);
t8 = cos(t13);
t23 = t7 * pkin(4) - t8 * qJ(5);
t3 = g(1) * t18 + g(2) * t16;
t14 = -qJ(4) - pkin(6);
t22 = t16 * t14 + t18 * t26 + t24;
t21 = -t18 * t14 + t16 * t26 + t25;
t17 = cos(qJ(3));
t19 = g(3) * t15 + t17 * t29;
t1 = -g(3) * t7 - t29 * t8;
t2 = [0, -t29, t3, t29, -t3, -g(1) * t24 - g(2) * t25, 0, 0, 0, 0, 0, -t3 * t15, -t3 * t17, -t29, -g(1) * t22 - g(2) * t21, -t3 * t7, -t29, t3 * t8, -g(1) * (t23 * t18 + t22) - g(2) * (t23 * t16 + t21); 0, 0, 0, 0, 0, t29, 0, 0, 0, 0, 0, 0, 0, 0, t29, 0, 0, 0, t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, g(3) * t17 - t15 * t29, 0, t19 * pkin(3), -t1, 0, -g(3) * t8 + t29 * t7, -g(3) * (-t23 - t26) + t29 * (pkin(3) * t17 + pkin(4) * t8 + qJ(5) * t7); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, 0, 0, 0, -t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1;];
taug_reg = t2;
