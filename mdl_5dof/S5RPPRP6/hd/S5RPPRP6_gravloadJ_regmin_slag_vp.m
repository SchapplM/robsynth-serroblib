% Calculate minimal parameter regressor of gravitation load for
% S5RPPRP6
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
% taug_reg [5x21]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPPRP6_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP6_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP6_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP6_gravloadJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:55:16
% EndTime: 2019-12-31 17:55:17
% DurationCPUTime: 0.11s
% Computational Cost: add. (86->37), mult. (118->36), div. (0->0), fcn. (109->6), ass. (0->19)
t19 = sin(qJ(1));
t20 = cos(qJ(1));
t28 = -g(1) * t19 + g(2) * t20;
t25 = t20 * pkin(1) + t19 * qJ(2);
t24 = g(2) * t25;
t6 = g(1) * t20 + g(2) * t19;
t15 = pkin(7) + qJ(4);
t10 = cos(t15);
t9 = sin(t15);
t23 = -t9 * pkin(4) + t10 * qJ(5);
t16 = sin(pkin(7));
t21 = pkin(3) * t16 - t23;
t18 = -pkin(6) - qJ(3);
t12 = t20 * qJ(2);
t4 = t6 * t10;
t3 = t6 * t9;
t2 = -g(3) * t9 - t28 * t10;
t1 = g(3) * t10 - t28 * t9;
t5 = [0, -t28, t6, t28, -t6, -g(1) * (-t19 * pkin(1) + t12) - t24, -t6 * t16, -t6 * cos(pkin(7)), -t28, -g(1) * (t12 + (-pkin(1) - qJ(3)) * t19) - g(2) * (t20 * qJ(3) + t25), 0, 0, 0, 0, 0, -t3, -t4, -t3, -t28, t4, -g(1) * t12 - t24 + (-g(1) * t21 + g(2) * t18) * t20 + (-g(1) * (-pkin(1) + t18) - g(2) * t21) * t19; 0, 0, 0, 0, 0, t28, 0, 0, 0, t28, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, t1, -t2, 0, -t1, -g(3) * t23 + t28 * (pkin(4) * t10 + qJ(5) * t9); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2;];
taug_reg = t5;
