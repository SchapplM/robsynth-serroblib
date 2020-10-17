% Calculate minimal parameter regressor of gravitation load for
% S5RPPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
% 
% Output:
% taug_reg [5x16]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPPRP3_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP3_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP3_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP3_gravloadJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:51:02
% EndTime: 2019-12-31 17:51:02
% DurationCPUTime: 0.10s
% Computational Cost: add. (70->25), mult. (68->31), div. (0->0), fcn. (59->6), ass. (0->16)
t12 = sin(qJ(4));
t20 = pkin(4) * t12;
t15 = cos(qJ(1));
t10 = qJ(1) + pkin(7);
t7 = sin(t10);
t8 = cos(t10);
t19 = t15 * pkin(1) + t8 * pkin(2) + t7 * qJ(3);
t13 = sin(qJ(1));
t18 = -t13 * pkin(1) + t8 * qJ(3);
t2 = -g(1) * t8 - g(2) * t7;
t1 = g(1) * t7 - g(2) * t8;
t17 = g(1) * t13 - g(2) * t15;
t14 = cos(qJ(4));
t16 = g(3) * t12 - t1 * t14;
t11 = -qJ(5) - pkin(6);
t3 = [0, t17, g(1) * t15 + g(2) * t13, t17 * pkin(1), -t1, t2, -g(1) * (-t7 * pkin(2) + t18) - g(2) * t19, 0, 0, 0, 0, 0, t2 * t12, t2 * t14, t1, -g(1) * (t8 * t20 + (-pkin(2) + t11) * t7 + t18) - g(2) * (-t8 * t11 + t7 * t20 + t19); 0, 0, 0, -g(3), 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, -t1, 0, 0, 0, 0, 0, 0, 0, 0, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, g(3) * t14 + t1 * t12, 0, t16 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2;];
taug_reg = t3;
