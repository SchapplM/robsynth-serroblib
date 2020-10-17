% Calculate minimal parameter regressor of gravitation load for
% S4PRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% 
% Output:
% taug_reg [4x14]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4PRRR5_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR5_gravloadJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR5_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR5_gravloadJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:33:46
% EndTime: 2019-12-31 16:33:46
% DurationCPUTime: 0.09s
% Computational Cost: add. (59->17), mult. (82->26), div. (0->0), fcn. (86->8), ass. (0->20)
t8 = sin(pkin(7));
t9 = cos(pkin(7));
t14 = g(1) * t9 + g(2) * t8;
t7 = qJ(2) + qJ(3);
t5 = sin(t7);
t6 = cos(t7);
t3 = -g(3) * t6 + t14 * t5;
t20 = g(3) * t5;
t10 = sin(qJ(4));
t18 = t8 * t10;
t12 = cos(qJ(4));
t17 = t8 * t12;
t16 = t9 * t10;
t15 = t9 * t12;
t13 = cos(qJ(2));
t11 = sin(qJ(2));
t4 = t14 * t6 + t20;
t2 = t3 * t12;
t1 = t3 * t10;
t19 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -g(3) * t13 + t14 * t11, g(3) * t11 + t14 * t13, 0, t3, t4, 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t6 * t16 + t17) - g(2) * (-t6 * t18 - t15) + t10 * t20, -g(1) * (-t6 * t15 - t18) - g(2) * (-t6 * t17 + t16) + t12 * t20;];
taug_reg = t19;
