% Calculate inertial parameters regressor of gravitation load for
% S4PRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1]';
% 
% Output:
% taug_reg [4x(4*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4PRPR7_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR7_gravloadJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR7_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR7_gravloadJ_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:25:55
% EndTime: 2019-12-31 16:25:55
% DurationCPUTime: 0.12s
% Computational Cost: add. (43->31), mult. (102->40), div. (0->0), fcn. (100->6), ass. (0->19)
t11 = sin(qJ(2));
t8 = sin(pkin(6));
t9 = cos(pkin(6));
t15 = g(1) * t9 + g(2) * t8;
t26 = t15 * t11;
t13 = cos(qJ(2));
t2 = g(3) * t11 + t15 * t13;
t22 = t13 * pkin(2) + t11 * qJ(3);
t21 = pkin(2) * t11;
t19 = g(3) * t13;
t10 = sin(qJ(4));
t18 = t10 * t11;
t12 = cos(qJ(4));
t17 = t11 * t12;
t16 = qJ(3) * t13;
t4 = t9 * t16;
t3 = t8 * t16;
t1 = -t19 + t26;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, -t2, -g(1) * (-t9 * t21 + t4) - g(2) * (-t8 * t21 + t3) - g(3) * t22, 0, 0, 0, 0, 0, 0, -t2 * t10, -t2 * t12, t1, -g(1) * t4 - g(2) * t3 - g(3) * (t13 * pkin(5) + t22) + (pkin(2) + pkin(5)) * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t8 * t10 + t9 * t17) - g(2) * (t9 * t10 + t8 * t17) + t12 * t19, -g(1) * (-t8 * t12 - t9 * t18) - g(2) * (t9 * t12 - t8 * t18) - t10 * t19, 0, 0;];
taug_reg = t5;
