% Calculate inertial parameters regressor of gravitation load for
% S4PRRR3
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
% taug_reg [4x(4*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4PRRR3_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR3_gravloadJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR3_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR3_gravloadJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:31:39
% EndTime: 2019-12-31 16:31:39
% DurationCPUTime: 0.08s
% Computational Cost: add. (101->22), mult. (60->25), div. (0->0), fcn. (54->6), ass. (0->16)
t13 = pkin(7) + qJ(2);
t12 = qJ(3) + t13;
t8 = sin(t12);
t9 = cos(t12);
t18 = t9 * pkin(3) + t8 * pkin(6);
t17 = -t8 * pkin(3) + t9 * pkin(6);
t4 = g(1) * t9 + g(2) * t8;
t3 = g(1) * t8 - g(2) * t9;
t10 = sin(t13);
t11 = cos(t13);
t16 = g(1) * t10 - g(2) * t11;
t15 = cos(qJ(4));
t14 = sin(qJ(4));
t2 = t3 * t15;
t1 = t3 * t14;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, g(1) * t11 + g(2) * t10, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, t16 * pkin(2), 0, 0, 0, 0, 0, 0, t2, -t1, -t4, -g(1) * (-pkin(2) * t10 + t17) - g(2) * (pkin(2) * t11 + t18); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, 0, 0, 0, t2, -t1, -t4, -g(1) * t17 - g(2) * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t15 + t4 * t14, g(3) * t14 + t4 * t15, 0, 0;];
taug_reg = t5;
