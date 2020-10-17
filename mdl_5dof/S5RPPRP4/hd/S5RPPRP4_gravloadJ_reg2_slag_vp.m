% Calculate inertial parameters regressor of gravitation load for
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
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPPRP4_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP4_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP4_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP4_gravloadJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:52:20
% EndTime: 2019-12-31 17:52:21
% DurationCPUTime: 0.15s
% Computational Cost: add. (93->37), mult. (182->42), div. (0->0), fcn. (207->6), ass. (0->23)
t32 = cos(qJ(1));
t31 = sin(qJ(1));
t30 = t32 * pkin(1) + t31 * qJ(2);
t29 = cos(pkin(7));
t28 = sin(pkin(7));
t27 = t32 * pkin(2) + t30;
t7 = -t31 * t28 - t32 * t29;
t8 = t32 * t28 - t31 * t29;
t26 = g(1) * t8 - g(2) * t7;
t25 = g(1) * t7 + g(2) * t8;
t24 = -t31 * pkin(1) + t32 * qJ(2);
t23 = -t31 * pkin(2) + t24;
t21 = sin(qJ(4));
t22 = cos(qJ(4));
t1 = g(3) * t22 - t25 * t21;
t20 = -qJ(5) - pkin(6);
t14 = t22 * pkin(4) + pkin(3);
t10 = g(1) * t32 + g(2) * t31;
t9 = g(1) * t31 - g(2) * t32;
t4 = t26 * t22;
t3 = t26 * t21;
t2 = -g(3) * t21 - t25 * t22;
t5 = [0, 0, 0, 0, 0, 0, t9, t10, 0, 0, 0, 0, 0, 0, 0, 0, t9, 0, -t10, -g(1) * t24 - g(2) * t30, 0, 0, 0, 0, 0, 0, -t26, t25, 0, -g(1) * t23 - g(2) * t27, 0, 0, 0, 0, 0, 0, -t4, t3, -t25, -g(1) * (t8 * pkin(3) + t7 * pkin(6) + t23) - g(2) * (-t7 * pkin(3) + t8 * pkin(6) + t27), 0, 0, 0, 0, 0, 0, -t4, t3, -t25, -g(1) * (t8 * t14 - t7 * t20 + t23) - g(2) * (-t7 * t14 - t8 * t20 + t27); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26;];
taug_reg = t5;
