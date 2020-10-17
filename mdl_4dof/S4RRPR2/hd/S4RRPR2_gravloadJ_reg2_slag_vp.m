% Calculate inertial parameters regressor of gravitation load for
% S4RRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
% 
% Output:
% taug_reg [4x(4*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 18:16
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RRPR2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR2_gravloadJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPR2_gravloadJ_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-19 15:06:30
% EndTime: 2019-07-19 15:06:32
% DurationCPUTime: 0.11s
% Computational Cost: add. (115->29), mult. (96->33), div. (0->0), fcn. (98->6), ass. (0->22)
t19 = sin(qJ(1));
t29 = t19 * pkin(1);
t28 = cos(qJ(4));
t27 = sin(qJ(4));
t18 = qJ(1) + qJ(2);
t15 = sin(t18);
t16 = cos(t18);
t26 = t16 * pkin(2) + t15 * qJ(3);
t20 = cos(qJ(1));
t25 = t20 * pkin(1) + t26;
t11 = t16 * qJ(3);
t24 = -t15 * pkin(2) + t11;
t3 = -t15 * t27 - t16 * t28;
t4 = -t15 * t28 + t16 * t27;
t23 = g(1) * t4 - g(2) * t3;
t2 = g(1) * t3 + g(2) * t4;
t22 = t11 + (-pkin(2) - pkin(3)) * t15;
t21 = g(1) * t19 - g(2) * t20;
t12 = t16 * pkin(3);
t6 = g(1) * t16 + g(2) * t15;
t5 = g(1) * t15 - g(2) * t16;
t1 = [0, 0, 0, 0, 0, 0, t21, g(1) * t20 + g(2) * t19, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, t21 * pkin(1), 0, 0, 0, 0, 0, 0, t5, 0, -t6, -g(1) * (t24 - t29) - g(2) * t25, 0, 0, 0, 0, 0, 0, -t23, t2, 0, -g(1) * (t22 - t29) - g(2) * (t12 + t25); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, 0, 0, 0, t5, 0, -t6, -g(1) * t24 - g(2) * t26, 0, 0, 0, 0, 0, 0, -t23, t2, 0, -g(1) * t22 - g(2) * (t12 + t26); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, -t2, 0, 0;];
taug_reg  = t1;
