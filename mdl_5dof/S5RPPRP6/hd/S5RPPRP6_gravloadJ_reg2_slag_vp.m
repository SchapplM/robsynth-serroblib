% Calculate inertial parameters regressor of gravitation load for
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
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPPRP6_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP6_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP6_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP6_gravloadJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:55:16
% EndTime: 2019-12-31 17:55:17
% DurationCPUTime: 0.14s
% Computational Cost: add. (98->40), mult. (136->39), div. (0->0), fcn. (125->6), ass. (0->21)
t22 = sin(qJ(1));
t23 = cos(qJ(1));
t33 = -g(1) * t22 + g(2) * t23;
t19 = sin(pkin(7));
t32 = pkin(3) * t19;
t29 = t23 * pkin(1) + t22 * qJ(2);
t15 = t23 * qJ(2);
t28 = -t22 * pkin(1) + t15;
t6 = g(1) * t23 + g(2) * t22;
t18 = pkin(7) + qJ(4);
t12 = sin(t18);
t13 = cos(t18);
t26 = t12 * pkin(4) - t13 * qJ(5);
t21 = -pkin(6) - qJ(3);
t25 = t22 * t21 + t23 * t32 + t28;
t24 = -t23 * t21 + t22 * t32 + t29;
t4 = t6 * t13;
t3 = t6 * t12;
t2 = -g(3) * t12 - t33 * t13;
t1 = g(3) * t13 - t12 * t33;
t5 = [0, 0, 0, 0, 0, 0, -t33, t6, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, -t6, -g(1) * t28 - g(2) * t29, 0, 0, 0, 0, 0, 0, -t6 * t19, -t6 * cos(pkin(7)), -t33, -g(1) * (t15 + (-pkin(1) - qJ(3)) * t22) - g(2) * (t23 * qJ(3) + t29), 0, 0, 0, 0, 0, 0, -t3, -t4, -t33, -g(1) * t25 - g(2) * t24, 0, 0, 0, 0, 0, 0, -t3, -t33, t4, -g(1) * (t26 * t23 + t25) - g(2) * (t26 * t22 + t24); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, t1, 0, 0, 0, 0, 0, 0, 0, 0, -t2, 0, -t1, g(3) * t26 + t33 * (pkin(4) * t13 + qJ(5) * t12); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2;];
taug_reg = t5;
