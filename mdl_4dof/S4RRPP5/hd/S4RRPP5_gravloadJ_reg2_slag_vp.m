% Calculate inertial parameters regressor of gravitation load for
% S4RRPP5
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
% Datum: 2019-12-31 17:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RRPP5_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP5_gravloadJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP5_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP5_gravloadJ_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t16 = sin(qJ(2));
t18 = cos(qJ(2));
t26 = t18 * pkin(2) + t16 * qJ(3);
t17 = sin(qJ(1));
t19 = cos(qJ(1));
t5 = g(1) * t19 + g(2) * t17;
t32 = t5 * t16;
t31 = pkin(2) * t16;
t30 = g(1) * t17;
t25 = t19 * pkin(1) + t17 * pkin(5);
t24 = qJ(3) * t18;
t23 = t18 * qJ(4);
t22 = t26 * t19 + t25;
t21 = -g(2) * t19 + t30;
t20 = -pkin(1) - t26;
t13 = t19 * pkin(5);
t8 = t19 * t24;
t6 = t17 * t24;
t4 = t21 * t18;
t3 = t21 * t16;
t2 = g(3) * t16 + t5 * t18;
t1 = -g(3) * t18 + t32;
t7 = [0, 0, 0, 0, 0, 0, t21, t5, 0, 0, 0, 0, 0, 0, 0, 0, t4, -t3, -t5, -g(1) * (-t17 * pkin(1) + t13) - g(2) * t25, 0, 0, 0, 0, 0, 0, -t5, -t4, t3, -g(1) * t13 - g(2) * t22 - t20 * t30, 0, 0, 0, 0, 0, 0, -t5, t3, t4, -g(1) * (t19 * pkin(3) + t13) - g(2) * (t19 * t23 + t22) + (-g(1) * (t20 - t23) - g(2) * pkin(3)) * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, -t2, -g(1) * (-t19 * t31 + t8) - g(2) * (-t17 * t31 + t6) - g(3) * t26, 0, 0, 0, 0, 0, 0, 0, -t2, t1, -g(1) * t8 - g(2) * t6 - g(3) * (t23 + t26) + (pkin(2) + qJ(4)) * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2;];
taug_reg = t7;
