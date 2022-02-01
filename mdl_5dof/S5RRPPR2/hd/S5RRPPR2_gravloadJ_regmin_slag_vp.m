% Calculate minimal parameter regressor of gravitation load for
% S5RRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% 
% Output:
% taug_reg [5x17]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 10:06
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPPR2_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR2_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:05:51
% EndTime: 2022-01-20 10:05:52
% DurationCPUTime: 0.11s
% Computational Cost: add. (149->32), mult. (106->49), div. (0->0), fcn. (110->10), ass. (0->32)
t21 = qJ(1) + qJ(2);
t18 = sin(t21);
t35 = pkin(2) * t18;
t34 = g(3) * sin(pkin(9));
t25 = sin(qJ(1));
t33 = t25 * pkin(1);
t23 = cos(pkin(9));
t24 = sin(qJ(5));
t32 = t23 * t24;
t26 = cos(qJ(5));
t31 = t23 * t26;
t17 = pkin(8) + t21;
t14 = sin(t17);
t15 = cos(t17);
t19 = cos(t21);
t16 = pkin(2) * t19;
t30 = t15 * pkin(3) + t14 * qJ(4) + t16;
t29 = g(1) * t14 - g(2) * t15;
t9 = g(1) * t18 - g(2) * t19;
t28 = -t14 * pkin(3) + t15 * qJ(4) - t35;
t27 = cos(qJ(1));
t20 = t27 * pkin(1);
t10 = g(1) * t19 + g(2) * t18;
t8 = -g(1) * t15 - g(2) * t14;
t7 = t29 * t23;
t6 = t14 * t24 + t15 * t31;
t5 = t14 * t26 - t15 * t32;
t4 = -t14 * t31 + t15 * t24;
t3 = t14 * t32 + t15 * t26;
t2 = -g(1) * t4 - g(2) * t6;
t1 = -g(1) * t3 - g(2) * t5;
t11 = [0, g(1) * t25 - g(2) * t27, g(1) * t27 + g(2) * t25, 0, t9, t10, -g(1) * (-t33 - t35) - g(2) * (t16 + t20), t7, t8, -g(1) * (t28 - t33) - g(2) * (t20 + t30), 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, t9, t10, t9 * pkin(2), t7, t8, -g(1) * t28 - g(2) * t30, 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, 0, -g(3), 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t5 + g(2) * t3 + t24 * t34, g(1) * t6 - g(2) * t4 + t26 * t34;];
taug_reg = t11;
