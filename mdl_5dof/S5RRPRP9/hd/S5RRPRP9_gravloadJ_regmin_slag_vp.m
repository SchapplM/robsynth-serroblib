% Calculate minimal parameter regressor of gravitation load for
% S5RRPRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% taug_reg [5x25]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPRP9_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP9_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP9_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP9_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:07:28
% EndTime: 2019-12-31 20:07:28
% DurationCPUTime: 0.23s
% Computational Cost: add. (204->62), mult. (314->90), div. (0->0), fcn. (331->8), ass. (0->39)
t23 = cos(pkin(8));
t14 = t23 * pkin(3) + pkin(2);
t24 = -pkin(7) - qJ(3);
t25 = sin(qJ(2));
t27 = cos(qJ(2));
t49 = t27 * t14 - t25 * t24;
t26 = sin(qJ(1));
t28 = cos(qJ(1));
t34 = g(1) * t28 + g(2) * t26;
t9 = -g(3) * t27 + t34 * t25;
t48 = g(1) * t26;
t45 = g(3) * t25;
t42 = t26 * t27;
t21 = pkin(8) + qJ(4);
t15 = sin(t21);
t40 = t28 * t15;
t16 = cos(t21);
t39 = t28 * t16;
t22 = sin(pkin(8));
t38 = t28 * t22;
t37 = t28 * t23;
t36 = t28 * pkin(1) + t26 * pkin(6);
t5 = t15 * t42 + t39;
t7 = -t26 * t16 + t27 * t40;
t35 = g(1) * t5 - g(2) * t7;
t33 = -g(2) * t28 + t48;
t32 = t27 * pkin(2) + t25 * qJ(3);
t30 = pkin(4) * t16 + qJ(5) * t15 + t14;
t1 = g(1) * t7 + g(2) * t5 + t15 * t45;
t6 = t16 * t42 - t40;
t8 = t26 * t15 + t27 * t39;
t29 = g(1) * t8 + g(2) * t6 + t16 * t45;
t18 = t28 * pkin(6);
t11 = t33 * t25;
t10 = t34 * t27 + t45;
t4 = t9 * t16;
t3 = t9 * t15;
t2 = g(1) * t6 - g(2) * t8;
t12 = [0, t33, t34, 0, 0, 0, 0, 0, t33 * t27, -t11, -g(1) * (-t23 * t42 + t38) - g(2) * (t26 * t22 + t27 * t37), -g(1) * (t22 * t42 + t37) - g(2) * (t26 * t23 - t27 * t38), t11, -g(1) * t18 - g(2) * (t32 * t28 + t36) - (-pkin(1) - t32) * t48, 0, 0, 0, 0, 0, t2, -t35, t2, t11, t35, -g(1) * (pkin(3) * t38 - t6 * pkin(4) - t5 * qJ(5) + t18) - g(2) * (t8 * pkin(4) + t7 * qJ(5) + t49 * t28 + t36) + (-g(1) * (-pkin(1) - t49) - g(2) * pkin(3) * t22) * t26; 0, 0, 0, 0, 0, 0, 0, 0, t9, t10, t9 * t23, -t9 * t22, -t10, -g(3) * t32 + t34 * (pkin(2) * t25 - qJ(3) * t27), 0, 0, 0, 0, 0, t4, -t3, t4, -t10, t3, (-g(3) * t30 + t34 * t24) * t27 + (g(3) * t24 + t34 * t30) * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t29, t1, 0, -t29, -g(1) * (-t7 * pkin(4) + t8 * qJ(5)) - g(2) * (-t5 * pkin(4) + t6 * qJ(5)) - (-pkin(4) * t15 + qJ(5) * t16) * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t12;
