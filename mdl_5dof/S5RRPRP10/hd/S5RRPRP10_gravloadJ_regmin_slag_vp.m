% Calculate minimal parameter regressor of gravitation load for
% S5RRPRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% 
% Output:
% taug_reg [5x23]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPRP10_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP10_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP10_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP10_gravloadJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t23 = sin(qJ(2));
t15 = t23 * qJ(3);
t26 = cos(qJ(2));
t36 = t26 * pkin(2) + t15;
t24 = sin(qJ(1));
t27 = cos(qJ(1));
t9 = g(1) * t27 + g(2) * t24;
t25 = cos(qJ(4));
t37 = t27 * t25;
t22 = sin(qJ(4));
t42 = t24 * t22;
t3 = t23 * t37 - t42;
t45 = g(3) * t26;
t38 = t27 * t22;
t41 = t24 * t25;
t5 = t23 * t41 + t38;
t54 = -g(1) * t3 - g(2) * t5 + t25 * t45;
t2 = g(3) * t23 + t9 * t26;
t51 = pkin(2) * t23;
t50 = pkin(4) * t22;
t49 = g(1) * t24;
t43 = t22 * t26;
t21 = -qJ(5) - pkin(7);
t40 = t26 * t21;
t39 = t26 * t27;
t35 = qJ(3) * t26;
t34 = pkin(4) * t43;
t32 = t23 * t38;
t31 = pkin(2) * t39 + t24 * pkin(6) + (pkin(1) + t15) * t27;
t30 = -g(2) * t27 + t49;
t29 = -pkin(1) - t36;
t18 = t27 * pkin(6);
t14 = t25 * pkin(4) + pkin(3);
t12 = t27 * t35;
t10 = t24 * t35;
t8 = t30 * t26;
t7 = t30 * t23;
t6 = -t23 * t42 + t37;
t4 = t32 + t41;
t1 = t9 * t23 - t45;
t11 = [0, t30, t9, 0, 0, 0, 0, 0, t8, -t7, -t9, -t8, t7, -g(1) * t18 - g(2) * t31 - t29 * t49, 0, 0, 0, 0, 0, -g(1) * t6 - g(2) * t4, g(1) * t5 - g(2) * t3, t8, -g(1) * (t27 * t14 + t18) - g(2) * (pkin(4) * t32 - t21 * t39 + t31) + (-g(1) * (-t23 * t50 + t29 + t40) - g(2) * t14) * t24; 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, -t1, -t2, -g(1) * (-t27 * t51 + t12) - g(2) * (-t24 * t51 + t10) - g(3) * t36, 0, 0, 0, 0, 0, -t2 * t22, -t2 * t25, t1, -g(1) * (t27 * t34 + t12) - g(2) * (t24 * t34 + t10) - g(3) * (t36 - t40) + (-g(3) * t50 + t9 * (pkin(2) - t21)) * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, 0, 0, 0, 0, 0, 0, 0, 0, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, g(1) * t4 - g(2) * t6 - g(3) * t43, 0, t54 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2;];
taug_reg = t11;
