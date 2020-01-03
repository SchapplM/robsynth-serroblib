% Calculate inertial parameters regressor of gravitation load for
% S5RRPRP11
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
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:14
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPRP11_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP11_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP11_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP11_gravloadJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t31 = sin(qJ(1));
t34 = cos(qJ(1));
t13 = g(1) * t34 + g(2) * t31;
t30 = sin(qJ(2));
t64 = t13 * t30;
t33 = cos(qJ(2));
t63 = -g(3) * t30 - t13 * t33;
t45 = qJ(3) * t33;
t15 = t31 * t45;
t17 = t34 * t45;
t20 = t30 * qJ(3);
t24 = t33 * pkin(2);
t48 = t24 + t20;
t60 = pkin(2) + pkin(7);
t61 = t60 * t64 - g(1) * t17 - g(2) * t15 - g(3) * (t33 * pkin(7) + t48);
t59 = pkin(2) * t30;
t58 = g(1) * t31;
t54 = g(3) * t33;
t29 = sin(qJ(4));
t53 = t31 * t29;
t32 = cos(qJ(4));
t52 = t31 * t32;
t51 = t33 * t34;
t50 = t34 * t29;
t49 = t34 * t32;
t25 = t34 * pkin(6);
t47 = t34 * pkin(3) + t25;
t46 = t34 * pkin(1) + t31 * pkin(6);
t44 = -pkin(1) - t20;
t43 = pkin(2) * t51 + t34 * t20 + t46;
t7 = -t30 * t49 + t53;
t9 = t30 * t52 + t50;
t42 = g(1) * t9 + g(2) * t7;
t41 = -g(2) * t34 + t58;
t39 = t31 * pkin(3) + pkin(7) * t51 + t43;
t1 = g(1) * t7 - g(2) * t9 + t32 * t54;
t10 = -t30 * t53 + t49;
t8 = t30 * t50 + t52;
t38 = -g(1) * t8 + g(2) * t10 + t29 * t54;
t36 = (-t60 * t33 + t44) * t58;
t12 = t41 * t33;
t11 = t41 * t30;
t5 = -t54 + t64;
t4 = t63 * t32;
t3 = t63 * t29;
t2 = -g(1) * t10 - g(2) * t8;
t6 = [0, 0, 0, 0, 0, 0, t41, t13, 0, 0, 0, 0, 0, 0, 0, 0, t12, -t11, -t13, -g(1) * (-t31 * pkin(1) + t25) - g(2) * t46, 0, 0, 0, 0, 0, 0, -t13, -t12, t11, -g(1) * t25 - g(2) * t43 - (t44 - t24) * t58, 0, 0, 0, 0, 0, 0, t2, t42, t12, -g(1) * t47 - g(2) * t39 - t36, 0, 0, 0, 0, 0, 0, t2, t12, -t42, -g(1) * (t10 * pkin(4) + t9 * qJ(5) + t47) - g(2) * (t8 * pkin(4) + t7 * qJ(5) + t39) - t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, -t63, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, t63, -g(1) * (-t34 * t59 + t17) - g(2) * (-t31 * t59 + t15) - g(3) * t48, 0, 0, 0, 0, 0, 0, t3, t4, t5, t61, 0, 0, 0, 0, 0, 0, t3, t5, -t4, t63 * (pkin(4) * t29 - qJ(5) * t32) + t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, -t38, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, t38, -g(1) * (-t7 * pkin(4) + t8 * qJ(5)) - g(2) * (t9 * pkin(4) - t10 * qJ(5)) - (-pkin(4) * t32 - qJ(5) * t29) * t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t6;
