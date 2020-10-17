% Calculate inertial parameters regressor of gravitation load for
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
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPRP10_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP10_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP10_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP10_gravloadJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:11:05
% EndTime: 2019-12-31 20:11:07
% DurationCPUTime: 0.29s
% Computational Cost: add. (135->70), mult. (326->88), div. (0->0), fcn. (322->6), ass. (0->49)
t30 = sin(qJ(1));
t33 = cos(qJ(1));
t15 = g(1) * t33 + g(2) * t30;
t29 = sin(qJ(2));
t63 = t15 * t29;
t21 = t29 * qJ(3);
t32 = cos(qJ(2));
t43 = t32 * pkin(2) + t21;
t28 = sin(qJ(4));
t45 = t33 * t28;
t31 = cos(qJ(4));
t48 = t30 * t31;
t11 = t29 * t48 + t45;
t53 = g(3) * t32;
t44 = t33 * t31;
t49 = t30 * t28;
t9 = t29 * t44 - t49;
t1 = -g(1) * t9 - g(2) * t11 + t31 * t53;
t8 = g(3) * t29 + t15 * t32;
t60 = pkin(2) * t29;
t59 = pkin(4) * t28;
t58 = g(1) * t30;
t52 = t32 * pkin(7);
t50 = t28 * t32;
t27 = -qJ(5) - pkin(7);
t47 = t32 * t27;
t46 = t32 * t33;
t42 = t33 * pkin(1) + t30 * pkin(6);
t41 = qJ(3) * t32;
t40 = pkin(4) * t50;
t38 = t29 * t45;
t37 = pkin(2) * t46 + t33 * t21 + t42;
t36 = -g(2) * t33 + t58;
t35 = -pkin(1) - t43;
t24 = t33 * pkin(6);
t20 = t31 * pkin(4) + pkin(3);
t18 = t33 * t41;
t16 = t30 * t41;
t14 = t36 * t32;
t13 = t36 * t29;
t12 = -t29 * t49 + t44;
t10 = t38 + t48;
t7 = -t53 + t63;
t6 = t8 * t31;
t5 = t8 * t28;
t4 = -g(1) * t12 - g(2) * t10;
t3 = g(1) * t11 - g(2) * t9;
t2 = g(1) * t10 - g(2) * t12 - g(3) * t50;
t17 = [0, 0, 0, 0, 0, 0, t36, t15, 0, 0, 0, 0, 0, 0, 0, 0, t14, -t13, -t15, -g(1) * (-t30 * pkin(1) + t24) - g(2) * t42, 0, 0, 0, 0, 0, 0, -t15, -t14, t13, -g(1) * t24 - g(2) * t37 - t35 * t58, 0, 0, 0, 0, 0, 0, t4, t3, t14, -g(1) * (t33 * pkin(3) + t24) - g(2) * (pkin(7) * t46 + t37) + (-g(1) * (t35 - t52) - g(2) * pkin(3)) * t30, 0, 0, 0, 0, 0, 0, t4, t3, t14, -g(1) * (t33 * t20 + t24) - g(2) * (pkin(4) * t38 - t27 * t46 + t37) + (-g(1) * (-t29 * t59 + t35 + t47) - g(2) * t20) * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, -t8, -g(1) * (-t33 * t60 + t18) - g(2) * (-t30 * t60 + t16) - g(3) * t43, 0, 0, 0, 0, 0, 0, -t5, -t6, t7, -g(1) * t18 - g(2) * t16 - g(3) * (t43 + t52) + (pkin(2) + pkin(7)) * t63, 0, 0, 0, 0, 0, 0, -t5, -t6, t7, -g(1) * (t33 * t40 + t18) - g(2) * (t30 * t40 + t16) - g(3) * (t43 - t47) + (-g(3) * t59 + t15 * (pkin(2) - t27)) * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8;];
taug_reg = t17;
