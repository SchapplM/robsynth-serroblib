% Calculate inertial parameters regressor of gravitation load for
% S5RRPRP8
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
% Datum: 2019-12-31 20:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPRP8_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP8_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP8_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP8_gravloadJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:04:27
% EndTime: 2019-12-31 20:04:28
% DurationCPUTime: 0.29s
% Computational Cost: add. (147->67), mult. (356->88), div. (0->0), fcn. (368->6), ass. (0->45)
t34 = sin(qJ(1));
t37 = cos(qJ(1));
t18 = g(1) * t37 + g(2) * t34;
t33 = sin(qJ(2));
t58 = t18 * t33;
t25 = t33 * qJ(3);
t36 = cos(qJ(2));
t45 = t36 * pkin(2) + t25;
t56 = g(1) * t34;
t35 = cos(qJ(4));
t32 = sin(qJ(4));
t50 = t33 * t32;
t15 = t36 * t35 + t50;
t53 = g(3) * t15;
t52 = t36 * pkin(3);
t49 = t33 * t37;
t24 = t35 * pkin(4) + pkin(3);
t48 = t36 * t24;
t47 = t36 * t32;
t46 = t36 * t37;
t44 = t37 * pkin(1) + t34 * pkin(6);
t43 = qJ(3) * t36;
t42 = pkin(4) * t50;
t41 = t32 * t46;
t40 = pkin(2) * t46 + t37 * t25 + t44;
t17 = -g(2) * t37 + t56;
t39 = -t33 * t35 + t47;
t38 = -pkin(1) - t45;
t11 = -t35 * t49 + t41;
t9 = t39 * t34;
t2 = g(1) * t11 + g(2) * t9 + t53;
t10 = t15 * t34;
t12 = t15 * t37;
t4 = g(1) * t12 + g(2) * t10 - g(3) * t39;
t31 = -qJ(5) - pkin(7);
t28 = t37 * pkin(6);
t22 = t37 * t43;
t20 = t34 * t43;
t14 = t17 * t36;
t13 = t17 * t33;
t8 = g(3) * t33 + t18 * t36;
t7 = -g(3) * t36 + t58;
t6 = g(1) * t10 - g(2) * t12;
t5 = -g(1) * t9 + g(2) * t11;
t1 = [0, 0, 0, 0, 0, 0, t17, t18, 0, 0, 0, 0, 0, 0, 0, 0, t14, -t13, -t18, -g(1) * (-t34 * pkin(1) + t28) - g(2) * t44, 0, 0, 0, 0, 0, 0, t14, -t18, t13, -g(1) * t28 - g(2) * t40 - t38 * t56, 0, 0, 0, 0, 0, 0, t6, t5, t18, -g(1) * (-t37 * pkin(7) + t28) - g(2) * (pkin(3) * t46 + t40) + (-g(1) * (t38 - t52) + g(2) * pkin(7)) * t34, 0, 0, 0, 0, 0, 0, t6, t5, t18, -g(1) * (t37 * t31 + t28) - g(2) * (t24 * t46 + t37 * t42 + t40) + (-g(1) * (t38 - t42 - t48) - g(2) * t31) * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, 0, 0, 0, t7, 0, -t8, -g(1) * (-pkin(2) * t49 + t22) - g(2) * (-t34 * t33 * pkin(2) + t20) - g(3) * t45, 0, 0, 0, 0, 0, 0, -t2, -t4, 0, -g(1) * t22 - g(2) * t20 - g(3) * (t45 + t52) + (pkin(2) + pkin(3)) * t58, 0, 0, 0, 0, 0, 0, -t2, -t4, 0, -g(1) * (pkin(4) * t41 + t22) - g(2) * (t34 * pkin(4) * t47 + t20) - g(3) * (t45 + t48) + (-g(3) * pkin(4) * t32 + t18 * (pkin(2) + t24)) * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t4, 0, 0, 0, 0, 0, 0, 0, 0, t2, t4, 0, (t18 * t39 + t53) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17;];
taug_reg = t1;
