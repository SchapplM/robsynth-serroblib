% Calculate inertial parameters regressor of gravitation load for
% S5RRRPP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRPP7_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP7_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP7_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP7_gravloadJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:05:54
% EndTime: 2019-12-31 21:05:55
% DurationCPUTime: 0.31s
% Computational Cost: add. (184->79), mult. (464->105), div. (0->0), fcn. (489->6), ass. (0->53)
t37 = sin(qJ(1));
t40 = cos(qJ(1));
t19 = g(1) * t40 + g(2) * t37;
t36 = sin(qJ(2));
t71 = t19 * t36;
t70 = -pkin(3) - pkin(4);
t69 = g(1) * t37;
t30 = t36 * pkin(7);
t39 = cos(qJ(2));
t32 = t39 * pkin(2);
t35 = sin(qJ(3));
t66 = t35 * t36;
t65 = t36 * t37;
t38 = cos(qJ(3));
t64 = t36 * t38;
t63 = t36 * t40;
t62 = t37 * t39;
t61 = t38 * t39;
t60 = t39 * t40;
t59 = t40 * t35;
t58 = t40 * t38;
t57 = t32 + t30;
t56 = t40 * pkin(1) + t37 * pkin(6);
t55 = qJ(4) * t35;
t54 = qJ(5) * t40;
t53 = -pkin(1) - t32;
t52 = -pkin(2) - t55;
t14 = t35 * t62 + t58;
t15 = t37 * t61 - t59;
t51 = -t14 * pkin(3) + qJ(4) * t15;
t16 = -t37 * t38 + t39 * t59;
t17 = t35 * t37 + t39 * t58;
t50 = -t16 * pkin(3) + qJ(4) * t17;
t49 = pkin(3) * t61 + t39 * t55 + t57;
t48 = pkin(2) * t60 + pkin(7) * t63 + t56;
t33 = t40 * pkin(6);
t47 = -t15 * pkin(3) - t14 * qJ(4) + t33;
t4 = g(1) * t14 - g(2) * t16;
t46 = -g(2) * t40 + t69;
t43 = t17 * pkin(3) + qJ(4) * t16 + t48;
t42 = (t53 - t30) * t69;
t2 = g(1) * t16 + g(2) * t14 + g(3) * t66;
t41 = g(1) * t17 + g(2) * t15 + g(3) * t64;
t8 = -g(3) * t39 + t71;
t25 = pkin(7) * t60;
t22 = pkin(7) * t62;
t20 = qJ(4) * t64;
t18 = g(1) * t65 - g(2) * t63;
t9 = g(3) * t36 + t19 * t39;
t7 = t8 * t38;
t6 = t8 * t35;
t5 = g(1) * t15 - g(2) * t17;
t1 = [0, 0, 0, 0, 0, 0, t46, t19, 0, 0, 0, 0, 0, 0, 0, 0, t46 * t39, -t18, -t19, -g(1) * (-t37 * pkin(1) + t33) - g(2) * t56, 0, 0, 0, 0, 0, 0, t5, -t4, t18, -g(1) * t33 - g(2) * t48 - t42, 0, 0, 0, 0, 0, 0, t5, t18, t4, -g(1) * t47 - g(2) * t43 - t42, 0, 0, 0, 0, 0, 0, t5, t4, -t18, -g(1) * (-t15 * pkin(4) + t47) - g(2) * (pkin(4) * t17 - t36 * t54 + t43) - ((-pkin(7) + qJ(5)) * t36 + t53) * t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, t9, 0, 0, 0, 0, 0, 0, 0, 0, t7, -t6, -t9, -g(1) * (-pkin(2) * t63 + t25) - g(2) * (-pkin(2) * t65 + t22) - g(3) * t57, 0, 0, 0, 0, 0, 0, t7, -t9, t6, -g(1) * t25 - g(2) * t22 - g(3) * t49 + (pkin(3) * t38 - t52) * t71, 0, 0, 0, 0, 0, 0, t7, t6, t9, -g(1) * (-t39 * t54 + t25) - g(2) * (-qJ(5) * t62 + t22) - g(3) * (pkin(4) * t61 + t49) + (g(3) * qJ(5) + t19 * (-t38 * t70 - t52)) * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t41, 0, 0, 0, 0, 0, 0, 0, 0, t2, 0, -t41, -g(1) * t50 - g(2) * t51 - g(3) * (-pkin(3) * t66 + t20), 0, 0, 0, 0, 0, 0, t2, -t41, 0, -g(1) * (-pkin(4) * t16 + t50) - g(2) * (-pkin(4) * t14 + t51) - g(3) * (t66 * t70 + t20); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8;];
taug_reg = t1;
