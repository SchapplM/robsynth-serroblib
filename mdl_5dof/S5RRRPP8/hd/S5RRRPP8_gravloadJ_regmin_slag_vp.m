% Calculate minimal parameter regressor of gravitation load for
% S5RRRPP8
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
% taug_reg [5x25]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRPP8_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP8_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP8_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP8_gravloadJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t33 = sin(qJ(2));
t34 = sin(qJ(1));
t37 = cos(qJ(1));
t44 = g(1) * t37 + g(2) * t34;
t64 = t44 * t33;
t63 = g(1) * t34;
t27 = t33 * pkin(7);
t36 = cos(qJ(2));
t29 = t36 * pkin(2);
t32 = sin(qJ(3));
t60 = t32 * t33;
t35 = cos(qJ(3));
t59 = t33 * t35;
t58 = t33 * t37;
t57 = t34 * t36;
t56 = t35 * t36;
t55 = t36 * t37;
t54 = t37 * t32;
t53 = t37 * t35;
t52 = -pkin(3) - qJ(5);
t51 = qJ(4) * t32;
t50 = -pkin(1) - t29;
t49 = -pkin(2) - t51;
t13 = t32 * t57 + t53;
t14 = t34 * t56 - t54;
t48 = -t13 * pkin(3) + t14 * qJ(4);
t15 = -t34 * t35 + t36 * t54;
t16 = t34 * t32 + t36 * t53;
t47 = -t15 * pkin(3) + t16 * qJ(4);
t46 = pkin(3) * t56 + t36 * t51 + t27 + t29;
t45 = -t14 * pkin(3) + t37 * pkin(6) - t13 * qJ(4);
t4 = g(1) * t13 - g(2) * t15;
t5 = g(1) * t14 - g(2) * t16;
t43 = -g(2) * t37 + t63;
t41 = t37 * pkin(1) + pkin(2) * t55 + t16 * pkin(3) + t34 * pkin(6) + pkin(7) * t58 + t15 * qJ(4);
t2 = g(1) * t15 + g(2) * t13 + g(3) * t60;
t39 = g(1) * t16 + g(2) * t14 + g(3) * t59;
t38 = -g(3) * t36 + t64;
t23 = pkin(7) * t55;
t20 = pkin(7) * t57;
t18 = qJ(4) * t59;
t17 = t43 * t33;
t8 = g(3) * t33 + t44 * t36;
t7 = t38 * t35;
t6 = t38 * t32;
t1 = [0, t43, t44, 0, 0, 0, 0, 0, t43 * t36, -t17, 0, 0, 0, 0, 0, t5, -t4, t17, -t5, t4, -g(1) * t45 - g(2) * t41 - (t50 - t27) * t63, t17, t4, t5, -g(1) * (-t14 * qJ(5) + t45) - g(2) * (pkin(4) * t58 + t16 * qJ(5) + t41) - ((-pkin(4) - pkin(7)) * t33 + t50) * t63; 0, 0, 0, 0, 0, 0, 0, 0, t38, t8, 0, 0, 0, 0, 0, t7, -t6, -t8, -t7, t6, -g(1) * t23 - g(2) * t20 - g(3) * t46 + (pkin(3) * t35 - t49) * t64, -t8, t6, t7, -g(1) * (pkin(4) * t55 + t23) - g(2) * (pkin(4) * t57 + t20) - g(3) * (qJ(5) * t56 + t46) + (-g(3) * pkin(4) + t44 * (-t52 * t35 - t49)) * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t39, 0, -t2, -t39, -g(1) * t47 - g(2) * t48 - g(3) * (-pkin(3) * t60 + t18), 0, -t39, t2, -g(1) * (-t15 * qJ(5) + t47) - g(2) * (-t13 * qJ(5) + t48) - g(3) * (t52 * t60 + t18); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, 0, 0, 0, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t39;];
taug_reg = t1;
