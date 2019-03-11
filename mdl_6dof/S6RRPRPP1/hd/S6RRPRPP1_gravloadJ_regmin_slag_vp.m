% Calculate minimal parameter regressor of gravitation load for
% S6RRPRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3,theta5]';
% 
% Output:
% taug_reg [6x25]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRPP1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP1_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPP1_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t38 = cos(qJ(4));
t23 = t38 * pkin(4) + pkin(3);
t32 = qJ(2) + pkin(9);
t26 = sin(t32);
t28 = cos(t32);
t33 = -qJ(5) - pkin(8);
t68 = t28 * t23 - t26 * t33;
t37 = sin(qJ(1));
t40 = cos(qJ(1));
t16 = g(1) * t40 + g(2) * t37;
t44 = -g(3) * t28 + t16 * t26;
t64 = g(3) * t26;
t67 = t16 * t28 + t64;
t31 = qJ(4) + pkin(10);
t25 = sin(t31);
t61 = t37 * t25;
t27 = cos(t31);
t60 = t37 * t27;
t35 = sin(qJ(4));
t59 = t37 * t35;
t58 = t37 * t38;
t57 = t40 * t25;
t56 = t40 * t27;
t34 = -qJ(3) - pkin(7);
t55 = t40 * t34;
t54 = t40 * t35;
t53 = t40 * t38;
t52 = t28 * t54;
t39 = cos(qJ(2));
t29 = t39 * pkin(2);
t24 = t29 + pkin(1);
t51 = t40 * t24 - t37 * t34;
t50 = t29 + t68;
t15 = g(1) * t37 - g(2) * t40;
t36 = sin(qJ(2));
t49 = -pkin(2) * t36 - t28 * t33;
t48 = pkin(5) * t27 + qJ(6) * t25;
t8 = t28 * t59 + t53;
t3 = t28 * t61 + t56;
t5 = t28 * t57 - t60;
t46 = g(1) * t5 + g(2) * t3 + t25 * t64;
t45 = pkin(4) * t59 + t68 * t40 + t51;
t43 = -g(3) * t39 + t16 * t36;
t42 = pkin(4) * t54 - t55 + (-t24 - t68) * t37;
t21 = pkin(4) * t58;
t11 = t28 * t53 + t59;
t10 = -t52 + t58;
t9 = -t28 * t58 + t54;
t7 = t15 * t26;
t6 = t28 * t56 + t61;
t4 = t28 * t60 - t57;
t1 = [0, t15, t16, 0, 0, 0, 0, 0, t15 * t39, -t15 * t36, -t16, -g(1) * (-t37 * t24 - t55) - g(2) * t51, 0, 0, 0, 0, 0, -g(1) * t9 - g(2) * t11, -g(1) * t8 - g(2) * t10, t7, -g(1) * t42 - g(2) * t45, g(1) * t4 - g(2) * t6, t7, g(1) * t3 - g(2) * t5, -g(1) * (-t4 * pkin(5) - t3 * qJ(6) + t42) - g(2) * (t6 * pkin(5) + t5 * qJ(6) + t45); 0, 0, 0, 0, 0, 0, 0, 0, t43, g(3) * t36 + t16 * t39, 0, t43 * pkin(2), 0, 0, 0, 0, 0, t44 * t38, -t44 * t35, -t67, -g(3) * t50 + t16 * (t23 * t26 - t49) t44 * t27, -t67, t44 * t25, -g(3) * (t48 * t28 + t50) + t16 * (-(-t23 - t48) * t26 - t49); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, 0, 0, 0, 0, 0, 0, 0, 0, -t15, 0, 0, 0, -t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t10 + g(2) * t8 + t35 * t64, g(1) * t11 - g(2) * t9 + t38 * t64, 0, -g(1) * t21 + (g(2) * t53 + t67 * t35) * pkin(4), t46, 0, -g(1) * t6 - g(2) * t4 - t27 * t64, -g(1) * (-pkin(4) * t52 - t5 * pkin(5) + t6 * qJ(6) + t21) - g(2) * (-t8 * pkin(4) - t3 * pkin(5) + t4 * qJ(6)) - (-pkin(4) * t35 - pkin(5) * t25 + qJ(6) * t27) * t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44, 0, 0, 0, -t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46;];
taug_reg  = t1;
