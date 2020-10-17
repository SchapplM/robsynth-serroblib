% Calculate minimal parameter regressor of gravitation load for
% S6RRPPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3]';
% 
% Output:
% taug_reg [6x29]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPPPR5_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR5_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR5_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR5_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 08:54:16
% EndTime: 2019-05-06 08:54:18
% DurationCPUTime: 0.34s
% Computational Cost: add. (183->80), mult. (480->110), div. (0->0), fcn. (525->8), ass. (0->55)
t38 = cos(qJ(2));
t35 = sin(qJ(2));
t36 = sin(qJ(1));
t39 = cos(qJ(1));
t50 = g(1) * t39 + g(2) * t36;
t71 = t50 * t35;
t8 = -g(3) * t38 + t71;
t70 = g(1) * t36;
t67 = g(3) * t35;
t28 = t38 * pkin(2);
t33 = cos(pkin(9));
t65 = t33 * t38;
t64 = t35 * t36;
t63 = t35 * t39;
t62 = t36 * t38;
t61 = t38 * t39;
t32 = sin(pkin(9));
t60 = t39 * t32;
t59 = t39 * t33;
t26 = t35 * qJ(3);
t58 = t28 + t26;
t57 = qJ(3) * t38;
t56 = qJ(4) * t32;
t55 = -pkin(1) - t28;
t54 = -pkin(2) - t56;
t53 = pkin(3) * t65 + t38 * t56 + t58;
t52 = pkin(2) * t61 + t36 * pkin(7) + (pkin(1) + t26) * t39;
t12 = t32 * t62 + t59;
t13 = t33 * t62 - t60;
t29 = t39 * pkin(7);
t51 = -t13 * pkin(3) - t12 * qJ(4) + t29;
t14 = -t36 * t33 + t38 * t60;
t4 = g(1) * t12 - g(2) * t14;
t15 = t36 * t32 + t38 * t59;
t5 = g(1) * t13 - g(2) * t15;
t49 = -g(2) * t39 + t70;
t34 = sin(qJ(6));
t37 = cos(qJ(6));
t48 = t12 * t37 + t13 * t34;
t47 = t12 * t34 - t13 * t37;
t46 = t32 * t37 + t33 * t34;
t45 = t32 * t34 - t33 * t37;
t43 = g(3) * t45;
t42 = t15 * pkin(3) + t14 * qJ(4) + t52;
t40 = (t55 - t26) * t70;
t21 = t39 * t57;
t18 = t36 * t57;
t16 = g(1) * t64 - g(2) * t63;
t9 = t50 * t38 + t67;
t7 = t8 * t33;
t6 = t8 * t32;
t3 = t14 * t37 + t15 * t34;
t2 = -t14 * t34 + t15 * t37;
t1 = -g(1) * t14 - g(2) * t12 - t32 * t67;
t10 = [0, t49, t50, 0, 0, 0, 0, 0, t49 * t38, -t16, t5, -t4, t16, -g(1) * t29 - g(2) * t52 - t40, t16, -t5, t4, -g(1) * t51 - g(2) * t42 - t40, t4, -t16, t5, -g(1) * (-t13 * qJ(5) + t51) - g(2) * (pkin(4) * t63 + t15 * qJ(5) + t42) - ((-pkin(4) - qJ(3)) * t35 + t55) * t70, 0, 0, 0, 0, 0, g(1) * t48 - g(2) * t3, -g(1) * t47 - g(2) * t2; 0, 0, 0, 0, 0, 0, 0, 0, t8, t9, t7, -t6, -t9, -g(1) * (-pkin(2) * t63 + t21) - g(2) * (-pkin(2) * t64 + t18) - g(3) * t58, -t9, -t7, t6, -g(1) * t21 - g(2) * t18 - g(3) * t53 + (pkin(3) * t33 - t54) * t71, t6, t9, t7, -g(1) * (pkin(4) * t61 + t21) - g(2) * (pkin(4) * t62 + t18) - g(3) * (qJ(5) * t65 + t53) + (-g(3) * pkin(4) + t50 * (-(-pkin(3) - qJ(5)) * t33 - t54)) * t35, 0, 0, 0, 0, 0, t8 * t46, t38 * t43 - t45 * t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, 0, 0, 0, -t8, 0, 0, 0, -t8, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, 0, 0, t1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t15 - g(2) * t13 - t33 * t67, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t2 + g(2) * t47 + t35 * t43, g(1) * t3 + g(2) * t48 + t46 * t67;];
taug_reg  = t10;
