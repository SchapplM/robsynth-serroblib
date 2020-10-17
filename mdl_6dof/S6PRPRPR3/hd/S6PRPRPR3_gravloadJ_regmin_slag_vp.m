% Calculate minimal parameter regressor of gravitation load for
% S6PRPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
% 
% Output:
% taug_reg [6x23]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRPRPR3_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR3_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR3_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR3_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 22:36:19
% EndTime: 2019-05-04 22:36:20
% DurationCPUTime: 0.34s
% Computational Cost: add. (299->76), mult. (791->131), div. (0->0), fcn. (1013->12), ass. (0->47)
t30 = sin(pkin(11));
t37 = sin(qJ(2));
t40 = cos(qJ(2));
t51 = cos(pkin(11));
t46 = -t37 * t30 + t40 * t51;
t31 = sin(pkin(10));
t33 = cos(pkin(10));
t34 = cos(pkin(6));
t55 = t34 * t40;
t66 = -t31 * t55 - t33 * t37;
t23 = -t40 * t30 - t37 * t51;
t42 = t46 * t34;
t12 = t33 * t23 - t31 * t42;
t32 = sin(pkin(6));
t19 = t46 * t32;
t9 = t31 * t23 + t33 * t42;
t65 = -g(1) * t12 - g(2) * t9 - g(3) * t19;
t61 = t31 * t37;
t36 = sin(qJ(4));
t60 = t32 * t36;
t39 = cos(qJ(4));
t59 = t32 * t39;
t58 = t32 * t40;
t56 = t34 * t37;
t35 = sin(qJ(6));
t54 = t35 * t36;
t38 = cos(qJ(6));
t53 = t36 * t38;
t49 = t33 * t55;
t21 = t23 * t34;
t10 = -t33 * t21 + t31 * t46;
t11 = -t31 * t21 - t33 * t46;
t20 = t23 * t32;
t14 = -t20 * t36 - t34 * t39;
t4 = t10 * t36 + t33 * t59;
t6 = -t11 * t36 - t31 * t59;
t45 = g(1) * t6 + g(2) * t4 + g(3) * t14;
t15 = -t20 * t39 + t34 * t36;
t5 = t10 * t39 - t33 * t60;
t7 = -t11 * t39 + t31 * t60;
t44 = g(1) * t7 + g(2) * t5 + g(3) * t15;
t41 = -g(1) * t66 - g(3) * t58;
t24 = pkin(2) * t49;
t18 = -g(3) * t34 + (-g(1) * t31 + g(2) * t33) * t32;
t3 = t65 * t39;
t2 = t65 * t36;
t1 = [-g(3), 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, -g(2) * (t49 - t61) + t41, -g(1) * (t31 * t56 - t33 * t40) - g(2) * (-t31 * t40 - t33 * t56) + g(3) * t32 * t37, -g(2) * t24 + (g(2) * t61 + t41) * pkin(2), 0, 0, 0, 0, 0, t3, -t2, g(1) * t11 - g(2) * t10 + g(3) * t20, -t3, t2, -g(1) * (t66 * pkin(2) - t11 * pkin(8)) - g(2) * (-pkin(2) * t61 + pkin(8) * t10 + t24) - g(3) * (pkin(2) * t58 - t20 * pkin(8)) + t65 * (pkin(4) * t39 + qJ(5) * t36 + pkin(3)) 0, 0, 0, 0, 0, -g(1) * (-t11 * t38 + t12 * t54) - g(2) * (t10 * t38 + t9 * t54) - g(3) * (t19 * t54 - t20 * t38) -g(1) * (t11 * t35 + t12 * t53) - g(2) * (-t10 * t35 + t9 * t53) - g(3) * (t19 * t53 + t20 * t35); 0, 0, 0, 0, t18, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, t44, 0, -t45, -t44, -g(1) * (-t6 * pkin(4) + t7 * qJ(5)) - g(2) * (-t4 * pkin(4) + t5 * qJ(5)) - g(3) * (-t14 * pkin(4) + t15 * qJ(5)) 0, 0, 0, 0, 0, -t44 * t35, -t44 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t45, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t12 * t35 + t6 * t38) - g(2) * (t9 * t35 + t4 * t38) - g(3) * (t14 * t38 + t19 * t35) -g(1) * (t12 * t38 - t6 * t35) - g(2) * (-t4 * t35 + t9 * t38) - g(3) * (-t14 * t35 + t19 * t38);];
taug_reg  = t1;
