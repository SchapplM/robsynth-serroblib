% Calculate minimal parameter regressor of gravitation load for
% S6RRPPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
% 
% Output:
% taug_reg [6x32]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPPRR7_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR7_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR7_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR7_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 11:07:01
% EndTime: 2019-05-06 11:07:02
% DurationCPUTime: 0.33s
% Computational Cost: add. (213->79), mult. (555->132), div. (0->0), fcn. (675->10), ass. (0->47)
t34 = sin(qJ(2));
t35 = sin(qJ(1));
t38 = cos(qJ(2));
t39 = cos(qJ(1));
t50 = cos(pkin(6));
t48 = t39 * t50;
t20 = t34 * t48 + t35 * t38;
t32 = sin(qJ(6));
t36 = cos(qJ(6));
t19 = t35 * t34 - t38 * t48;
t33 = sin(qJ(5));
t37 = cos(qJ(5));
t31 = sin(pkin(6));
t55 = t31 * t39;
t9 = t19 * t37 + t33 * t55;
t63 = -t20 * t36 + t9 * t32;
t62 = t20 * t32 + t9 * t36;
t61 = g(3) * t31;
t58 = t31 * t34;
t57 = t31 * t35;
t56 = t31 * t38;
t54 = t32 * t37;
t53 = t36 * t37;
t52 = pkin(2) * t56 + qJ(3) * t58;
t51 = qJ(4) * t31;
t49 = t35 * t50;
t47 = -t19 * pkin(2) + t20 * qJ(3);
t21 = t39 * t34 + t38 * t49;
t22 = -t34 * t49 + t39 * t38;
t46 = -t21 * pkin(2) + t22 * qJ(3);
t6 = g(1) * t19 - g(2) * t21;
t7 = g(1) * t20 - g(2) * t22;
t45 = g(1) * t39 + g(2) * t35;
t44 = g(1) * t35 - g(2) * t39;
t43 = t39 * pkin(1) + t22 * pkin(2) + pkin(8) * t57 + t21 * qJ(3);
t8 = -t19 * t33 + t37 * t55;
t11 = -t21 * t33 - t37 * t57;
t42 = g(1) * t11 + g(2) * t8 + g(3) * (t33 * t56 - t50 * t37);
t41 = -t35 * pkin(1) - t20 * pkin(2) + pkin(8) * t55 - t19 * qJ(3);
t3 = -g(1) * t21 - g(2) * t19 + g(3) * t56;
t40 = g(1) * t22 + g(2) * t20 + g(3) * t58;
t23 = t45 * t31;
t18 = t50 * t33 + t37 * t56;
t12 = t21 * t37 - t33 * t57;
t2 = t12 * t36 + t22 * t32;
t1 = -t12 * t32 + t22 * t36;
t4 = [0, t44, t45, 0, 0, 0, 0, 0, t7, -t6, t7, -t23, t6, -g(1) * t41 - g(2) * t43, t6, -t7, t23, -g(1) * (-t20 * pkin(3) - t39 * t51 + t41) - g(2) * (t22 * pkin(3) - t35 * t51 + t43) 0, 0, 0, 0, 0, g(1) * t9 - g(2) * t12, g(1) * t8 - g(2) * t11, 0, 0, 0, 0, 0, g(1) * t62 - g(2) * t2, -g(1) * t63 - g(2) * t1; 0, 0, 0, 0, 0, 0, 0, 0, -t3, t40, -t3, 0, -t40, -g(1) * t46 - g(2) * t47 - g(3) * t52, -t40, t3, 0, -g(1) * (-t21 * pkin(3) + t46) - g(2) * (-t19 * pkin(3) + t47) - g(3) * (pkin(3) * t56 + t52) 0, 0, 0, 0, 0, -t40 * t37, t40 * t33, 0, 0, 0, 0, 0, -g(1) * (-t21 * t32 + t22 * t53) - g(2) * (-t19 * t32 + t20 * t53) - (t32 * t38 + t34 * t53) * t61, -g(1) * (-t21 * t36 - t22 * t54) - g(2) * (-t19 * t36 - t20 * t54) - (-t34 * t54 + t36 * t38) * t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, 0, 0, 0, t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3) * t50 + t44 * t31, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t42, g(1) * t12 + g(2) * t9 - g(3) * t18, 0, 0, 0, 0, 0, -t42 * t36, t42 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t1 + g(2) * t63 - g(3) * (t18 * t32 + t36 * t58) g(1) * t2 + g(2) * t62 - g(3) * (t18 * t36 - t32 * t58);];
taug_reg  = t4;
