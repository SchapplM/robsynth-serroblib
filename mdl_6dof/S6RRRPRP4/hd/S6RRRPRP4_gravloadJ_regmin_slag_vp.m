% Calculate minimal parameter regressor of gravitation load for
% S6RRRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
% 
% Output:
% taug_reg [6x32]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRPRP4_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP4_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP4_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP4_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 07:47:27
% EndTime: 2019-05-07 07:47:28
% DurationCPUTime: 0.31s
% Computational Cost: add. (328->89), mult. (431->106), div. (0->0), fcn. (438->8), ass. (0->55)
t31 = qJ(2) + qJ(3);
t28 = sin(t31);
t23 = t28 * qJ(4);
t29 = cos(t31);
t51 = t29 * pkin(3) + t23;
t34 = sin(qJ(1));
t37 = cos(qJ(1));
t16 = g(1) * t37 + g(2) * t34;
t6 = g(3) * t28 + t16 * t29;
t67 = pkin(3) + pkin(9);
t33 = sin(qJ(2));
t66 = pkin(2) * t33;
t65 = pkin(3) * t28;
t61 = g(3) * t29;
t24 = t29 * pkin(9);
t32 = sin(qJ(5));
t60 = t29 * t32;
t59 = t29 * t37;
t58 = t34 * t32;
t35 = cos(qJ(5));
t57 = t34 * t35;
t56 = t37 * t32;
t55 = t37 * t35;
t38 = -pkin(8) - pkin(7);
t54 = t37 * t38;
t50 = qJ(4) * t29;
t17 = t34 * t50;
t48 = pkin(5) * t60;
t53 = t34 * t48 + t17;
t19 = t37 * t50;
t52 = t37 * t48 + t19;
t49 = qJ(6) * t35;
t36 = cos(qJ(2));
t30 = t36 * pkin(2);
t27 = t30 + pkin(1);
t47 = pkin(3) * t59 + (t23 + t27) * t37;
t46 = t29 * t49;
t45 = t28 * t32 * pkin(5) + t24 + t51;
t11 = t28 * t57 + t56;
t9 = -t28 * t55 + t58;
t44 = g(1) * t11 + g(2) * t9;
t43 = -t65 - t66;
t42 = g(1) * t34 - g(2) * t37;
t41 = -t27 - t51;
t1 = g(1) * t9 - g(2) * t11 + t35 * t61;
t10 = t28 * t56 + t57;
t12 = -t28 * t58 + t55;
t40 = -g(1) * t10 + g(2) * t12 + g(3) * t60;
t8 = t42 * t29;
t7 = t42 * t28;
t5 = t16 * t28 - t61;
t4 = t6 * t35;
t3 = t6 * t32;
t2 = -g(1) * t12 - g(2) * t10;
t13 = [0, t42, t16, 0, 0, 0, 0, 0, t42 * t36, -t42 * t33, 0, 0, 0, 0, 0, t8, -t7, -t16, -t8, t7, g(1) * t54 - g(2) * t47 + (-g(1) * t41 + g(2) * t38) * t34, 0, 0, 0, 0, 0, t2, t44, t2, t8, -t44, -g(1) * (t37 * pkin(4) + t12 * pkin(5) + t11 * qJ(6) - t54) - g(2) * (t10 * pkin(5) + pkin(9) * t59 + t9 * qJ(6) + t47) + (-g(1) * (t41 - t24) - g(2) * (pkin(4) - t38)) * t34; 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t36 + t16 * t33, g(3) * t33 + t16 * t36, 0, 0, 0, 0, 0, t5, t6, 0, -t5, -t6, -g(1) * (t43 * t37 + t19) - g(2) * (t43 * t34 + t17) - g(3) * (t30 + t51) 0, 0, 0, 0, 0, -t3, -t4, -t3, t5, t4, -g(1) * t52 - g(2) * t53 - g(3) * (-t28 * t49 + t30 + t45) + t16 * (t67 * t28 + t46 + t66); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, -t5, -t6, -g(1) * (-t37 * t65 + t19) - g(2) * (-t34 * t65 + t17) - g(3) * t51, 0, 0, 0, 0, 0, -t3, -t4, -t3, t5, t4, -g(1) * (-t37 * t46 + t52) - g(2) * (-t34 * t46 + t53) - g(3) * t45 + (g(3) * t49 + t16 * t67) * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, -t40, t1, 0, t40, -g(1) * (-t9 * pkin(5) + t10 * qJ(6)) - g(2) * (t11 * pkin(5) - t12 * qJ(6)) - (-pkin(5) * t35 - qJ(6) * t32) * t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t13;
