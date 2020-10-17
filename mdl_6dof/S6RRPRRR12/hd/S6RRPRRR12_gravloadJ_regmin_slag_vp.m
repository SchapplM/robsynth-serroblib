% Calculate minimal parameter regressor of gravitation load for
% S6RRPRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6]';
% 
% Output:
% taug_reg [6x35]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 14:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRRR12_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR12_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR12_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR12_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 00:44:27
% EndTime: 2019-05-07 00:44:28
% DurationCPUTime: 0.36s
% Computational Cost: add. (309->84), mult. (599->142), div. (0->0), fcn. (744->12), ass. (0->50)
t31 = sin(qJ(2));
t32 = sin(qJ(1));
t35 = cos(qJ(2));
t36 = cos(qJ(1));
t47 = cos(pkin(6));
t44 = t36 * t47;
t19 = t31 * t44 + t32 * t35;
t29 = sin(qJ(6));
t33 = cos(qJ(6));
t18 = t32 * t31 - t35 * t44;
t27 = qJ(4) + qJ(5);
t25 = sin(t27);
t26 = cos(t27);
t28 = sin(pkin(6));
t50 = t28 * t36;
t40 = -t18 * t25 + t26 * t50;
t59 = t19 * t33 + t29 * t40;
t58 = -t19 * t29 + t33 * t40;
t57 = g(3) * t28;
t54 = t25 * t29;
t53 = t25 * t33;
t52 = t28 * t32;
t51 = t28 * t35;
t49 = t29 * t31;
t48 = t31 * t33;
t30 = sin(qJ(4));
t34 = cos(qJ(4));
t46 = t18 * t30 - t34 * t50;
t45 = t32 * t47;
t20 = t36 * t31 + t35 * t45;
t43 = g(1) * t18 - g(2) * t20;
t21 = -t31 * t45 + t36 * t35;
t42 = g(1) * t19 - g(2) * t21;
t41 = g(1) * t36 + g(2) * t32;
t10 = t18 * t26 + t25 * t50;
t39 = t18 * t34 + t30 * t50;
t8 = t20 * t26 - t25 * t52;
t38 = g(1) * t8 + g(2) * t10 + g(3) * (-t47 * t25 - t26 * t51);
t7 = -g(1) * t20 - g(2) * t18 + g(3) * t51;
t37 = g(1) * t21 + g(2) * t19 + t31 * t57;
t17 = -t25 * t51 + t47 * t26;
t14 = t20 * t30 + t34 * t52;
t13 = t20 * t34 - t30 * t52;
t9 = t20 * t25 + t26 * t52;
t6 = t21 * t29 + t9 * t33;
t5 = t21 * t33 - t9 * t29;
t4 = g(1) * t9 - g(2) * t40 + g(3) * t17;
t2 = t38 * t33;
t1 = t38 * t29;
t3 = [0, g(1) * t32 - g(2) * t36, t41, 0, 0, 0, 0, 0, t42, -t43, -t41 * t28, -t42, t43, -g(1) * (-t32 * pkin(1) - t19 * pkin(2) + pkin(8) * t50 - t18 * qJ(3)) - g(2) * (t36 * pkin(1) + t21 * pkin(2) + pkin(8) * t52 + t20 * qJ(3)) 0, 0, 0, 0, 0, g(1) * t46 - g(2) * t14, g(1) * t39 - g(2) * t13, 0, 0, 0, 0, 0, -g(1) * t40 - g(2) * t9, g(1) * t10 - g(2) * t8, 0, 0, 0, 0, 0, -g(1) * t58 - g(2) * t6, g(1) * t59 - g(2) * t5; 0, 0, 0, 0, 0, 0, 0, 0, -t7, t37, 0, t7, -t37, -g(1) * (-t20 * pkin(2) + t21 * qJ(3)) - g(2) * (-t18 * pkin(2) + t19 * qJ(3)) - (pkin(2) * t35 + qJ(3) * t31) * t57, 0, 0, 0, 0, 0, -t37 * t30, -t37 * t34, 0, 0, 0, 0, 0, -t37 * t25, -t37 * t26, 0, 0, 0, 0, 0, -g(1) * (-t20 * t29 + t21 * t53) - g(2) * (-t18 * t29 + t19 * t53) - (t25 * t48 + t29 * t35) * t57, -g(1) * (-t20 * t33 - t21 * t54) - g(2) * (-t18 * t33 - t19 * t54) - (-t25 * t49 + t33 * t35) * t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t13 - g(2) * t39 - g(3) * (-t47 * t30 - t34 * t51) g(1) * t14 + g(2) * t46 - g(3) * (t30 * t51 - t47 * t34) 0, 0, 0, 0, 0, -t38, t4, 0, 0, 0, 0, 0, -t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38, t4, 0, 0, 0, 0, 0, -t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t5 - g(2) * t59 - g(3) * (-t17 * t29 + t28 * t48) g(1) * t6 - g(2) * t58 - g(3) * (-t17 * t33 - t28 * t49);];
taug_reg  = t3;
