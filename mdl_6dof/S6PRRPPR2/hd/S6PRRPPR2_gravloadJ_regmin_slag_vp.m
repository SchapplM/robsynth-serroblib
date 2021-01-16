% Calculate minimal parameter regressor of gravitation load for
% S6PRRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
% 
% Output:
% taug_reg [6x26]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 02:24
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRPPR2_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR2_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR2_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR2_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 02:21:17
% EndTime: 2021-01-16 02:21:19
% DurationCPUTime: 0.42s
% Computational Cost: add. (316->77), mult. (590->124), div. (0->0), fcn. (707->12), ass. (0->54)
t30 = qJ(3) + pkin(11);
t28 = sin(t30);
t29 = cos(t30);
t64 = cos(pkin(6));
t32 = sin(pkin(6));
t36 = sin(qJ(2));
t68 = t32 * t36;
t11 = t28 * t68 - t64 * t29;
t31 = sin(pkin(10));
t63 = cos(pkin(10));
t50 = t64 * t63;
t74 = cos(qJ(2));
t16 = t31 * t74 + t36 * t50;
t56 = t32 * t63;
t43 = -t16 * t28 - t29 * t56;
t57 = t31 * t64;
t14 = t36 * t57 - t63 * t74;
t69 = t31 * t32;
t49 = t14 * t28 + t29 * t69;
t88 = g(1) * t49 + g(2) * t43 - g(3) * t11;
t35 = sin(qJ(3));
t38 = cos(qJ(3));
t87 = -t35 * t68 + t64 * t38;
t67 = t32 * t38;
t86 = t14 * t35 + t31 * t67;
t12 = t64 * t28 + t29 * t68;
t48 = -t14 * t29 + t28 * t69;
t7 = t16 * t29 - t28 * t56;
t85 = -g(1) * t48 - g(2) * t7 - g(3) * t12;
t78 = pkin(4) * t29;
t75 = g(3) * t32;
t34 = sin(qJ(6));
t71 = t28 * t34;
t37 = cos(qJ(6));
t70 = t28 * t37;
t17 = t63 * t36 + t74 * t57;
t27 = t38 * pkin(3) + pkin(2);
t33 = qJ(4) + pkin(8);
t66 = -t14 * t33 - t17 * t27;
t59 = t32 * t74;
t65 = t27 * t59 + t33 * t68;
t60 = t28 * t74;
t15 = t31 * t36 - t74 * t50;
t58 = -t15 * t27 + t16 * t33;
t54 = t86 * pkin(3);
t52 = -qJ(5) * t28 - t78;
t51 = t87 * pkin(3);
t44 = -t16 * t35 - t38 * t56;
t42 = -g(1) * t14 + g(2) * t16 + g(3) * t68;
t40 = t44 * pkin(3);
t39 = g(1) * t17 + g(2) * t15 - g(3) * t59;
t2 = t39 * t29;
t1 = t39 * t28;
t3 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, t39, t42, 0, 0, 0, 0, 0, t39 * t38, -t39 * t35, t2, -t1, -t42, -g(1) * t66 - g(2) * t58 - g(3) * t65, -t42, -t2, t1, -g(1) * (t52 * t17 + t66) - g(2) * (t52 * t15 + t58) - g(3) * ((qJ(5) * t60 + t74 * t78) * t32 + t65), 0, 0, 0, 0, 0, -g(1) * (-t14 * t37 - t17 * t71) - g(2) * (-t15 * t71 + t16 * t37) - (t34 * t60 + t36 * t37) * t75, -g(1) * (t14 * t34 - t17 * t70) - g(2) * (-t15 * t70 - t16 * t34) - (-t34 * t36 + t37 * t60) * t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t86 - g(2) * t44 - g(3) * t87, -g(1) * (t14 * t38 - t35 * t69) - g(2) * (-t16 * t38 + t35 * t56) - g(3) * (-t64 * t35 - t36 * t67), -t88, -t85, 0, -g(1) * t54 - g(2) * t40 - g(3) * t51, 0, t88, t85, -g(1) * (pkin(4) * t49 + qJ(5) * t48 + t54) - g(2) * (pkin(4) * t43 + t7 * qJ(5) + t40) - g(3) * (-t11 * pkin(4) + t12 * qJ(5) + t51), 0, 0, 0, 0, 0, t85 * t34, t85 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t39, 0, 0, 0, -t39, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t88, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39 * t34 + t37 * t88, -t34 * t88 + t39 * t37;];
taug_reg = t3;
