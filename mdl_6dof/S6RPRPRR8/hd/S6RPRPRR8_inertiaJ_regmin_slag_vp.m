% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
% 
% Output:
% MM_reg [((6+1)*6/2)x29]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRPRR8_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR8_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR8_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 19:21:37
% EndTime: 2019-05-05 19:21:39
% DurationCPUTime: 0.62s
% Computational Cost: add. (566->95), mult. (1041->152), div. (0->0), fcn. (1285->8), ass. (0->62)
t45 = sin(pkin(10));
t46 = cos(pkin(10));
t51 = cos(qJ(3));
t66 = sin(qJ(3));
t28 = -t45 * t66 + t46 * t51;
t29 = -t45 * t51 - t46 * t66;
t76 = (t28 * t46 - t29 * t45) * pkin(3);
t26 = t28 ^ 2;
t27 = t29 ^ 2;
t75 = t27 + t26;
t47 = sin(qJ(6));
t48 = sin(qJ(5));
t49 = cos(qJ(6));
t50 = cos(qJ(5));
t31 = t47 * t48 - t49 * t50;
t65 = t28 * t31;
t74 = 0.2e1 * t65;
t39 = -t46 * pkin(3) - pkin(4);
t34 = -t50 * pkin(5) + t39;
t73 = 0.2e1 * t34;
t72 = 2 * qJ(2);
t71 = t29 * pkin(5);
t70 = t47 * pkin(5);
t69 = t49 * pkin(5);
t40 = t66 * pkin(3) + qJ(2);
t16 = -t29 * pkin(4) - t28 * pkin(8) + t40;
t52 = -pkin(1) - pkin(7);
t58 = t66 * t52;
t33 = -t66 * qJ(4) + t58;
t57 = (-qJ(4) + t52) * t51;
t19 = t46 * t33 + t45 * t57;
t61 = t50 * t19;
t5 = t61 + (-pkin(9) * t28 + t16) * t48;
t68 = t49 * t5;
t38 = t45 * pkin(3) + pkin(8);
t67 = pkin(9) + t38;
t32 = t47 * t50 + t49 * t48;
t9 = t28 * t32;
t10 = t32 * t29;
t64 = t48 * t28;
t63 = t48 * t29;
t62 = t48 * t50;
t60 = t50 * t28;
t22 = t50 * t29;
t6 = t50 * t16 - t48 * t19;
t4 = -pkin(9) * t60 + t6 - t71;
t1 = t49 * t4 - t47 * t5;
t17 = t45 * t33 - t46 * t57;
t56 = t17 * t28 + t19 * t29;
t55 = t28 * t39 + t29 * t38;
t44 = t50 ^ 2;
t43 = t48 ^ 2;
t25 = t67 * t50;
t24 = t67 * t48;
t20 = t31 * t29;
t15 = -t47 * t24 + t49 * t25;
t14 = -t49 * t24 - t47 * t25;
t12 = -t49 * t22 + t47 * t63;
t8 = pkin(5) * t64 + t17;
t7 = t48 * t16 + t61;
t2 = t47 * t4 + t68;
t3 = [1, 0, 0, -2 * pkin(1), t72, pkin(1) ^ 2 + qJ(2) ^ 2, t51 ^ 2, -0.2e1 * t51 * t66, 0, 0, 0, t66 * t72, t51 * t72, 0.2e1 * t56, t17 ^ 2 + t19 ^ 2 + t40 ^ 2, t44 * t26, -0.2e1 * t26 * t62, -0.2e1 * t28 * t22, 0.2e1 * t28 * t63, t27, 0.2e1 * t17 * t64 - 0.2e1 * t6 * t29, 0.2e1 * t17 * t60 + 0.2e1 * t7 * t29, t65 ^ 2, t9 * t74, t29 * t74, 0.2e1 * t9 * t29, t27, -0.2e1 * t1 * t29 + 0.2e1 * t8 * t9, 0.2e1 * t2 * t29 - 0.2e1 * t65 * t8; 0, 0, 0, 1, 0, -pkin(1), 0, 0, 0, 0, 0, 0, 0, -t75, -t56, 0, 0, 0, 0, 0, -t75 * t48, -t75 * t50, 0, 0, 0, 0, 0, -t10 * t29 - t28 * t9, t12 * t29 + t28 * t65; 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, t75, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, t51, -t66, 0, t51 * t52, -t58, -t76 (-t17 * t46 + t19 * t45) * pkin(3), t48 * t60 (-t43 + t44) * t28, -t63, -t22, 0, -t17 * t50 + t55 * t48, t17 * t48 + t55 * t50, -t65 * t32, t31 * t65 - t32 * t9, -t10, t20, 0, -t14 * t29 + t8 * t31 + t34 * t9, t15 * t29 + t8 * t32 - t34 * t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, -t66, 0, t76, 0, 0, 0, 0, 0, t60, -t64, 0, 0, 0, 0, 0, -t65, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t45 ^ 2 + t46 ^ 2) * pkin(3) ^ 2, t43, 0.2e1 * t62, 0, 0, 0, -0.2e1 * t39 * t50, 0.2e1 * t39 * t48, t32 ^ 2, -0.2e1 * t32 * t31, 0, 0, 0, t31 * t73, t32 * t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, 0, 0, 0, 0, 0, -t22, t63, 0, 0, 0, 0, 0, t20, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, -t64, -t29, t6, -t7, 0, 0, -t65, -t9, -t29, -t29 * t69 + t1, -t68 + (-t4 + t71) * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63, t22, 0, 0, 0, 0, 0, t10, -t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, t50, 0, -t48 * t38, -t50 * t38, 0, 0, t32, -t31, 0, t14, -t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, -t48, 0, 0, 0, 0, 0, -t31, -t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t69, -0.2e1 * t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t65, -t9, -t29, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, -t31, 0, t14, -t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, -t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t69, -t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t3;
