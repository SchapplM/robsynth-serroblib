% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPRRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
% 
% Output:
% MM_reg [((6+1)*6/2)x29]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRRPP5_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP5_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP5_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t44 = sin(pkin(9));
t45 = cos(pkin(9));
t47 = sin(qJ(3));
t73 = cos(qJ(3));
t28 = t73 * t44 + t47 * t45;
t81 = -0.2e1 * t28;
t46 = sin(qJ(4));
t37 = t46 * qJ(5);
t48 = cos(qJ(4));
t68 = t48 * pkin(4) + t37;
t49 = pkin(4) + pkin(5);
t80 = t48 * t49 + t37;
t66 = t48 * qJ(5);
t79 = t49 * t46 - t66;
t36 = -t45 * pkin(2) - pkin(1);
t78 = 0.2e1 * t36;
t77 = -0.2e1 * t46;
t76 = 0.2e1 * t48;
t27 = t47 * t44 - t73 * t45;
t75 = pkin(8) * t27;
t23 = t27 * pkin(4);
t74 = t46 * pkin(8);
t11 = t27 * pkin(3) - t28 * pkin(8) + t36;
t69 = pkin(7) + qJ(2);
t30 = t69 * t44;
t31 = t69 * t45;
t16 = -t47 * t30 + t73 * t31;
t7 = t46 * t11 + t48 * t16;
t18 = t46 * t27;
t72 = t46 * t28;
t71 = t46 * t48;
t19 = t48 * t27;
t20 = t48 * t28;
t6 = t48 * t11 - t46 * t16;
t67 = t44 ^ 2 + t45 ^ 2;
t42 = t46 ^ 2;
t43 = t48 ^ 2;
t34 = t42 + t43;
t65 = t48 * qJ(6);
t64 = t27 * t81;
t29 = -pkin(3) - t68;
t22 = t27 * qJ(5);
t3 = t22 + t7;
t63 = 0.2e1 * t22 + t7;
t4 = -t23 - t6;
t15 = t73 * t30 + t47 * t31;
t62 = -pkin(3) * t28 - t75;
t55 = t28 * t65 - t4;
t1 = -t27 * pkin(5) - t55;
t17 = qJ(6) * t72;
t2 = t17 + t3;
t61 = -t1 * t48 + t2 * t46;
t60 = t3 * t48 + t4 * t46;
t59 = t3 * t46 - t4 * t48;
t58 = -t28 * t29 + t75;
t57 = pkin(4) * t46 - t66;
t32 = (pkin(8) - qJ(6)) * t46;
t38 = t48 * pkin(8);
t33 = t38 - t65;
t56 = -t48 * t32 + t46 * t33;
t52 = qJ(5) ^ 2;
t51 = 0.2e1 * qJ(5);
t25 = t28 ^ 2;
t24 = t48 * pkin(5) - t29;
t12 = t34 * t28;
t8 = t57 * t28 + t15;
t5 = -t79 * t28 - t15;
t9 = [1, 0, 0, 0.2e1 * pkin(1) * t45, -0.2e1 * pkin(1) * t44, 0.2e1 * t67 * qJ(2), t67 * qJ(2) ^ 2 + pkin(1) ^ 2, t25, t64, 0, 0, 0, t27 * t78, t28 * t78, t43 * t25, -0.2e1 * t25 * t71, 0.2e1 * t27 * t20, t46 * t64, t27 ^ 2, 0.2e1 * t15 * t72 + 0.2e1 * t6 * t27, 0.2e1 * t15 * t20 - 0.2e1 * t7 * t27, -0.2e1 * t4 * t27 + 0.2e1 * t8 * t72, t59 * t81, -0.2e1 * t8 * t20 + 0.2e1 * t3 * t27, t3 ^ 2 + t4 ^ 2 + t8 ^ 2, -0.2e1 * t1 * t27 - 0.2e1 * t5 * t72, 0.2e1 * t2 * t27 + 0.2e1 * t5 * t20, 0.2e1 * t61 * t28, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, -t45, t44, 0, -pkin(1), 0, 0, 0, 0, 0, t27, t28, 0, 0, 0, 0, 0, t19, -t18, t19, -t12, t18, t59, t19, t18, t12, t61; 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, 0, 0, 0, t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, -t27, 0, -t15, -t16, t46 * t20 (-t42 + t43) * t28, t18, t19, 0, -t15 * t48 + t62 * t46, t15 * t46 + t62 * t48, -t58 * t46 - t8 * t48, t60, -t8 * t46 + t58 * t48, t60 * pkin(8) + t8 * t29, -t24 * t72 - t32 * t27 + t5 * t48, t24 * t20 + t33 * t27 + t5 * t46, -t1 * t46 - t2 * t48 + t56 * t28, t1 * t32 + t2 * t33 + t5 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t42, 0.2e1 * t71, 0, 0, 0, pkin(3) * t76, pkin(3) * t77, -0.2e1 * t29 * t48, 0.2e1 * t34 * pkin(8), t29 * t77, t34 * pkin(8) ^ 2 + t29 ^ 2, t24 * t76, 0.2e1 * t24 * t46, -0.2e1 * t32 * t46 - 0.2e1 * t33 * t48, t24 ^ 2 + t32 ^ 2 + t33 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, -t72, t27, t6, -t7, -t4 + t23, -t68 * t28, t63, -t4 * pkin(4) + t3 * qJ(5) (pkin(5) + t49) * t27 + t55, t17 + t63, t80 * t28, t2 * qJ(5) - t1 * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, -t46, t48, 0, t46, t68, t48, t46, 0, t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, t48, 0, -t74, -t38, -t74, -t57, t38, -t57 * pkin(8), -t32, t33, t79, t33 * qJ(5) - t32 * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4), 0, t51, pkin(4) ^ 2 + t52, 0.2e1 * t49, t51, 0, t49 ^ 2 + t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t27, t20, 0, t4, -t27, 0, -t20, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48, 0, 0, 0, -t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, 0, t74, 0, 0, -t46, t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(4), -1, 0, 0, -t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t72, t20, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, t46, 0, t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t9;
