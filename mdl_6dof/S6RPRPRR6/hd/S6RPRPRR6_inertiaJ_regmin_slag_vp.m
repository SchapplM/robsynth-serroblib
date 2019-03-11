% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
% 
% Output:
% MM_reg [((6+1)*6/2)x32]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRPRR6_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR6_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR6_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t58 = sin(pkin(11));
t60 = cos(pkin(11));
t63 = sin(qJ(5));
t66 = cos(qJ(5));
t42 = t63 * t58 - t66 * t60;
t52 = -t60 * pkin(4) - pkin(3);
t36 = t42 * pkin(5) + t52;
t89 = 0.2e1 * t36;
t59 = sin(pkin(10));
t61 = cos(pkin(10));
t64 = sin(qJ(3));
t80 = cos(qJ(3));
t43 = t64 * t59 - t80 * t61;
t88 = -0.2e1 * t43;
t87 = 0.2e1 * t43;
t86 = 0.2e1 * t52;
t53 = -t61 * pkin(2) - pkin(1);
t85 = 0.2e1 * t53;
t84 = t43 * pkin(5);
t62 = sin(qJ(6));
t83 = t62 * pkin(5);
t65 = cos(qJ(6));
t82 = t65 * pkin(5);
t44 = t66 * t58 + t63 * t60;
t45 = t80 * t59 + t64 * t61;
t23 = t44 * t45;
t27 = t43 * pkin(3) - t45 * qJ(4) + t53;
t75 = pkin(7) + qJ(2);
t47 = t75 * t59;
t49 = t75 * t61;
t34 = -t64 * t47 + t80 * t49;
t17 = t58 * t27 + t60 * t34;
t77 = t58 * t45;
t12 = -pkin(8) * t77 + t17;
t16 = t60 * t27 - t58 * t34;
t76 = t60 * t45;
t9 = t43 * pkin(4) - pkin(8) * t76 + t16;
t7 = t66 * t12 + t63 * t9;
t5 = -t23 * pkin(9) + t7;
t81 = t65 * t5;
t29 = -t62 * t42 + t65 * t44;
t79 = t29 * t43;
t78 = t44 * t43;
t74 = pkin(8) + qJ(4);
t73 = t58 ^ 2 + t60 ^ 2;
t72 = t59 ^ 2 + t61 ^ 2;
t24 = t42 * t45;
t6 = -t63 * t12 + t66 * t9;
t4 = t24 * pkin(9) + t6 + t84;
t1 = t65 * t4 - t62 * t5;
t46 = t74 * t58;
t48 = t74 * t60;
t31 = -t66 * t46 - t63 * t48;
t32 = t80 * t47 + t64 * t49;
t71 = -pkin(3) * t45 - qJ(4) * t43;
t70 = t16 * t60 + t17 * t58;
t69 = -t16 * t58 + t17 * t60;
t33 = -t63 * t46 + t66 * t48;
t19 = pkin(4) * t77 + t32;
t39 = t43 ^ 2;
t35 = t42 * t43;
t28 = t65 * t42 + t62 * t44;
t22 = -t42 * pkin(9) + t33;
t21 = -t44 * pkin(9) + t31;
t18 = t28 * t43;
t15 = t23 * pkin(5) + t19;
t14 = -t62 * t23 - t65 * t24;
t13 = t65 * t23 - t62 * t24;
t11 = t62 * t21 + t65 * t22;
t10 = t65 * t21 - t62 * t22;
t2 = t62 * t4 + t81;
t3 = [1, 0, 0, 0.2e1 * pkin(1) * t61, -0.2e1 * pkin(1) * t59, 0.2e1 * t72 * qJ(2), t72 * qJ(2) ^ 2 + pkin(1) ^ 2, t45 ^ 2, t45 * t88, 0, 0, 0, t43 * t85, t45 * t85, 0.2e1 * t16 * t43 + 0.2e1 * t32 * t77, -0.2e1 * t17 * t43 + 0.2e1 * t32 * t76, -0.2e1 * t70 * t45, t16 ^ 2 + t17 ^ 2 + t32 ^ 2, t24 ^ 2, 0.2e1 * t24 * t23, -t24 * t87, t23 * t88, t39, 0.2e1 * t19 * t23 + 0.2e1 * t6 * t43, -0.2e1 * t19 * t24 - 0.2e1 * t7 * t43, t14 ^ 2, -0.2e1 * t14 * t13, t14 * t87, t13 * t88, t39, 0.2e1 * t1 * t43 + 0.2e1 * t15 * t13, 0.2e1 * t15 * t14 - 0.2e1 * t2 * t43; 0, 0, 0, -t61, t59, 0, -pkin(1), 0, 0, 0, 0, 0, t43, t45, t60 * t43, -t58 * t43, -t73 * t45, t70, 0, 0, 0, 0, 0, -t35, -t78, 0, 0, 0, 0, 0, -t18, -t79; 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, -t43, 0, -t32, -t34, -t32 * t60 + t71 * t58, t32 * t58 + t71 * t60, t69, -t32 * pkin(3) + t69 * qJ(4), -t24 * t44, -t44 * t23 + t24 * t42, t78, -t35, 0, t19 * t42 + t52 * t23 + t31 * t43, t19 * t44 - t52 * t24 - t33 * t43, t14 * t29, -t29 * t13 - t14 * t28, t79, -t18, 0, t10 * t43 + t36 * t13 + t15 * t28, -t11 * t43 + t36 * t14 + t15 * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(3) * t60, -0.2e1 * pkin(3) * t58, 0.2e1 * t73 * qJ(4), t73 * qJ(4) ^ 2 + pkin(3) ^ 2, t44 ^ 2, -0.2e1 * t44 * t42, 0, 0, 0, t42 * t86, t44 * t86, t29 ^ 2, -0.2e1 * t29 * t28, 0, 0, 0, t28 * t89, t29 * t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77, t76, 0, t32, 0, 0, 0, 0, 0, t23, -t24, 0, 0, 0, 0, 0, t13, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t60, t58, 0, -pkin(3), 0, 0, 0, 0, 0, t42, t44, 0, 0, 0, 0, 0, t28, t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, -t23, t43, t6, -t7, 0, 0, t14, -t13, t43, t43 * t82 + t1, -t81 + (-t4 - t84) * t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t42, -t44, 0, 0, 0, 0, 0, -t28, -t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, -t42, 0, t31, -t33, 0, 0, t29, -t28, 0, t10, -t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t82, -0.2e1 * t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, -t13, t43, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28, -t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, -t28, 0, t10, -t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t82, -t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t3;
