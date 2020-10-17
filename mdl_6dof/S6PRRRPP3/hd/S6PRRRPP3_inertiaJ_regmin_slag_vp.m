% Calculate minimal parameter regressor of joint inertia matrix for
% S6PRRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% MM_reg [((6+1)*6/2)x26]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PRRRPP3_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP3_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP3_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 06:53:54
% EndTime: 2019-05-05 06:53:55
% DurationCPUTime: 0.65s
% Computational Cost: add. (466->117), mult. (997->195), div. (0->0), fcn. (1081->8), ass. (0->70)
t44 = sin(qJ(4));
t47 = cos(qJ(4));
t41 = cos(pkin(6));
t45 = sin(qJ(3));
t48 = cos(qJ(3));
t63 = sin(pkin(6));
t59 = t63 * sin(qJ(2));
t19 = t41 * t45 + t48 * t59;
t56 = t63 * cos(qJ(2));
t8 = t19 * t44 + t47 * t56;
t9 = t19 * t47 - t44 * t56;
t61 = t8 * t44 + t9 * t47;
t33 = qJ(5) * t47;
t87 = -pkin(4) * t44 + t33;
t42 = pkin(4) + qJ(6);
t65 = t44 * qJ(5);
t86 = -t42 * t47 - t65;
t85 = -0.2e1 * t44;
t84 = -0.2e1 * t45;
t83 = 0.2e1 * t45;
t82 = 0.2e1 * t48;
t81 = pkin(3) * t47;
t79 = pkin(8) * t44;
t78 = pkin(9) * t48;
t18 = -t41 * t48 + t45 * t59;
t76 = t18 * t44;
t75 = t18 * t47;
t73 = t44 * t45;
t72 = t44 * t47;
t71 = t44 * t48;
t31 = t47 * t45;
t70 = t47 * t48;
t69 = t9 * qJ(5);
t35 = t45 * pkin(8);
t68 = pkin(4) * t73 + t35;
t24 = -t48 * pkin(3) - t45 * pkin(9) - pkin(2);
t67 = pkin(8) * t71 - t47 * t24;
t15 = pkin(8) * t70 + t44 * t24;
t38 = t44 ^ 2;
t40 = t47 ^ 2;
t66 = t38 + t40;
t64 = t48 * qJ(5);
t62 = t45 * t82;
t37 = t48 * pkin(4);
t13 = t37 + t67;
t60 = t18 ^ 2 + t8 ^ 2 + t9 ^ 2;
t58 = t8 * t31 - t9 * t73;
t55 = -t47 * pkin(4) - t65;
t23 = -pkin(3) + t55;
t57 = -t23 * t45 - t78;
t12 = t64 - t15;
t54 = -t12 * t47 + t13 * t44;
t53 = -pkin(5) * t73 + t15;
t1 = t18 * t73 + t8 * t48;
t2 = -t18 * t31 - t9 * t48;
t52 = -pkin(5) * t31 - t13;
t50 = qJ(5) ^ 2;
t49 = 0.2e1 * qJ(5);
t39 = t45 ^ 2;
t36 = t47 * pkin(9);
t34 = t44 * pkin(9);
t32 = -0.2e1 * t64;
t26 = t47 * pkin(5) + t36;
t25 = t44 * pkin(5) + t34;
t20 = -pkin(3) + t86;
t16 = -t45 * t33 + t68;
t11 = (qJ(6) * t44 - t33) * t45 + t68;
t5 = t53 - t64;
t3 = t48 * qJ(6) - t52;
t4 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, 0, 0, 0, t60; 0, 0, t56, -t59, 0, 0, 0, 0, 0, t48 * t56, -t45 * t56, 0, 0, 0, 0, 0, t1, -t2, t58, -t1, t2, -t9 * t12 + t8 * t13 + t18 * t16, t58, t2, t1, t18 * t11 + t8 * t3 + t9 * t5; 0, 1, 0, 0, t39, t62, 0, 0, 0, pkin(2) * t82, pkin(2) * t84, t40 * t39, -0.2e1 * t39 * t72, t70 * t84, t44 * t62, t48 ^ 2, 0.2e1 * t39 * t79 + 0.2e1 * t48 * t67, 0.2e1 * t39 * pkin(8) * t47 + 0.2e1 * t15 * t48 (t12 * t44 + t13 * t47) * t83, -0.2e1 * t13 * t48 - 0.2e1 * t16 * t73, 0.2e1 * t12 * t48 - 0.2e1 * t16 * t31, t12 ^ 2 + t13 ^ 2 + t16 ^ 2 (t3 * t47 - t44 * t5) * t83, -0.2e1 * t11 * t31 - 0.2e1 * t5 * t48, 0.2e1 * t11 * t73 + 0.2e1 * t3 * t48, t11 ^ 2 + t3 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, -t19, 0, 0, 0, 0, 0, -t75, t76, t61, t75, -t76, t61 * pkin(9) + t18 * t23, t61, -t76, -t75, t18 * t20 + t8 * t25 + t9 * t26; 0, 0, 0, 0, 0, 0, t45, t48, 0, -t35, -t48 * pkin(8), t44 * t31 (-t38 + t40) * t45, -t71, -t70, 0, -pkin(8) * t31 + (-pkin(3) * t45 + t78) * t44, pkin(9) * t70 + (t79 - t81) * t45, t54, t16 * t47 + t44 * t57, -t16 * t44 + t47 * t57, pkin(9) * t54 + t16 * t23 (t25 * t45 + t5) * t47 + (-t26 * t45 + t3) * t44, -t11 * t44 - t20 * t31 - t26 * t48, -t11 * t47 + t20 * t73 + t25 * t48, t11 * t20 + t3 * t25 + t5 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t38, 0.2e1 * t72, 0, 0, 0, 0.2e1 * t81, pkin(3) * t85, 0.2e1 * t66 * pkin(9), 0.2e1 * t23 * t47, t23 * t85, pkin(9) ^ 2 * t66 + t23 ^ 2, 0.2e1 * t25 * t44 + 0.2e1 * t26 * t47, t20 * t85, -0.2e1 * t20 * t47, t20 ^ 2 + t25 ^ 2 + t26 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, -t9, 0, t8, t9, -t8 * pkin(4) + t69, 0, t9, -t8, -t8 * t42 + t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, -t73, -t48, -t67, -t15, t55 * t45, 0.2e1 * t37 + t67, t32 + t15, -t13 * pkin(4) - t12 * qJ(5), t86 * t45, t32 + t53 (-qJ(6) - t42) * t48 + t52, t5 * qJ(5) - t3 * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, t47, 0, -t34, -t36, t87, t34, t36, t87 * pkin(9), -t42 * t44 + t33, t26, -t25, t26 * qJ(5) - t25 * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -0.2e1 * pkin(4), t49, pkin(4) ^ 2 + t50, 0, t49, 0.2e1 * t42, t42 ^ 2 + t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, 0, 0, 0, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, -t48, 0, t13, t31, 0, t48, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, 0, 0, t34, t44, 0, 0, t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(4), 0, 0, -1, -t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t73, -t48, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, 0, 0, t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, qJ(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t4;
