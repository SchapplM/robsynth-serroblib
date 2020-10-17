% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
% 
% Output:
% MM_reg [((6+1)*6/2)x28]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRPRRP5_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP5_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP5_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 17:53:56
% EndTime: 2019-05-06 17:53:59
% DurationCPUTime: 0.92s
% Computational Cost: add. (1489->142), mult. (3626->274), div. (0->0), fcn. (4099->10), ass. (0->96)
t57 = sin(pkin(11));
t58 = sin(pkin(6));
t59 = cos(pkin(11));
t63 = sin(qJ(2));
t66 = cos(qJ(2));
t35 = (t57 * t66 + t59 * t63) * t58;
t60 = cos(pkin(6));
t62 = sin(qJ(4));
t65 = cos(qJ(4));
t26 = t35 * t62 - t60 * t65;
t104 = -0.2e1 * t26;
t103 = 0.2e1 * t62;
t102 = pkin(1) * t63;
t64 = cos(qJ(5));
t101 = pkin(4) * t64;
t45 = t60 * t66 * pkin(1);
t79 = pkin(8) + qJ(3);
t93 = t58 * t63;
t28 = t60 * pkin(2) - t79 * t93 + t45;
t74 = t60 * t102;
t92 = t58 * t66;
t31 = t79 * t92 + t74;
t17 = t57 * t28 + t59 * t31;
t14 = t60 * pkin(9) + t17;
t34 = t57 * t93 - t59 * t92;
t41 = (-pkin(2) * t66 - pkin(1)) * t58;
t20 = t34 * pkin(3) - t35 * pkin(9) + t41;
t10 = -t62 * t14 + t65 * t20;
t6 = -t34 * pkin(4) - t10;
t61 = sin(qJ(5));
t100 = t6 * t61;
t99 = t6 * t64;
t98 = t61 * pkin(5);
t27 = t35 * t65 + t60 * t62;
t19 = t27 * t64 + t34 * t61;
t97 = t19 * t65;
t47 = t57 * pkin(2) + pkin(9);
t96 = t47 * t61;
t52 = t58 ^ 2;
t95 = t52 * t66;
t53 = t61 ^ 2;
t94 = t53 * t62;
t91 = t61 * t26;
t90 = t61 * t62;
t89 = t61 * t64;
t88 = t61 * t65;
t87 = t62 * t26;
t86 = t62 * t34;
t85 = t62 * t47;
t84 = t64 * t26;
t83 = t64 * t62;
t82 = t64 * t65;
t18 = t27 * t61 - t34 * t64;
t81 = t65 * t18;
t80 = t65 * t47;
t78 = -qJ(6) - pkin(10);
t77 = qJ(6) * t62;
t76 = 0.2e1 * t58 * t60;
t75 = t65 * t103;
t73 = t19 * t90;
t72 = t61 * t87;
t71 = t26 * t83;
t70 = t64 * t80;
t48 = -t59 * pkin(2) - pkin(3);
t11 = t65 * t14 + t62 * t20;
t7 = t34 * pkin(10) + t11;
t16 = t59 * t28 - t57 * t31;
t13 = -t60 * pkin(3) - t16;
t9 = t26 * pkin(4) - t27 * pkin(10) + t13;
t3 = -t61 * t7 + t64 * t9;
t1 = t26 * pkin(5) - t19 * qJ(6) + t3;
t4 = t61 * t9 + t64 * t7;
t2 = -t18 * qJ(6) + t4;
t69 = -t1 * t61 + t2 * t64;
t42 = t78 * t61;
t43 = t78 * t64;
t68 = -t42 * t61 - t43 * t64;
t56 = t65 ^ 2;
t55 = t64 ^ 2;
t54 = t62 ^ 2;
t51 = -t64 * pkin(5) - pkin(4);
t50 = t55 * t62;
t49 = t55 * t54;
t40 = -t65 * pkin(4) - t62 * pkin(10) + t48;
t39 = pkin(8) * t92 + t74;
t38 = -pkin(8) * t93 + t45;
t37 = (t47 + t98) * t62;
t36 = t64 * t40;
t32 = t65 * t34;
t25 = t61 * t40 + t70;
t24 = -t61 * t80 + t36;
t22 = t70 + (t40 - t77) * t61;
t21 = -t64 * t77 + t36 + (-pkin(5) - t96) * t65;
t15 = t18 * t83;
t5 = t18 * pkin(5) + t6;
t8 = [1, 0, 0, t52 * t63 ^ 2, 0.2e1 * t63 * t95, t63 * t76, t66 * t76, t60 ^ 2, 0.2e1 * pkin(1) * t95 + 0.2e1 * t38 * t60, -0.2e1 * t52 * t102 - 0.2e1 * t39 * t60, -0.2e1 * t16 * t35 - 0.2e1 * t17 * t34, t16 ^ 2 + t17 ^ 2 + t41 ^ 2, t27 ^ 2, t27 * t104, 0.2e1 * t27 * t34, t34 * t104, t34 ^ 2, 0.2e1 * t10 * t34 + 0.2e1 * t13 * t26, -0.2e1 * t11 * t34 + 0.2e1 * t13 * t27, t19 ^ 2, -0.2e1 * t19 * t18, 0.2e1 * t19 * t26, t18 * t104, t26 ^ 2, 0.2e1 * t6 * t18 + 0.2e1 * t3 * t26, 0.2e1 * t6 * t19 - 0.2e1 * t4 * t26, -0.2e1 * t1 * t19 - 0.2e1 * t2 * t18, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, t93, t92, t60, t38, -t39 (-t34 * t57 - t35 * t59) * pkin(2) (t16 * t59 + t17 * t57) * pkin(2), t27 * t62, t27 * t65 - t87, t86, t32, 0, -t13 * t65 + t48 * t26 - t34 * t85, t13 * t62 + t48 * t27 - t34 * t80, t19 * t83, -t15 - t73, t71 - t97, -t72 + t81, -t26 * t65, t24 * t26 - t3 * t65 + (t18 * t47 + t100) * t62, -t25 * t26 + t4 * t65 + (t19 * t47 + t99) * t62, -t22 * t18 - t21 * t19 + (-t1 * t64 - t2 * t61) * t62, t1 * t21 + t2 * t22 + t5 * t37; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t57 ^ 2 + t59 ^ 2) * pkin(2) ^ 2, t54, t75, 0, 0, 0, -0.2e1 * t48 * t65, t48 * t103, t49, -0.2e1 * t54 * t89, -0.2e1 * t62 * t82, t61 * t75, t56, -0.2e1 * t24 * t65 + 0.2e1 * t54 * t96, 0.2e1 * t54 * t47 * t64 + 0.2e1 * t25 * t65 (-t21 * t64 - t22 * t61) * t103, t21 ^ 2 + t22 ^ 2 + t37 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, 0, 0, 0, 0, 0, t32, -t86, 0, 0, 0, 0, 0, -t72 - t81, -t71 - t97, -t15 + t73, -t5 * t65 + t69 * t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37 * t65 + (-t21 * t61 + t22 * t64) * t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53 * t54 + t49 + t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, -t26, t34, t10, -t11, t19 * t61, -t61 * t18 + t19 * t64, t91, t84, 0, -pkin(4) * t18 - pkin(10) * t91 - t99, -pkin(4) * t19 - pkin(10) * t84 + t100, t43 * t18 - t42 * t19 + t69, t1 * t42 - t2 * t43 + t5 * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, t65, 0, -t85, -t80, t61 * t83, t50 - t94, -t88, -t82, 0, -t47 * t83 + (-pkin(4) * t62 + pkin(10) * t65) * t61, pkin(10) * t82 + (t96 - t101) * t62 (-t42 * t62 + t22) * t64 + (t43 * t62 - t21) * t61, t21 * t42 - t22 * t43 + t37 * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, -t62, 0, 0, 0, 0, 0, t82, -t88, t50 + t94, -t65 * t51 + t68 * t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t53, 0.2e1 * t89, 0, 0, 0, 0.2e1 * t101, -0.2e1 * pkin(4) * t61, 0.2e1 * t68, t42 ^ 2 + t43 ^ 2 + t51 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, -t18, t26, t3, -t4, -pkin(5) * t19, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t83, -t90, -t65, t24, -t25, -pkin(5) * t83, t21 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t90, -t83, 0, -pkin(5) * t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, t64, 0, -t61 * pkin(10), -t64 * pkin(10), -t98, t42 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, pkin(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t8;
