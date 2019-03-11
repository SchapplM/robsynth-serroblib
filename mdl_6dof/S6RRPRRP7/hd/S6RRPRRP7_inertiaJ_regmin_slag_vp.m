% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRPRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x32]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRPRRP7_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP7_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP7_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t45 = sin(qJ(5));
t40 = t45 ^ 2;
t48 = cos(qJ(5));
t43 = t48 ^ 2;
t70 = t40 + t43;
t97 = t70 * pkin(9);
t46 = sin(qJ(4));
t47 = sin(qJ(2));
t49 = cos(qJ(4));
t50 = cos(qJ(2));
t18 = t47 * t46 + t50 * t49;
t19 = -t50 * t46 + t47 * t49;
t55 = t18 * t46 + t19 * t49;
t96 = t55 * t48;
t51 = -pkin(2) - pkin(3);
t25 = t49 * qJ(3) + t46 * t51;
t23 = -pkin(9) + t25;
t63 = t70 * t23;
t95 = 0.2e1 * t18;
t94 = 0.2e1 * t19;
t93 = -0.2e1 * t45;
t92 = -0.2e1 * t47;
t91 = 0.2e1 * t48;
t90 = 0.2e1 * t50;
t89 = pkin(9) * t18;
t88 = t18 * pkin(5);
t87 = t45 * pkin(9);
t86 = t48 * pkin(9);
t36 = t47 * pkin(7);
t29 = -t47 * pkin(8) + t36;
t37 = t50 * pkin(7);
t30 = -t50 * pkin(8) + t37;
t12 = -t49 * t29 + t46 * t30;
t58 = pkin(5) * t45 - t48 * qJ(6);
t6 = t58 * t19 + t12;
t85 = t6 * t45;
t84 = t6 * t48;
t24 = t46 * qJ(3) - t49 * t51;
t22 = pkin(4) + t24;
t83 = pkin(4) + t22;
t13 = t46 * t29 + t49 * t30;
t27 = -t50 * pkin(2) - t47 * qJ(3) - pkin(1);
t16 = t50 * pkin(3) - t27;
t9 = t18 * pkin(4) - t19 * pkin(9) + t16;
t5 = t48 * t13 + t45 * t9;
t82 = t12 * t45;
t81 = t12 * t48;
t80 = t18 * t23;
t79 = t45 * t18;
t78 = t45 * t19;
t77 = t45 * t23;
t76 = t45 * t46;
t75 = t45 * t48;
t74 = t48 * t18;
t15 = t48 * t19;
t73 = t48 * t23;
t72 = t48 * t46;
t59 = -t48 * pkin(5) - t45 * qJ(6);
t26 = -pkin(4) + t59;
t14 = t24 - t26;
t71 = -t14 + t26;
t42 = t47 ^ 2;
t69 = t50 ^ 2 + t42;
t68 = t18 * qJ(6);
t67 = -0.2e1 * t19 * t18;
t66 = -0.2e1 * t75;
t65 = t45 * t15;
t64 = t45 * t13 - t48 * t9;
t20 = t70 * t46;
t62 = -pkin(4) * t19 - t89;
t2 = t68 + t5;
t3 = t64 - t88;
t1 = t2 * t48 + t3 * t45;
t61 = -t19 * t26 + t89;
t60 = -t47 * pkin(2) + t50 * qJ(3);
t57 = t14 * t19 - t80;
t56 = t19 * t22 - t80;
t54 = t55 * t45;
t33 = t49 * t48;
t32 = t49 * t45;
t31 = 0.2e1 * t75;
t17 = t19 ^ 2;
t10 = (t40 - t43) * t19;
t4 = [1, 0, 0, t42, t47 * t90, 0, 0, 0, pkin(1) * t90, pkin(1) * t92, -0.2e1 * t27 * t50, 0.2e1 * t69 * pkin(7), t27 * t92, t69 * pkin(7) ^ 2 + t27 ^ 2, t17, t67, 0, 0, 0, t16 * t95, t16 * t94, t43 * t17, t17 * t66, t15 * t95, t45 * t67, t18 ^ 2, 0.2e1 * t12 * t78 - 0.2e1 * t18 * t64, 0.2e1 * t12 * t15 - 0.2e1 * t5 * t18, -0.2e1 * t3 * t18 + 0.2e1 * t6 * t78 (-t2 * t45 + t3 * t48) * t94, -0.2e1 * t15 * t6 + 0.2e1 * t2 * t18, t2 ^ 2 + t3 ^ 2 + t6 ^ 2; 0, 0, 0, 0, 0, t47, t50, 0, -t36, -t37, -t36, t60, t37, t60 * pkin(7), 0, 0, -t19, t18, 0, t12, t13, -t65, t10, -t79, -t74, 0, t56 * t45 + t81, t48 * t56 - t82, t45 * t57 + t84, -t1, -t48 * t57 + t85, t1 * t23 + t6 * t14; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(2), 0, 0.2e1 * qJ(3), pkin(2) ^ 2 + qJ(3) ^ 2, 0, 0, 0, 0, 1, 0.2e1 * t24, 0.2e1 * t25, t40, t31, 0, 0, 0, t22 * t91, t22 * t93, t14 * t91, -0.2e1 * t63, 0.2e1 * t14 * t45, t23 ^ 2 * t70 + t14 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, 0, t36, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t54, -t96, -t54, 0, t96, t1 * t46 - t6 * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(2), 0, 0, 0, 0, 0, -t49, t46, 0, 0, 0, 0, 0, -t33, t32, -t33, -t20, -t32, -t14 * t49 + t20 * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46 ^ 2 * t70 + t49 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, -t18, 0, -t12, -t13, t65, -t10, t79, t74, 0, t62 * t45 - t81, t48 * t62 + t82, -t45 * t61 - t84, t1, t48 * t61 - t85, pkin(9) * t1 + t6 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -t24, -t25, -t40, t66, 0, 0, 0, -t83 * t48, t83 * t45, t71 * t48, t63 - t97, t71 * t45, pkin(9) * t63 + t14 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, -t46, 0, 0, 0, 0, 0, t33, -t32, t33, t20, t32, pkin(9) * t20 - t49 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t40, t31, 0, 0, 0, pkin(4) * t91, pkin(4) * t93, -0.2e1 * t26 * t48, 0.2e1 * t97, t26 * t93, pkin(9) ^ 2 * t70 + t26 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, -t78, t18, -t64, -t5, -t64 + 0.2e1 * t88, t59 * t19, 0.2e1 * t68 + t5, -t3 * pkin(5) + t2 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t45, -t48, 0, -t77, -t73, -t77, t58, t73, -t58 * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t76, -t72, -t76, 0, t72, -t58 * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, t48, 0, -t87, -t86, -t87, -t58, t86, -t58 * pkin(9); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(5), 0, 0.2e1 * qJ(6), pkin(5) ^ 2 + qJ(6) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, t15, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t45, 0, t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, 0, t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t4;
