% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRPPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta4]';
% 
% Output:
% MM_reg [((6+1)*6/2)x32]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRPPRR11_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR11_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR11_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t68 = sin(pkin(11));
t70 = cos(pkin(11));
t74 = sin(qJ(5));
t89 = t74 * t70;
t99 = cos(qJ(5));
t44 = t68 * t99 + t89;
t40 = t44 ^ 2;
t82 = t99 * t70;
t43 = t68 * t74 - t82;
t41 = t43 ^ 2;
t109 = -t40 - t41;
t108 = -2 * pkin(2);
t71 = cos(pkin(6));
t69 = sin(pkin(6));
t77 = cos(qJ(2));
t94 = t69 * t77;
t34 = t68 * t71 + t70 * t94;
t35 = -t68 * t94 + t70 * t71;
t18 = t34 * t99 + t35 * t74;
t107 = -0.2e1 * t18;
t106 = -0.2e1 * t43;
t105 = 0.2e1 * t44;
t104 = 0.2e1 * t69;
t103 = 2 * qJ(3);
t75 = sin(qJ(2));
t102 = pkin(1) * t75;
t101 = pkin(1) * t77;
t72 = -pkin(2) - qJ(4);
t100 = -pkin(9) + t72;
t19 = -t34 * t74 + t35 * t99;
t58 = t69 * t75;
t73 = sin(qJ(6));
t76 = cos(qJ(6));
t15 = t19 * t76 + t58 * t73;
t98 = t15 * t73;
t97 = t43 * t73;
t96 = t43 * t76;
t63 = t69 ^ 2;
t95 = t63 * t77;
t93 = t71 * t77;
t92 = t73 * t18;
t91 = t73 * t44;
t90 = t73 * t76;
t17 = t76 * t18;
t36 = t76 * t44;
t51 = pkin(8) * t58;
t83 = -pkin(2) - t101;
t23 = pkin(3) * t58 + t51 + (-qJ(4) + t83) * t71;
t81 = -qJ(3) * t75 - pkin(1);
t27 = (t72 * t77 + t81) * t69;
t13 = t23 * t68 + t27 * t70;
t38 = pkin(8) * t94 + t102 * t71;
t50 = t68 ^ 2 + t70 ^ 2;
t56 = pkin(4) * t68 + qJ(3);
t88 = t43 * t105;
t87 = 0.2e1 * t58;
t86 = t44 * t58;
t85 = t68 * t58;
t84 = t70 * t58;
t60 = t71 * qJ(3);
t29 = -t60 - t38;
t12 = t23 * t70 - t27 * t68;
t28 = pkin(3) * t94 - t29;
t80 = pkin(5) * t43 - pkin(10) * t44;
t79 = t12 * t70 + t13 * t68;
t8 = pkin(4) * t58 - pkin(9) * t35 + t12;
t9 = -pkin(9) * t34 + t13;
t5 = -t74 * t9 + t8 * t99;
t6 = t74 * t8 + t9 * t99;
t16 = pkin(4) * t34 + t28;
t78 = qJ(3) ^ 2;
t67 = t76 ^ 2;
t66 = t73 ^ 2;
t57 = t63 * t75 ^ 2;
t46 = t100 * t68;
t42 = t50 * t72;
t37 = pkin(1) * t93 - t51;
t33 = t43 * t58;
t32 = t71 * t83 + t51;
t30 = (-pkin(2) * t77 + t81) * t69;
t26 = t100 * t89 + t46 * t99;
t25 = -t100 * t82 + t46 * t74;
t20 = pkin(5) * t44 + pkin(10) * t43 + t56;
t14 = t19 * t73 - t58 * t76;
t11 = t20 * t73 + t26 * t76;
t10 = t20 * t76 - t26 * t73;
t7 = pkin(5) * t18 - pkin(10) * t19 + t16;
t4 = pkin(10) * t58 + t6;
t3 = -pkin(5) * t58 - t5;
t2 = t4 * t76 + t7 * t73;
t1 = -t4 * t73 + t7 * t76;
t21 = [1, 0, 0, t57, 0.2e1 * t75 * t95, t71 * t87, t93 * t104, t71 ^ 2, 0.2e1 * pkin(1) * t95 + 0.2e1 * t37 * t71, -0.2e1 * t102 * t63 - 0.2e1 * t38 * t71 (-t29 * t77 + t32 * t75) * t104, 0.2e1 * t30 * t94 + 0.2e1 * t32 * t71, -0.2e1 * t29 * t71 - 0.2e1 * t30 * t58, t29 ^ 2 + t30 ^ 2 + t32 ^ 2, 0.2e1 * t12 * t58 + 0.2e1 * t28 * t34, -0.2e1 * t13 * t58 + 0.2e1 * t28 * t35, -0.2e1 * t12 * t35 - 0.2e1 * t13 * t34, t12 ^ 2 + t13 ^ 2 + t28 ^ 2, t19 ^ 2, t19 * t107, t19 * t87, t58 * t107, t57, 0.2e1 * t16 * t18 + 0.2e1 * t5 * t58, 0.2e1 * t16 * t19 - 0.2e1 * t58 * t6, t15 ^ 2, -0.2e1 * t15 * t14, 0.2e1 * t15 * t18, t14 * t107, t18 ^ 2, 0.2e1 * t1 * t18 + 0.2e1 * t14 * t3, 0.2e1 * t15 * t3 - 0.2e1 * t18 * t2; 0, 0, 0, 0, 0, t58, t94, t71, t37, -t38 (-pkin(2) * t75 + qJ(3) * t77) * t69, t51 + (t108 - t101) * t71, 0.2e1 * t60 + t38, -pkin(2) * t32 - qJ(3) * t29, qJ(3) * t34 + t28 * t68 + t72 * t84, qJ(3) * t35 + t28 * t70 - t72 * t85 (-t35 * t72 - t12) * t70 + (-t34 * t72 - t13) * t68, t28 * qJ(3) + t72 * t79, -t19 * t43, t18 * t43 - t19 * t44, -t33, -t86, 0, t16 * t44 + t18 * t56 - t25 * t58, -t16 * t43 + t19 * t56 - t26 * t58, -t15 * t96 (t14 * t76 + t98) * t43, t15 * t44 - t17 * t43, -t14 * t44 + t43 * t92, t18 * t44, t1 * t44 + t10 * t18 + t14 * t25 - t3 * t97, -t11 * t18 + t15 * t25 - t2 * t44 - t3 * t96; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, t108, t103, pkin(2) ^ 2 + t78, t68 * t103, t70 * t103, -0.2e1 * t42, t50 * (t72 ^ 2) + t78, t41, t88, 0, 0, 0, t56 * t105, t56 * t106, t67 * t41, -0.2e1 * t41 * t90, t36 * t106, t73 * t88, t40, 0.2e1 * t10 * t44 - 0.2e1 * t25 * t97, -0.2e1 * t11 * t44 - 0.2e1 * t25 * t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, t71, 0, t32, t84, -t85, -t34 * t68 - t35 * t70, t79, 0, 0, 0, 0, 0, -t33, -t86, 0, 0, 0, 0, 0, t14 * t43 - t18 * t91, t15 * t43 - t18 * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(2), 0, 0, -t50, t42, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t109 * t73, t109 * t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, t50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, t35, 0, t28, 0, 0, 0, 0, 0, t18, t19, 0, 0, 0, 0, 0, t17, -t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68, t70, 0, qJ(3), 0, 0, 0, 0, 0, t44, -t43, 0, 0, 0, 0, 0, t36, -t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, -t18, t58, t5, -t6, t98, -t14 * t73 + t15 * t76, t92, t17, 0, -pkin(5) * t14 - pkin(10) * t92 - t3 * t76, -pkin(5) * t15 - pkin(10) * t17 + t3 * t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43, -t44, 0, -t25, -t26, -t43 * t90 (t66 - t67) * t43, t91, t36, 0, -t25 * t76 + t73 * t80, t25 * t73 + t76 * t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43, -t44, 0, 0, 0, 0, 0, -t96, t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t66, 0.2e1 * t90, 0, 0, 0, 0.2e1 * pkin(5) * t76, -0.2e1 * pkin(5) * t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, -t14, t18, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t96, t97, t44, t10, -t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t91, -t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t76, -t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73, t76, 0, -t73 * pkin(10), -t76 * pkin(10); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t21;
