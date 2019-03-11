% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RRPPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta5]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x27]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRPPPR2_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR2_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR2_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR2_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:12:07
% EndTime: 2019-03-09 08:12:09
% DurationCPUTime: 0.86s
% Computational Cost: add. (1434->140), mult. (3292->266), div. (0->0), fcn. (3071->8), ass. (0->92)
t69 = sin(pkin(10));
t71 = cos(pkin(10));
t73 = sin(qJ(6));
t75 = cos(qJ(6));
t108 = -t73 * t69 + t75 * t71;
t50 = t108 * qJD(6);
t97 = t69 ^ 2 + t71 ^ 2;
t57 = t97 * qJD(5);
t107 = 2 * qJD(4);
t70 = sin(pkin(9));
t72 = cos(pkin(9));
t74 = sin(qJ(2));
t76 = cos(qJ(2));
t51 = t70 * t74 - t72 * t76;
t106 = pkin(8) * t51;
t65 = -t72 * pkin(2) - pkin(3);
t62 = -qJ(5) + t65;
t105 = -pkin(8) + t62;
t94 = t76 * qJD(2);
t95 = t74 * qJD(2);
t48 = -t70 * t95 + t72 * t94;
t104 = t69 * t48;
t53 = t70 * t76 + t72 * t74;
t47 = t53 * qJD(2);
t103 = t71 * t47;
t102 = t71 * t48;
t99 = pkin(3) + qJ(5);
t98 = -qJ(3) - pkin(7);
t66 = pkin(2) * t95;
t79 = -t48 * qJ(4) - t53 * qJD(4) + t66;
t13 = t51 * qJD(5) + t99 * t47 + t79;
t90 = qJD(2) * t98;
t43 = t76 * qJD(3) + t74 * t90;
t44 = -t74 * qJD(3) + t76 * t90;
t31 = t70 * t43 - t72 * t44;
t22 = t48 * pkin(4) + t31;
t7 = t71 * t13 + t69 * t22;
t92 = -t76 * pkin(2) - pkin(1);
t80 = -t53 * qJ(4) + t92;
t28 = t99 * t51 + t80;
t58 = t98 * t74;
t59 = t98 * t76;
t38 = -t72 * t58 - t70 * t59;
t33 = t53 * pkin(4) + t38;
t11 = t71 * t28 + t69 * t33;
t63 = t70 * pkin(2) + qJ(4);
t96 = t63 * qJD(4);
t93 = -0.2e1 * pkin(1) * qJD(2);
t91 = -pkin(5) * t71 - pkin(4);
t18 = t71 * t22;
t6 = -t69 * t13 + t18;
t3 = t6 * t71 + t7 * t69;
t89 = -t6 * t69 + t7 * t71;
t30 = t71 * t33;
t8 = t53 * pkin(5) + t30 + (-t28 - t106) * t69;
t9 = t71 * t106 + t11;
t88 = t73 * t9 - t75 * t8;
t87 = t73 * t8 + t75 * t9;
t32 = t72 * t43 + t70 * t44;
t23 = -t47 * pkin(4) + t32;
t39 = t70 * t58 - t72 * t59;
t34 = -t51 * pkin(4) + t39;
t86 = t23 * t51 + t34 * t47;
t85 = t38 * t31 + t39 * t32;
t45 = t105 * t69;
t46 = t105 * t71;
t84 = t75 * t45 + t73 * t46;
t83 = t73 * t45 - t75 * t46;
t52 = t75 * t69 + t73 * t71;
t49 = t52 * qJD(6);
t19 = t108 * t48 - t49 * t53;
t20 = -t52 * t48 - t50 * t53;
t81 = -qJD(4) * t51 - t63 * t47;
t78 = qJD(5) * t53 - t48 * t62 - t81;
t77 = 0.2e1 * t31 * t53 - 0.2e1 * t32 * t51 + 0.2e1 * t38 * t48 - 0.2e1 * t39 * t47;
t56 = t69 * pkin(5) + t63;
t37 = t51 * pkin(3) + t80;
t36 = t52 * t51;
t35 = t108 * t51;
t26 = t47 * pkin(3) + t79;
t25 = -qJD(5) * t108 - t84 * qJD(6);
t24 = t52 * qJD(5) + t83 * qJD(6);
t21 = t91 * t51 + t39;
t16 = t91 * t47 + t32;
t15 = -t108 * t47 + t49 * t51;
t14 = t52 * t47 + t51 * t50;
t10 = -t69 * t28 + t30;
t5 = pkin(8) * t103 + t7;
t4 = t48 * pkin(5) + t18 + (-pkin(8) * t47 - t13) * t69;
t2 = -t87 * qJD(6) + t75 * t4 - t73 * t5;
t1 = t88 * qJD(6) - t73 * t4 - t75 * t5;
t12 = [0, 0, 0, 0.2e1 * t74 * t94, 0.2e1 * (-t74 ^ 2 + t76 ^ 2) * qJD(2), 0, 0, 0, t74 * t93, t76 * t93, t77, 0.2e1 * t92 * t66 + 0.2e1 * t85, t77, -0.2e1 * t26 * t51 - 0.2e1 * t37 * t47, -0.2e1 * t26 * t53 - 0.2e1 * t37 * t48, 0.2e1 * t37 * t26 + 0.2e1 * t85, 0.2e1 * t10 * t48 + 0.2e1 * t6 * t53 - 0.2e1 * t86 * t71, -0.2e1 * t11 * t48 - 0.2e1 * t7 * t53 + 0.2e1 * t86 * t69, 0.2e1 * t89 * t51 + 0.2e1 * (-t10 * t69 + t11 * t71) * t47, 0.2e1 * t10 * t6 + 0.2e1 * t11 * t7 + 0.2e1 * t34 * t23, 0.2e1 * t36 * t14, 0.2e1 * t14 * t35 - 0.2e1 * t36 * t15, 0.2e1 * t14 * t53 + 0.2e1 * t36 * t48, -0.2e1 * t15 * t53 + 0.2e1 * t35 * t48, 0.2e1 * t53 * t48, 0.2e1 * t21 * t15 - 0.2e1 * t16 * t35 + 0.2e1 * t2 * t53 - 0.2e1 * t88 * t48, 0.2e1 * t1 * t53 + 0.2e1 * t21 * t14 + 0.2e1 * t16 * t36 - 0.2e1 * t87 * t48; 0, 0, 0, 0, 0, t94, -t95, 0, -pkin(7) * t94, pkin(7) * t95 (-t47 * t70 - t48 * t72) * pkin(2) (-t31 * t72 + t32 * t70) * pkin(2), t65 * t48 + t81, t31, t32, t39 * qJD(4) + t31 * t65 + t32 * t63, t23 * t69 - t71 * t78, t23 * t71 + t69 * t78, -t3, t34 * qJD(4) + t23 * t63 + t3 * t62 + (-t10 * t71 - t11 * t69) * qJD(5), t108 * t14 - t36 * t49, -t108 * t15 - t14 * t52 - t49 * t35 - t36 * t50, t19, t20, 0, -qJD(4) * t35 + t56 * t15 + t16 * t52 + t21 * t50 + t25 * t53 - t83 * t48, qJD(4) * t36 + t108 * t16 + t56 * t14 - t21 * t49 + t24 * t53 - t84 * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t107, 0.2e1 * t96, t69 * t107, t71 * t107, 0.2e1 * t57, -0.2e1 * t62 * t57 + 0.2e1 * t96, -0.2e1 * t108 * t49, -0.2e1 * t108 * t50 + 0.2e1 * t49 * t52, 0, 0, 0, 0.2e1 * qJD(4) * t52 + 0.2e1 * t56 * t50, 0.2e1 * qJD(4) * t108 - 0.2e1 * t56 * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66, 0, -t47, -t48, t26, -t104, -t102, t97 * t47, t89, 0, 0, 0, 0, 0, t20, -t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, 0, 0, t31, t102, -t104, 0, t3, 0, 0, 0, 0, 0, t19, t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t57, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t103, t69 * t47, 0, t23, 0, 0, 0, 0, 0, t15, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4), 0, 0, 0, 0, 0, t50, -t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, -t15, t48, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49, -t50, 0, t25, t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t50, t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49, -t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t12;
