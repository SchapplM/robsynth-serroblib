% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S4RRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,d2,d3,d4]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x24]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4RRRR6_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR6_inertiaDJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR6_inertiaDJ_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4RRRR6_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:30:51
% EndTime: 2019-12-31 17:30:55
% DurationCPUTime: 1.04s
% Computational Cost: add. (831->156), mult. (2587->337), div. (0->0), fcn. (2340->8), ass. (0->98)
t53 = sin(qJ(3));
t117 = -0.4e1 * t53;
t54 = sin(qJ(2));
t114 = pkin(1) * t54;
t51 = cos(pkin(4));
t50 = sin(pkin(4));
t57 = cos(qJ(2));
t109 = t50 * t57;
t91 = pkin(6) * t109;
t30 = t91 + (pkin(7) + t114) * t51;
t31 = (-pkin(2) * t57 - pkin(7) * t54 - pkin(1)) * t50;
t56 = cos(qJ(3));
t116 = t56 * t30 + t53 * t31;
t55 = cos(qJ(4));
t48 = t55 ^ 2;
t52 = sin(qJ(4));
t105 = t52 ^ 2 - t48;
t79 = t105 * qJD(4);
t32 = (pkin(2) * t54 - pkin(7) * t57) * t50 * qJD(2);
t102 = qJD(2) * t57;
t103 = qJD(2) * t54;
t84 = t50 * t103;
t33 = -t51 * pkin(1) * t102 + pkin(6) * t84;
t8 = -qJD(3) * t116 + t56 * t32 + t53 * t33;
t115 = 0.2e1 * t50;
t113 = pkin(7) * t50;
t112 = t52 * pkin(7);
t110 = t50 * t54;
t35 = t53 * t110 - t51 * t56;
t101 = qJD(3) * t35;
t83 = t50 * t102;
t19 = t56 * t83 - t101;
t36 = t56 * t110 + t51 * t53;
t20 = t55 * t109 + t36 * t52;
t11 = -t20 * qJD(4) + t19 * t55 + t52 * t84;
t111 = t11 * t52;
t108 = t53 * t55;
t107 = t55 * t56;
t47 = t53 ^ 2;
t104 = -t56 ^ 2 + t47;
t100 = qJD(3) * t53;
t99 = qJD(3) * t55;
t98 = qJD(3) * t56;
t97 = qJD(3) * t57;
t96 = qJD(4) * t52;
t95 = qJD(4) * t55;
t94 = qJD(4) * t56;
t93 = t56 * t112;
t92 = pkin(7) * t107;
t90 = -0.2e1 * pkin(2) * qJD(3);
t89 = -0.2e1 * pkin(3) * qJD(4);
t88 = t52 * t109;
t45 = t50 ^ 2;
t87 = t45 * t102;
t86 = t52 * t94;
t85 = t55 * t94;
t82 = t52 * t95;
t81 = t53 * t98;
t80 = t55 * t98;
t78 = t104 * qJD(3);
t77 = 0.2e1 * t81;
t76 = t54 * t87;
t75 = t52 * t80;
t74 = -t56 * pkin(3) - t53 * pkin(8);
t73 = pkin(3) * t53 - pkin(8) * t56;
t29 = pkin(6) * t110 + (-pkin(1) * t57 - pkin(2)) * t51;
t13 = t35 * pkin(3) - t36 * pkin(8) + t29;
t15 = -pkin(8) * t109 + t116;
t4 = t52 * t13 + t55 * t15;
t21 = t36 * t55 - t88;
t72 = -t20 * t55 - t21 * t52;
t70 = -t53 * t30 + t56 * t31;
t41 = -pkin(2) + t74;
t28 = t52 * t41 + t92;
t14 = pkin(3) * t109 - t70;
t6 = -pkin(3) * t84 - t8;
t68 = t14 * t95 + t6 * t52;
t67 = t14 * t96 - t6 * t55;
t66 = t73 * t52;
t18 = t36 * qJD(3) + t53 * t83;
t65 = t52 * t18 + t35 * t95;
t64 = -t55 * t18 + t35 * t96;
t7 = t30 * t100 - t31 * t98 - t53 * t32 + t56 * t33;
t63 = t56 * t103 + t53 * t97;
t62 = t53 * t103 - t56 * t97;
t61 = -t53 * t96 + t80;
t60 = t53 * t99 + t86;
t34 = (t51 * t114 + t91) * qJD(2);
t59 = pkin(8) * t84 - t7;
t58 = t18 * pkin(3) - t19 * pkin(8) + t34;
t27 = t55 * t41 - t93;
t17 = -t28 * qJD(4) + (t53 * t112 + t55 * t73) * qJD(3);
t16 = t60 * pkin(7) - qJD(3) * t66 - t41 * t95;
t10 = -qJD(4) * t88 + t19 * t52 + t36 * t95 - t55 * t84;
t3 = t55 * t13 - t52 * t15;
t2 = -t4 * qJD(4) - t52 * t59 + t55 * t58;
t1 = -t13 * t95 + t15 * t96 - t52 * t58 - t55 * t59;
t5 = [0, 0, 0, 0.2e1 * t76, 0.2e1 * (-t54 ^ 2 + t57 ^ 2) * t45 * qJD(2), 0.2e1 * t51 * t83, -0.2e1 * t51 * t84, 0, -0.2e1 * t45 * pkin(1) * t103 - 0.2e1 * t34 * t51, -0.2e1 * pkin(1) * t87 + 0.2e1 * t33 * t51, 0.2e1 * t36 * t19, -0.2e1 * t36 * t18 - 0.2e1 * t19 * t35, (t36 * t103 - t19 * t57) * t115, (-t35 * t103 + t18 * t57) * t115, -0.2e1 * t76, 0.2e1 * t29 * t18 + 0.2e1 * t34 * t35 + 0.2e1 * (t70 * t103 - t8 * t57) * t50, 0.2e1 * t29 * t19 + 0.2e1 * t34 * t36 + 0.2e1 * (-t103 * t116 - t7 * t57) * t50, 0.2e1 * t21 * t11, -0.2e1 * t21 * t10 - 0.2e1 * t11 * t20, 0.2e1 * t11 * t35 + 0.2e1 * t21 * t18, -0.2e1 * t10 * t35 - 0.2e1 * t20 * t18, 0.2e1 * t35 * t18, 0.2e1 * t14 * t10 + 0.2e1 * t3 * t18 + 0.2e1 * t2 * t35 + 0.2e1 * t6 * t20, 0.2e1 * t1 * t35 + 0.2e1 * t14 * t11 - 0.2e1 * t4 * t18 + 0.2e1 * t6 * t21; 0, 0, 0, 0, 0, t83, -t84, 0, -t34, t33, t19 * t53 + t36 * t98, -t53 * t18 + t19 * t56 + (-t35 * t56 - t36 * t53) * qJD(3), t62 * t50, t63 * t50, 0, -pkin(2) * t18 + t29 * t100 - t62 * t113 - t34 * t56, -pkin(2) * t19 - t63 * t113 + t29 * t98 + t34 * t53, t108 * t11 + t21 * t61, t72 * t98 + (-t10 * t55 - t111 + (t20 * t52 - t21 * t55) * qJD(4)) * t53, (t35 * t99 - t11) * t56 + (qJD(3) * t21 - t64) * t53, (-t101 * t52 + t10) * t56 + (-qJD(3) * t20 - t65) * t53, t100 * t35 - t18 * t56, t17 * t35 + t27 * t18 + (-t2 + (pkin(7) * t20 + t14 * t52) * qJD(3)) * t56 + (pkin(7) * t10 + qJD(3) * t3 + t68) * t53, t16 * t35 - t28 * t18 + (-t1 + (pkin(7) * t21 + t14 * t55) * qJD(3)) * t56 + (pkin(7) * t11 - qJD(3) * t4 - t67) * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77, -0.2e1 * t78, 0, 0, 0, t53 * t90, t56 * t90, -0.2e1 * t47 * t82 + 0.2e1 * t48 * t81, t75 * t117 + 0.2e1 * t47 * t79, 0.2e1 * t104 * t99 + 0.2e1 * t53 * t86, -0.2e1 * t52 * t78 + 0.2e1 * t53 * t85, -0.2e1 * t81, 0.2e1 * t27 * t100 - 0.2e1 * t17 * t56 + 0.2e1 * (t47 * t95 + t52 * t77) * pkin(7), -0.2e1 * t28 * t100 - 0.2e1 * t16 * t56 + 0.2e1 * (-t47 * t96 + t55 * t77) * pkin(7); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, -t18, t84, t8, t7, t21 * t95 + t111, qJD(4) * t72 - t52 * t10 + t11 * t55, t65, -t64, 0, -pkin(3) * t10 - pkin(8) * t65 + t67, -pkin(3) * t11 + pkin(8) * t64 + t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t98, -t100, 0, -pkin(7) * t98, pkin(7) * t100, -t53 * t79 + t75, -t105 * t98 + t82 * t117, t100 * t52 - t85, t60, 0, (pkin(8) * t107 + (-t55 * pkin(3) + t112) * t53) * qJD(4) + (t52 * t74 - t92) * qJD(3), (pkin(7) * t108 + t66) * qJD(4) + (t55 * t74 + t93) * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t82, -0.2e1 * t79, 0, 0, 0, t52 * t89, t55 * t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, -t10, t18, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, -t52 * t98 - t53 * t95, t100, t17, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t95, -t96, 0, -pkin(8) * t95, pkin(8) * t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t5;
