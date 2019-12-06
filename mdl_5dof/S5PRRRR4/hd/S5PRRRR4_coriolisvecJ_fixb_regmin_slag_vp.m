% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5PRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% 
% Output:
% tauc_reg [5x21]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:08
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRRRR4_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR4_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR4_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR4_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:07:58
% EndTime: 2019-12-05 17:08:01
% DurationCPUTime: 0.61s
% Computational Cost: add. (653->112), mult. (1181->175), div. (0->0), fcn. (752->6), ass. (0->95)
t119 = pkin(7) + pkin(8);
t63 = sin(qJ(5));
t64 = sin(qJ(4));
t66 = cos(qJ(5));
t67 = cos(qJ(4));
t36 = t63 * t67 + t66 * t64;
t60 = qJD(2) + qJD(3);
t28 = t36 * t60;
t59 = qJD(4) + qJD(5);
t118 = qJD(5) - t59;
t100 = pkin(2) * qJD(2);
t65 = sin(qJ(3));
t82 = t65 * t100 + t119 * t60;
t21 = t67 * qJD(1) - t82 * t64;
t22 = t64 * qJD(1) + t82 * t67;
t68 = cos(qJ(3));
t116 = t68 * pkin(2);
t52 = t65 * pkin(2) + pkin(7);
t115 = -pkin(8) - t52;
t18 = t59 * t36;
t54 = -t67 * pkin(4) - pkin(3);
t88 = t68 * t100;
t29 = t54 * t60 - t88;
t99 = pkin(2) * qJD(3);
t84 = qJD(2) * t99;
t79 = t65 * t84;
t97 = t64 * qJD(4);
t87 = t60 * t97;
t30 = pkin(4) * t87 + t79;
t104 = t66 * t67;
t109 = t63 * t64;
t35 = -t104 + t109;
t114 = t29 * t18 + t30 * t35;
t17 = t59 * t35;
t113 = -t29 * t17 + t30 * t36;
t13 = t17 * t59;
t91 = t60 * t104;
t92 = t60 * t109;
t26 = -t91 + t92;
t112 = t28 * t26;
t111 = t59 * t68;
t110 = t60 * t64;
t107 = t65 * t67;
t106 = t66 * t22;
t69 = qJD(4) ^ 2;
t103 = t69 * t64;
t56 = t69 * t67;
t42 = -t60 * pkin(3) - t88;
t96 = t67 * qJD(4);
t102 = t42 * t96 + t64 * t79;
t101 = t64 ^ 2 - t67 ^ 2;
t95 = -qJD(2) - t60;
t94 = -qJD(3) + t60;
t93 = pkin(4) * t110;
t90 = t68 * t99;
t89 = pkin(4) * t97;
t86 = t60 * t96;
t85 = t68 * t97;
t16 = qJD(4) * pkin(4) + t21;
t83 = -pkin(4) * t59 - t16;
t81 = qJD(4) * t119;
t80 = qJD(4) * t115;
t78 = t68 * t84;
t76 = t94 * t100;
t75 = t95 * t99;
t11 = t21 * qJD(4) + t67 * t78;
t12 = -t22 * qJD(4) - t64 * t78;
t73 = -t63 * t11 + t66 * t12 - t29 * t28;
t72 = -t65 * t110 + t68 * t96;
t71 = -t42 * t60 - t78;
t9 = qJD(5) * t91 - t59 * t92 + t66 * t86;
t70 = t29 * t26 + (t118 * t22 - t12) * t63;
t58 = t60 ^ 2;
t57 = t67 * pkin(8);
t53 = -pkin(3) - t116;
t51 = t67 * pkin(7) + t57;
t50 = t119 * t64;
t45 = t54 - t116;
t40 = 0.2e1 * t64 * t86;
t39 = t65 * t99 + t89;
t38 = t67 * t81;
t37 = t64 * t81;
t34 = t67 * t52 + t57;
t33 = t115 * t64;
t31 = t42 * t97;
t25 = -0.2e1 * t101 * t60 * qJD(4);
t24 = -t64 * t90 + t67 * t80;
t23 = t64 * t80 + t67 * t90;
t14 = t18 * t59;
t10 = t18 * t60;
t5 = -t26 ^ 2 + t28 ^ 2;
t3 = t26 * t59 + t9;
t2 = -t28 * t17 + t9 * t36;
t1 = -t36 * t10 + t17 * t26 - t28 * t18 - t9 * t35;
t4 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t103, -t56, 0, 0, 0, 0, 0, -t14, t13; 0, 0, 0, 0, 0, t65 * t75, t68 * t75, t40, t25, t56, -t103, 0, t53 * t87 - t52 * t56 + t31 + (t95 * t107 - t85) * t99, t52 * t103 + t53 * t86 - t72 * t99 + t102, t2, t1, -t13, -t14, 0, t39 * t26 + t45 * t10 + (-t63 * t23 + t66 * t24 + (-t33 * t63 - t34 * t66) * qJD(5)) * t59 + t114, t39 * t28 + t45 * t9 - (t66 * t23 + t63 * t24 + (t33 * t66 - t34 * t63) * qJD(5)) * t59 + t113; 0, 0, 0, 0, 0, t65 * t76, t68 * t76, t40, t25, t56, -t103, 0, -pkin(3) * t87 - pkin(7) * t56 + t31 + (t94 * t107 + t85) * t100, -pkin(3) * t86 + pkin(7) * t103 + t72 * t100 + t102, t2, t1, -t13, -t14, 0, t26 * t89 + t54 * t10 + (t63 * t37 - t66 * t38 + (t50 * t63 - t51 * t66) * qJD(5)) * t59 + (t36 * t111 - t65 * t26) * t100 + t114, t28 * t89 + t54 * t9 - (-t66 * t37 - t63 * t38 + (-t50 * t66 - t51 * t63) * qJD(5)) * t59 + (-t35 * t111 - t65 * t28) * t100 + t113; 0, 0, 0, 0, 0, 0, 0, -t64 * t58 * t67, t101 * t58, 0, 0, 0, t71 * t64, t71 * t67, t112, t5, t3, 0, 0, -t26 * t93 - (-t63 * t21 - t106) * t59 + (t83 * t63 - t106) * qJD(5) + t73, -t28 * t93 + (t83 * qJD(5) + t21 * t59 - t11) * t66 + t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t112, t5, t3, 0, 0, t73 + t118 * (-t63 * t16 - t106), (-t118 * t16 - t11) * t66 + t70;];
tauc_reg = t4;
