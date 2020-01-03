% Calculate minimal parameter regressor of coriolis matrix for
% S4RRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
% 
% Output:
% cmat_reg [(4*%NQJ)%x21]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4RRPR10_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR10_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR10_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR10_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:11:57
% EndTime: 2019-12-31 17:11:59
% DurationCPUTime: 0.64s
% Computational Cost: add. (364->109), mult. (861->180), div. (0->0), fcn. (685->4), ass. (0->107)
t54 = sin(qJ(2));
t107 = t54 * qJ(3);
t56 = cos(qJ(2));
t57 = -pkin(2) - pkin(6);
t69 = -t56 * t57 + t107;
t53 = sin(qJ(4));
t49 = t53 ^ 2;
t55 = cos(qJ(4));
t51 = t55 ^ 2;
t36 = t49 - t51;
t88 = t56 * qJD(1);
t81 = t55 * t88;
t60 = t36 * qJD(2) + 0.2e1 * t53 * t81;
t87 = t56 * qJD(3);
t120 = t69 * qJD(2) - t87;
t119 = pkin(3) + pkin(5);
t106 = t56 * qJ(3);
t31 = t54 * pkin(2) - t106;
t19 = t54 * pkin(6) + t31;
t117 = t53 * t19;
t32 = t119 * t54;
t116 = t53 * t32;
t115 = t53 * t56;
t114 = t55 * t19;
t113 = t55 * t32;
t112 = t55 * t56;
t50 = t54 ^ 2;
t52 = t56 ^ 2;
t37 = t52 - t50;
t17 = -pkin(1) - t69;
t9 = t53 * t17 - t113;
t1 = (-t9 - t113) * t56 - t117 * t54;
t111 = t1 * qJD(1);
t10 = t55 * t17 + t116;
t2 = t114 * t54 + (t10 - t116) * t56;
t110 = t2 * qJD(1);
t33 = t119 * t56;
t3 = t33 * t112 - t9 * t54;
t109 = t3 * qJD(1);
t4 = -t10 * t54 - t33 * t115;
t108 = t4 * qJD(1);
t105 = qJD(1) * t54;
t104 = qJD(1) * t55;
t103 = qJD(3) * t54;
t102 = qJD(3) * t55;
t101 = qJD(4) * t54;
t100 = qJD(4) * t55;
t99 = qJD(4) * t57;
t70 = -t56 * pkin(2) - t107;
t28 = -pkin(1) + t70;
t11 = t28 * t56 + t31 * t54;
t98 = t11 * qJD(1);
t12 = -t28 * t54 + t31 * t56;
t97 = t12 * qJD(1);
t23 = t37 * t53;
t96 = t23 * qJD(1);
t25 = t37 * t55;
t95 = t25 * qJD(1);
t94 = t37 * qJD(1);
t93 = t50 * qJD(1);
t92 = t50 * qJD(3);
t91 = t53 * qJD(2);
t90 = t54 * qJD(2);
t89 = t55 * qJD(2);
t47 = t56 * qJD(2);
t86 = qJ(3) * qJD(4);
t85 = qJD(2) * qJ(3);
t84 = pkin(1) * t105;
t83 = pkin(1) * t88;
t82 = pkin(5) * t90;
t80 = t53 * t101;
t79 = t54 * t100;
t78 = t28 * t31 * qJD(1);
t77 = t28 * t105;
t76 = t53 * t47;
t42 = t54 * t47;
t41 = t54 * t88;
t75 = t53 * t100;
t74 = t53 * t89;
t72 = t56 * t74;
t71 = qJD(4) + t105;
t68 = -t93 - t101;
t67 = -t57 * t54 / 0.2e1 - t106 / 0.2e1;
t66 = t71 * t115;
t59 = t19 / 0.2e1 + t67;
t8 = t59 * t55;
t65 = -t8 * qJD(1) + t53 * t85;
t7 = t59 * t53;
t64 = -t7 * qJD(1) - t55 * t85;
t18 = (t51 / 0.2e1 - t49 / 0.2e1) * t56;
t63 = t18 * qJD(1) + t74;
t62 = t53 * t52 * t104 - t18 * qJD(2);
t24 = t36 * t52;
t61 = -t24 * qJD(1) + 0.2e1 * t72;
t58 = t70 * qJD(2) + t87;
t46 = pkin(5) * t47;
t45 = t47 / 0.2e1;
t44 = t55 * t47;
t43 = t54 * t104;
t40 = t53 * t105;
t22 = -t43 - t100;
t21 = -qJD(4) * t53 - t40;
t20 = t41 + t56 * qJD(4) / 0.2e1;
t15 = t18 * qJD(4);
t6 = -t33 * t53 - t114 / 0.2e1 + t67 * t55;
t5 = t33 * t55 - t117 / 0.2e1 + t67 * t53;
t13 = [0, 0, 0, t42, t37 * qJD(2), 0, 0, 0, -pkin(1) * t90, -pkin(1) * t47, 0, t12 * qJD(2) - t54 * t87, -t11 * qJD(2) + t92, (qJD(2) * t31 - t103) * t28, -t49 * t42 + t52 * t75, -t24 * qJD(4) - 0.2e1 * t54 * t72, -t23 * qJD(2) - t56 * t79, -t25 * qJD(2) + t56 * t80, t42, t1 * qJD(2) + t4 * qJD(4) + t53 * t92, -t2 * qJD(2) - t3 * qJD(4) + t55 * t92; 0, 0, 0, t41, t94, t47, -t90, 0, -t46 - t84, t82 - t83, t58, t46 + t97, -t82 - t98, t58 * pkin(5) + t78, -t15 + (-t49 * t88 + t74) * t54, -t60 * t54 + 0.2e1 * t56 * t75, t44 - t96, -t76 - t95, t20, t5 * qJD(4) - t120 * t55 - t32 * t91 + t111, t6 * qJD(4) + t120 * t53 - t32 * t89 - t110; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, -t41, t93, t46 - t77, 0, 0, 0, 0, 0, t53 * t93 + t44, t55 * t93 - t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, t61, -t71 * t112, t66, t45, t5 * qJD(2) - t10 * qJD(4) + t108, t6 * qJD(2) + t9 * qJD(4) - t109; 0, 0, 0, -t41, -t94, 0, 0, 0, t84, t83, 0, -t97, t98, -t78, t49 * t41 - t15, 0.2e1 * t55 * t66, -t80 + t96, -t79 + t95, -t20, t7 * qJD(4) - t111, t8 * qJD(4) + t110; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3), qJ(3) * qJD(3), -t75, t36 * qJD(4), 0, 0, 0, qJD(3) * t53 + t55 * t86, -t53 * t86 + t102; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2), t85, 0, 0, 0, 0, 0, t91, t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t63, t60, t21, t22, -t88 / 0.2e1, -t53 * t99 - t64, -t55 * t99 - t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, -t93, t77, 0, 0, 0, 0, 0, t68 * t53, t68 * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(2), -t85, 0, 0, 0, 0, 0, -t91, -t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t62, -t61, (t81 + t91) * t54, (-t53 * t88 + t89) * t54, t45, -t7 * qJD(2) + t53 * t103 - t108, -t8 * qJD(2) + t54 * t102 + t109; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63, -t60, t40, t43, t88 / 0.2e1, t64, t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t13;
