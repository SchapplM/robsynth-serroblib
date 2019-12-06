% Calculate minimal parameter regressor of coriolis matrix for
% S5RPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x18]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RPRPR3_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR3_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:51:35
% EndTime: 2019-12-05 17:51:39
% DurationCPUTime: 0.70s
% Computational Cost: add. (830->117), mult. (1846->166), div. (0->0), fcn. (1474->8), ass. (0->117)
t80 = sin(pkin(9));
t78 = t80 ^ 2;
t82 = cos(pkin(9));
t79 = t82 ^ 2;
t73 = t78 + t79;
t112 = qJD(1) + qJD(3);
t83 = sin(qJ(5));
t85 = cos(qJ(5));
t60 = (t83 ^ 2 - t85 ^ 2) * t78;
t147 = t112 * t60;
t61 = t73 * t83;
t27 = t112 * t61;
t62 = t73 * t85;
t28 = t112 * t62;
t146 = t112 * t73;
t99 = t112 * t82;
t84 = sin(qJ(3));
t96 = cos(pkin(8)) * pkin(1) + pkin(2);
t89 = t84 * t96;
t140 = cos(qJ(3));
t142 = pkin(1) * sin(pkin(8));
t98 = t140 * t142;
t56 = t98 + t89;
t54 = qJ(4) + t56;
t145 = qJ(4) + t54;
t144 = t54 / 0.2e1;
t55 = -t140 * t96 + t84 * t142;
t143 = t55 / 0.2e1;
t141 = qJ(4) / 0.2e1;
t133 = t82 * t83;
t64 = -t82 * pkin(4) - t80 * pkin(7) - pkin(3);
t31 = t55 + t64;
t16 = t54 * t133 - t85 * t31;
t139 = t16 * t82;
t132 = t82 * t85;
t17 = -t54 * t132 - t83 * t31;
t138 = t17 * t82;
t126 = qJ(4) * t82;
t43 = t83 * t126 - t85 * t64;
t137 = t43 * t82;
t44 = -t85 * t126 - t83 * t64;
t136 = t44 * t82;
t135 = t78 * t83;
t134 = t78 * t85;
t131 = t83 * t56;
t130 = t85 * t56;
t13 = (-t55 * t132 + t131) * t82 - t55 * t134;
t59 = t62 * qJD(4);
t129 = t13 * qJD(3) + t59;
t20 = t73 * t55;
t68 = t73 * qJD(4);
t128 = -t20 * qJD(3) + t68;
t12 = (t55 * t133 + t130) * t82 + t55 * t135;
t58 = t61 * qJD(4);
t127 = -t12 * qJD(3) + t58;
t5 = (-pkin(3) + t55) * t56 - t54 * t20;
t125 = t5 * qJD(1);
t6 = -t54 * t135 - t139;
t124 = t6 * qJD(1);
t38 = t54 * t134;
t7 = -t38 + t138;
t123 = t7 * qJD(1);
t122 = t78 * qJ(4);
t121 = qJD(4) * t82;
t120 = qJD(5) * t83;
t119 = qJD(5) * t85;
t118 = t12 * qJD(1);
t117 = t13 * qJD(1);
t18 = t73 * t54;
t116 = t18 * qJD(1);
t115 = t20 * qJD(1);
t114 = t55 * qJD(1);
t113 = t56 * qJD(1);
t53 = t56 * qJD(3);
t111 = t83 * t134;
t110 = t82 * t120;
t109 = t82 * t119;
t108 = t80 * t113;
t107 = t83 * t121;
t106 = t85 * t121;
t105 = t83 * t143;
t104 = t85 * t143;
t102 = t130 / 0.2e1;
t74 = t85 * t122;
t101 = -t38 / 0.2e1 - t74 / 0.2e1;
t97 = qJD(5) * t111;
t95 = t83 * t99;
t94 = t85 * t99;
t1 = t102 + (t105 + t44 / 0.2e1 + t17 / 0.2e1) * t82 + t101;
t24 = -t74 + t136;
t93 = -t1 * qJD(1) - t24 * qJD(3);
t2 = (-t56 / 0.2e1 + (t141 + t144) * t78) * t83 + (t104 + t43 / 0.2e1 + t16 / 0.2e1) * t82;
t23 = -t83 * t122 - t137;
t92 = -t2 * qJD(1) + t23 * qJD(3);
t67 = t73 * qJ(4);
t86 = t98 / 0.2e1 + t89 / 0.2e1;
t8 = t86 + t145 * (-t79 / 0.2e1 - t78 / 0.2e1);
t91 = t8 * qJD(1) - t67 * qJD(3);
t90 = -qJD(5) + t99;
t88 = t90 * t83;
t87 = t90 * t85;
t66 = t80 * t109;
t65 = t80 * t110;
t57 = t60 * qJD(5);
t52 = t55 * qJD(3);
t49 = t112 * t111;
t46 = t80 * t94;
t45 = t80 * t95;
t42 = t80 * t53;
t30 = t80 * t87;
t29 = t80 * t88;
t22 = t109 - t28;
t21 = t110 - t27;
t9 = t79 * t141 + t122 / 0.2e1 + t86 + t73 * t144;
t4 = -t136 / 0.2e1 - t138 / 0.2e1 + t82 * t105 + t102 - t101;
t3 = -t137 / 0.2e1 - t139 / 0.2e1 + t82 * t104 - t131 / 0.2e1 - t145 * t135 / 0.2e1;
t10 = [0, 0, 0, 0, 0, -t53, t52, -t82 * t53, t42, t128, t5 * qJD(3) + t18 * qJD(4), -t97, t57, t65, t66, 0, -t7 * qJD(5) + t127, t6 * qJD(5) + t129; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -t53 - t113, t52 + t114, -t56 * t99, t42 + t108, -t115 + t128, t125 + (-t56 * pkin(3) - qJ(4) * t20) * qJD(3) + t9 * qJD(4), -t97, t57, t65, t66, 0, t4 * qJD(5) - t118 + t127, t3 * qJD(5) + t117 + t129; 0, 0, 0, 0, 0, 0, 0, 0, 0, t146, t9 * qJD(3) + t116, 0, 0, 0, 0, 0, t27, t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49, t147, t29, t30, 0, t4 * qJD(3) + t17 * qJD(5) - t123, t3 * qJD(3) + t16 * qJD(5) + t124; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t80 * t119, t80 * t120; 0, 0, 0, 0, 0, t113, -t114, t82 * t113, -t108, t68 + t115, -t8 * qJD(4) - t125, -t97, t57, t65, t66, 0, -t1 * qJD(5) + t118 + t58, -t2 * qJD(5) - t117 + t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, t68, t67 * qJD(4), -t97, t57, t65, t66, 0, -t24 * qJD(5) + t58, t23 * qJD(5) + t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, t146, -t91, 0, 0, 0, 0, 0, t27, t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49, t147, t29, t30, 0, t44 * qJD(5) + t93, t43 * qJD(5) + t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t146, t8 * qJD(3) - t116, 0, 0, 0, 0, 0, t21, t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t146, t91, 0, 0, 0, 0, 0, t21, t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t88, t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, -t147, -t45, -t46, 0, t1 * qJD(3) - t107 + t123, t2 * qJD(3) - t106 - t124; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, -t147, -t45, -t46, 0, -t93 - t107, -t92 - t106; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t95, -t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t10;
