% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRRP5_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP5_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP5_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP5_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:40:57
% EndTime: 2019-12-31 18:40:59
% DurationCPUTime: 0.55s
% Computational Cost: add. (1091->132), mult. (2365->154), div. (0->0), fcn. (1254->6), ass. (0->91)
t62 = sin(qJ(4));
t58 = t62 ^ 2;
t64 = cos(qJ(4));
t59 = t64 ^ 2;
t121 = t58 + t59;
t51 = cos(pkin(8)) * pkin(1) + pkin(2);
t45 = t51 * qJD(1);
t63 = sin(qJ(3));
t65 = cos(qJ(3));
t115 = sin(pkin(8)) * pkin(1);
t88 = qJD(1) * t115;
t29 = t65 * t45 - t63 * t88;
t30 = t63 * t45 + t65 * t88;
t57 = qJD(1) + qJD(3);
t22 = t57 * pkin(7) + t30;
t104 = t62 * t22;
t90 = t64 * qJD(2);
t14 = t90 - t104;
t120 = qJD(5) - t14;
t27 = t29 * qJD(3);
t98 = qJD(4) * t90 + t64 * t27;
t2 = (qJD(5) - t104) * qJD(4) + t98;
t103 = t62 * t27;
t15 = t62 * qJD(2) + t64 * t22;
t6 = t15 * qJD(4) + t103;
t4 = t6 * t62;
t119 = t2 * t64 + t4;
t93 = qJD(4) * t62;
t5 = -t22 * t93 + t98;
t118 = t5 * t64 + t4 + (-t14 * t64 - t15 * t62) * qJD(4);
t11 = -qJD(4) * pkin(4) + t120;
t89 = qJD(4) * qJ(5);
t12 = t15 + t89;
t95 = t65 * t115 + t63 * t51;
t83 = -t63 * t115 + t65 * t51;
t117 = t30 * qJD(3);
t94 = -t58 + t59;
t116 = 0.2e1 * t94 * t57 * qJD(4);
t66 = qJD(4) ^ 2;
t114 = pkin(7) * t66;
t113 = t57 * pkin(3);
t112 = t6 * t64;
t111 = t29 * t57;
t110 = t30 * t57;
t31 = t83 * qJD(3);
t109 = t31 * t57;
t32 = t95 * qJD(3);
t108 = t32 * t57;
t36 = pkin(7) + t95;
t107 = t36 * t66;
t43 = -t64 * pkin(4) - t62 * qJ(5) - pkin(3);
t106 = t43 * t57;
t105 = t57 * t64;
t21 = -t29 - t113;
t92 = qJD(4) * t64;
t100 = t117 * t62 + t21 * t92;
t99 = t30 * t105 + t29 * t93;
t97 = t121 * t109;
t91 = t62 * qJD(5);
t38 = pkin(4) * t93 - t64 * t89 - t91;
t96 = t30 - t38;
t87 = t57 * t93;
t79 = pkin(4) * t62 - qJ(5) * t64;
t7 = t117 + (t79 * qJD(4) - t91) * t57;
t86 = -t7 - t114;
t85 = t14 + t104;
t23 = t43 - t83;
t84 = t23 * t57 - t31;
t81 = t64 * t87;
t80 = t121 * t111;
t78 = t11 * t62 + t12 * t64;
t77 = t14 * t62 - t15 * t64;
t76 = t107 + t108;
t35 = -pkin(3) - t83;
t75 = qJD(4) * (t35 * t57 - t31);
t13 = t32 + t38;
t74 = -t13 * t57 - t107 - t7;
t73 = t11 * t92 - t12 * t93 + t119;
t68 = (t11 * t64 - t12 * t62) * qJD(4) + t119;
t56 = t57 ^ 2;
t55 = t66 * t64;
t54 = t66 * t62;
t46 = t62 * t56 * t64;
t42 = -0.2e1 * t81;
t41 = 0.2e1 * t81;
t40 = t94 * t56;
t39 = t79 * t57;
t16 = t21 * t93;
t10 = -t29 + t106;
t8 = t10 * t93;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t117 - t108, -t27 - t109, 0, -t117 * t83 + t27 * t95 - t29 * t32 + t30 * t31, t41, t116, t55, t42, -t54, 0, t16 + t62 * t75 + (-t117 - t76) * t64, t62 * t76 + t64 * t75 + t100, t118 + t97, t117 * t35 + t118 * t36 + t21 * t32 - t31 * t77, t41, t55, -t116, 0, t54, t42, t64 * t74 + t84 * t93 + t8, t73 + t97, t74 * t62 + (-t10 - t84) * t92, t10 * t13 + t7 * t23 + t31 * t78 + t36 * t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t54, -t55, 0, -qJD(4) * t77 + t5 * t62 - t112, 0, 0, 0, 0, 0, 0, -t54, 0, t55, qJD(4) * t78 + t2 * t62 - t112; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t117 + t110, -t27 + t111, 0, 0, t41, t116, t55, t42, -t54, 0, -pkin(3) * t87 + t16 + (-t117 - t114) * t64 + t99, (-t110 + t114) * t62 + (t29 - t113) * t92 + t100, -t80 + t118, -pkin(3) * t117 + pkin(7) * t118 - t21 * t30 + t29 * t77, t41, t55, -t116, 0, t54, t42, t43 * t87 + t8 + (-t38 * t57 + t86) * t64 + t99, -t80 + t73, (-t10 - t29 - t106) * t92 + (t96 * t57 + t86) * t62, t68 * pkin(7) - t96 * t10 - t78 * t29 + t7 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46, -t40, 0, t46, 0, 0, (-t21 * t57 - t27) * t62, qJD(4) * t85 - t21 * t105 - t98, 0, 0, -t46, 0, t40, 0, 0, t46, -t103 + (-t10 * t62 + t39 * t64) * t57, 0, (t10 * t64 + t39 * t62) * t57 + (0.2e1 * qJD(5) - t85) * qJD(4) + t98, -t6 * pkin(4) + t2 * qJ(5) - t10 * t39 - t11 * t15 + t120 * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46, 0, -t58 * t56 - t66, (t10 * t57 + t27) * t62 + (-t12 + t15) * qJD(4);];
tauc_reg = t1;
