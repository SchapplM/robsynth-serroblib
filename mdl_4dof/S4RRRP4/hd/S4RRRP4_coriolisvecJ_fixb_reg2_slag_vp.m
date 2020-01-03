% Calculate inertial parameters regressor of coriolis joint torque vector for
% S4RRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% tauc_reg [4x(4*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RRRP4_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP4_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP4_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP4_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:15:41
% EndTime: 2019-12-31 17:15:44
% DurationCPUTime: 0.81s
% Computational Cost: add. (1192->168), mult. (3293->221), div. (0->0), fcn. (2150->4), ass. (0->112)
t139 = cos(qJ(3));
t88 = sin(qJ(3));
t89 = sin(qJ(2));
t90 = cos(qJ(2));
t64 = t139 * t89 + t88 * t90;
t122 = qJD(1) * t64;
t123 = t122 * qJ(4);
t140 = -pkin(6) - pkin(5);
t74 = t140 * t90;
t70 = qJD(1) * t74;
t53 = t88 * t70;
t126 = qJD(2) * pkin(2);
t73 = t140 * t89;
t68 = qJD(1) * t73;
t59 = t68 + t126;
t29 = t139 * t59 + t53;
t15 = t29 - t123;
t115 = t89 * t126;
t143 = 0.2e1 * t115;
t85 = qJD(2) + qJD(3);
t118 = qJD(1) * qJD(2);
t142 = -0.2e1 * t118;
t37 = -t139 * t74 + t88 * t73;
t141 = t122 ^ 2;
t111 = t139 * t90;
t101 = qJD(1) * t111;
t121 = qJD(1) * t89;
t112 = t88 * t121;
t49 = -t101 + t112;
t138 = t49 * t85;
t137 = t122 * t49;
t84 = -t90 * pkin(2) - pkin(1);
t72 = qJD(1) * t84;
t135 = t72 * t122;
t134 = t88 * t89;
t92 = qJD(1) ^ 2;
t133 = t90 * t92;
t91 = qJD(2) ^ 2;
t132 = t91 * t89;
t131 = t91 * t90;
t12 = t85 * pkin(3) + t15;
t130 = t12 - t15;
t107 = t139 * qJD(3);
t103 = pkin(2) * t107;
t34 = t85 * t64;
t26 = t34 * qJD(1);
t129 = -t88 * pkin(2) * t26 - t49 * t103;
t32 = t139 * t68 + t53;
t128 = t85 * t101;
t127 = t89 ^ 2 - t90 ^ 2;
t99 = t85 * t134;
t25 = qJD(1) * t99 - t128;
t125 = t25 * qJ(4);
t124 = t49 * qJ(4);
t120 = qJD(3) * t88;
t106 = t49 * pkin(3) + qJD(4);
t35 = t106 + t72;
t119 = qJD(4) + t35;
t117 = t89 * t133;
t116 = pkin(2) * t120;
t114 = pkin(2) * t121;
t113 = t139 * pkin(2);
t57 = t139 * t70;
t109 = t89 * t118;
t20 = pkin(2) * t109 + t26 * pkin(3);
t110 = qJD(2) * t140;
t100 = qJD(1) * t110;
t60 = t89 * t100;
t61 = t90 * t100;
t108 = t139 * t61 - t88 * t60;
t31 = -t88 * t68 + t57;
t36 = t139 * t73 + t88 * t74;
t105 = -t59 * t107 - t70 * t120 - t139 * t60 - t88 * t61;
t104 = pkin(1) * t142;
t102 = t90 * t109;
t98 = t72 * t49 + t105;
t30 = t88 * t59 - t57;
t97 = t26 * qJ(4) + t105;
t69 = t89 * t110;
t71 = t90 * t110;
t10 = t73 * t107 + t74 * t120 + t139 * t69 + t88 * t71;
t96 = -t85 * t112 + t128;
t95 = t119 * t49 + t97;
t9 = -qJD(3) * t30 + t108;
t11 = -t37 * qJD(3) + t139 * t71 - t88 * t69;
t94 = t9 + t125;
t93 = (t57 + (-pkin(2) * t85 - t59) * t88) * qJD(3) + t108;
t83 = t113 + pkin(3);
t75 = t85 * t103;
t63 = -t111 + t134;
t48 = t49 ^ 2;
t40 = t63 * pkin(3) + t84;
t38 = pkin(3) * t122 + t114;
t33 = -qJD(2) * t111 - t90 * t107 + t99;
t28 = t34 * t85;
t27 = t33 * t85;
t24 = t34 * pkin(3) + t115;
t22 = -t63 * qJ(4) + t37;
t21 = -t64 * qJ(4) + t36;
t19 = -t48 + t141;
t18 = -t123 + t32;
t17 = t31 + t124;
t16 = t30 - t124;
t13 = t96 + t138;
t7 = t26 * t63 + t49 * t34;
t6 = -t122 * t33 - t25 * t64;
t5 = t33 * qJ(4) - t64 * qJD(4) + t11;
t4 = -t34 * qJ(4) - t63 * qJD(4) + t10;
t3 = -qJD(4) * t122 + t94;
t2 = -t49 * qJD(4) - t97;
t1 = -t122 * t34 + t25 * t63 - t64 * t26 + t33 * t49;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t102, t127 * t142, t131, -0.2e1 * t102, -t132, 0, -pkin(5) * t131 + t89 * t104, pkin(5) * t132 + t90 * t104, 0, 0, t6, t1, -t27, t7, -t28, 0, t11 * t85 + t84 * t26 + t72 * t34 + (qJD(1) * t63 + t49) * t115, -t10 * t85 + t122 * t143 - t84 * t25 - t72 * t33, -t10 * t49 + t105 * t63 - t11 * t122 + t36 * t25 - t37 * t26 + t29 * t33 - t30 * t34 - t9 * t64, t30 * t10 - t105 * t37 + t29 * t11 + t72 * t143 + t9 * t36, t6, t1, -t27, t7, -t28, 0, t20 * t63 + t24 * t49 + t40 * t26 + t35 * t34 + t5 * t85, t122 * t24 + t20 * t64 - t40 * t25 - t35 * t33 - t4 * t85, t12 * t33 - t122 * t5 - t16 * t34 - t2 * t63 + t21 * t25 - t22 * t26 - t3 * t64 - t4 * t49, t12 * t5 + t16 * t4 + t2 * t22 + t20 * t40 + t3 * t21 + t35 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t117, t127 * t92, 0, t117, 0, 0, t92 * pkin(1) * t89, pkin(1) * t133, 0, 0, t137, t19, t13, -t137, 0, 0, -t49 * t114 - t31 * t85 - t135 + t93, -t114 * t122 + t32 * t85 - t75 + t98, t25 * t113 + (-t29 + t32) * t49 + (t30 + t31 + t116) * t122 + t129, -t29 * t31 - t30 * t32 + (-t72 * t121 + t139 * t9 - t105 * t88 + (t139 * t30 - t29 * t88) * qJD(3)) * pkin(2), t137, t19, t13, -t137, 0, 0, -t119 * t122 - t17 * t85 - t38 * t49 + t125 + t93, -t122 * t38 + t18 * t85 - t75 + t95, t83 * t25 + (-t12 + t18) * t49 + (t16 + t17 + t116) * t122 + t129, -t12 * t17 - t16 * t18 + t3 * t83 - t35 * t38 + (t2 * t88 + (-t12 * t88 + t139 * t16) * qJD(3)) * pkin(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t137, t19, t13, -t137, 0, 0, t30 * t85 - t135 + t9, t29 * t85 + t98, 0, 0, t137, t19, t13, -t137, 0, 0, t16 * t85 + (-t106 - t35) * t122 + t94, -t141 * pkin(3) + t15 * t85 + t95, t25 * pkin(3) - t130 * t49, t130 * t16 + (-t122 * t35 + t3) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t85 * t122 + t26, t96 - t138, -t48 - t141, t12 * t122 + t16 * t49 + t20;];
tauc_reg = t8;
