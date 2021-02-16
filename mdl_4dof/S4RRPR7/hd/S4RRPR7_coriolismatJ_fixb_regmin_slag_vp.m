% Calculate minimal parameter regressor of coriolis matrix for
% S4RRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% cmat_reg [(4*%NQJ)%x21]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 10:57
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4RRPR7_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR7_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR7_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR7_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 10:56:41
% EndTime: 2021-01-15 10:56:46
% DurationCPUTime: 1.14s
% Computational Cost: add. (1041->128), mult. (2289->228), div. (0->0), fcn. (2385->6), ass. (0->128)
t88 = sin(qJ(4));
t173 = 0.2e1 * t88;
t149 = cos(pkin(7));
t171 = -qJ(3) - pkin(5);
t91 = cos(qJ(2));
t156 = t171 * t91;
t89 = sin(qJ(2));
t73 = t171 * t89;
t87 = sin(pkin(7));
t45 = -t149 * t73 - t87 * t156;
t172 = t45 / 0.2e1;
t170 = -t149 * t156 + t87 * t73;
t165 = t170 * t88;
t90 = cos(qJ(4));
t164 = t170 * t90;
t67 = -t149 * t91 + t87 * t89;
t65 = t67 ^ 2;
t109 = t149 * t89;
t161 = t87 * t91;
t69 = t109 + t161;
t66 = t69 ^ 2;
t43 = t65 + t66;
t159 = t88 * t90;
t110 = 0.2e1 * t69 * t159;
t85 = t88 ^ 2;
t86 = t90 ^ 2;
t76 = t86 - t85;
t95 = qJD(1) * t110 - qJD(2) * t76;
t167 = -t69 / 0.2e1;
t166 = t89 * pkin(2);
t163 = t45 * t90;
t36 = pkin(3) * t69 + pkin(6) * t67 + t166;
t160 = t88 * t36;
t33 = t88 * t69;
t158 = t90 * t36;
t81 = -pkin(2) * t91 - pkin(1);
t93 = t67 * pkin(3) - t69 * pkin(6) + t81;
t12 = -t90 * t93 + t165;
t8 = t12 * t67 - t45 * t33;
t155 = qJD(1) * t8;
t13 = t88 * t93 + t164;
t9 = -t13 * t67 + t69 * t163;
t154 = qJD(1) * t9;
t153 = qJD(2) * pkin(2);
t1 = (-t12 + t165) * t69 + t158 * t67;
t152 = t1 * qJD(1);
t2 = (-t13 + t164) * t69 - t160 * t67;
t151 = t2 * qJD(1);
t11 = -t170 * t67 + t45 * t69;
t148 = qJD(1) * t11;
t123 = t66 - t65;
t14 = t123 * t88;
t147 = qJD(1) * t14;
t15 = t43 * t88;
t146 = qJD(1) * t15;
t16 = t123 * t90;
t145 = qJD(1) * t16;
t23 = t67 * t166 + t69 * t81;
t144 = qJD(1) * t23;
t24 = t69 * t166 - t67 * t81;
t143 = qJD(1) * t24;
t38 = t43 * t90;
t142 = qJD(1) * t38;
t141 = qJD(1) * t91;
t140 = qJD(2) * t90;
t139 = qJD(3) * t90;
t138 = qJD(4) * t88;
t137 = qJD(4) * t90;
t10 = t81 * t166;
t136 = t10 * qJD(1);
t92 = -t87 * t67 / 0.2e1 + t149 * t167;
t19 = (-t89 / 0.2e1 + t92) * pkin(2);
t135 = t19 * qJD(1);
t30 = t88 * t67;
t134 = t30 * qJD(1);
t133 = t33 * qJD(1);
t35 = t90 * t67;
t132 = t35 * qJD(1);
t131 = t43 * qJD(1);
t64 = t109 / 0.2e1 + t161 / 0.2e1;
t130 = t64 * qJD(1);
t129 = t67 * qJD(1);
t128 = t69 * qJD(1);
t127 = t69 * qJD(2);
t77 = -t89 ^ 2 + t91 ^ 2;
t126 = t77 * qJD(1);
t125 = t89 * qJD(2);
t124 = t91 * qJD(2);
t122 = pkin(1) * t89 * qJD(1);
t121 = pkin(1) * t141;
t120 = t67 * t128;
t119 = t86 * t128;
t118 = t88 * t140;
t117 = t67 * t137;
t116 = t67 * t127;
t115 = t88 * t137;
t114 = t69 * t138;
t113 = t89 * t141;
t112 = t90 * t128;
t111 = t172 - t45 / 0.2e1;
t108 = -qJD(4) - t129;
t107 = qJD(2) * t110;
t79 = pkin(2) * t87 + pkin(6);
t80 = -t149 * pkin(2) - pkin(3);
t105 = -t67 * t80 - t69 * t79;
t104 = t108 * t90;
t103 = t79 * t67 / 0.2e1 + t80 * t167;
t94 = t36 / 0.2e1 + t103;
t4 = t111 * t90 + t94 * t88;
t102 = -qJD(1) * t4 - t80 * t140;
t6 = t111 * t88 - t94 * t90;
t101 = -qJD(2) * t80 * t88 - qJD(1) * t6;
t100 = t69 * t104;
t29 = (t85 / 0.2e1 - t86 / 0.2e1) * t69;
t99 = -qJD(1) * t29 + t118;
t98 = qJD(4) * t64 + t120;
t97 = qJD(1) * t66 * t159 + qJD(2) * t29;
t37 = t76 * t66;
t96 = qJD(1) * t37 + t107;
t63 = t64 * qJD(2);
t62 = t90 * t127;
t26 = t30 * qJD(4);
t25 = t29 * qJD(4);
t18 = t166 / 0.2e1 + t92 * pkin(2);
t17 = -t134 - t138;
t7 = t158 / 0.2e1 - t103 * t90 + t172 * t173;
t5 = t163 / 0.2e1 + t90 * t172 - t160 / 0.2e1 + t103 * t88;
t3 = [0, 0, 0, t89 * t124, t77 * qJD(2), 0, 0, 0, -pkin(1) * t125, -pkin(1) * t124, t23 * qJD(2), t24 * qJD(2), qJD(3) * t43, qJD(2) * t10 + qJD(3) * t11, -t66 * t115 - t86 * t116, -qJD(4) * t37 + t67 * t107, qJD(2) * t16 - t67 * t114, -qJD(2) * t14 - t69 * t117, t116, qJD(2) * t1 + qJD(3) * t15 + qJD(4) * t9, qJD(2) * t2 + qJD(3) * t38 + qJD(4) * t8; 0, 0, 0, t113, t126, t124, -t125, 0, -pkin(5) * t124 - t122, pkin(5) * t125 - t121, -qJD(2) * t170 + t144, qJD(2) * t45 + t143, (t149 * t67 - t69 * t87) * t153, t136 + (-t149 * t170 - t45 * t87) * t153 + t18 * qJD(3), -t25 + (-t118 - t119) * t67, -0.2e1 * t90 * t114 + t95 * t67, t88 * t127 + t145, t62 - t147, t98, t152 + (t105 * t88 - t164) * qJD(2) + t7 * qJD(4), t151 + (t105 * t90 + t165) * qJD(2) + t5 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t131, qJD(2) * t18 + t148, 0, 0, 0, 0, 0, t146, t142; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t97, -t96, t108 * t33, t100, t63, qJD(2) * t7 - qJD(4) * t13 + t154, qJD(2) * t5 + qJD(4) * t12 + t155; 0, 0, 0, -t113, -t126, 0, 0, 0, t122, t121, -qJD(3) * t69 - t144, qJD(3) * t67 - t143, 0, qJD(3) * t19 - t136, t67 * t119 - t25, t100 * t173, qJD(4) * t35 - t145, -t26 + t147, -t98, qJD(4) * t6 - t69 * t139 - t152, qJD(3) * t33 + qJD(4) * t4 - t151; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t115, t76 * qJD(4), 0, 0, 0, t80 * t138, t80 * t137; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t128, t129, 0, t135, 0, 0, 0, 0, 0, -t112, t133; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t99, -t95, t132 + t137, t17, -t130, -t79 * t137 - t101, t138 * t79 - t102; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t127, -t67 * qJD(2), -t131, -qJD(2) * t19 - t148, 0, 0, 0, 0, 0, -t26 + t62 - t146, -qJD(2) * t33 - t117 - t142; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t128, -t129, 0, -t135, 0, 0, 0, 0, 0, t112, -t133; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, t104; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t97, t96, -qJD(2) * t35 + t88 * t120, qJD(2) * t30 + t67 * t112, t63, -qJD(2) * t6 + qJD(3) * t30 - t154, -qJD(2) * t4 + t139 * t67 - t155; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t99, t95, -t132, t134, t130, t101, t102; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t134, t90 * t129; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t3;
