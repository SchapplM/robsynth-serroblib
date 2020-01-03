% Calculate minimal parameter regressor of coriolis matrix for
% S5RRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x23]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RRPPR5_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR5_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR5_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR5_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:30:00
% EndTime: 2019-12-31 19:30:04
% DurationCPUTime: 1.52s
% Computational Cost: add. (1439->132), mult. (2816->186), div. (0->0), fcn. (3118->6), ass. (0->108)
t122 = qJD(2) - qJD(5);
t156 = cos(pkin(8));
t166 = -qJ(3) - pkin(6);
t96 = cos(qJ(2));
t86 = t166 * t96;
t108 = t156 * t86;
t94 = sin(qJ(2));
t85 = t166 * t94;
t92 = sin(pkin(8));
t176 = t92 * t85;
t186 = -t108 + t176;
t74 = -t156 * t96 + t92 * t94;
t187 = pkin(7) * t74 + t186;
t54 = t156 * t85 + t92 * t86;
t76 = t156 * t94 + t92 * t96;
t36 = pkin(7) * t76 + t54;
t93 = sin(qJ(5));
t95 = cos(qJ(5));
t209 = t122 * (t93 * t187 + t95 * t36);
t208 = t122 * (t95 * t187 - t93 * t36);
t103 = t93 * t74 + t95 * t76;
t190 = -t95 * t74 + t93 * t76;
t10 = t103 ^ 2 - t190 ^ 2;
t202 = t10 * qJD(1);
t145 = qJD(5) * t190;
t199 = qJD(2) * t190 - t145;
t198 = qJD(3) * t103;
t197 = qJD(3) * t190;
t196 = t103 * qJD(1);
t195 = t190 * qJD(1);
t73 = t76 ^ 2;
t191 = t74 ^ 2 + t73;
t194 = t191 * qJD(1);
t144 = qJD(5) * t103;
t193 = qJD(2) * t103 - t144;
t192 = qJD(3) * t191;
t181 = -t108 / 0.2e1;
t104 = -t186 * t74 - t54 * t76;
t185 = qJD(1) * t104;
t184 = qJD(3) * t104;
t182 = -t74 / 0.2e1;
t180 = -pkin(3) - pkin(4);
t178 = t76 * pkin(3);
t177 = t94 * pkin(2);
t165 = qJD(2) * pkin(2);
t109 = -t96 * pkin(2) - pkin(1);
t99 = t76 * qJ(4) - t109;
t40 = t74 * pkin(3) - t99;
t159 = t74 * qJ(4);
t106 = -t159 - t177;
t43 = -t106 + t178;
t5 = t40 * t43;
t162 = t5 * qJD(1);
t22 = t180 * t74 + t99;
t23 = t180 * t76 + t106;
t6 = -t103 * t22 + t190 * t23;
t161 = t6 * qJD(1);
t7 = t103 * t23 + t190 * t22;
t160 = t7 * qJD(1);
t13 = t40 * t76 + t43 * t74;
t155 = qJD(1) * t13;
t14 = t40 * t74 - t43 * t76;
t154 = qJD(1) * t14;
t150 = qJD(1) * t96;
t148 = qJD(4) * t76;
t147 = qJD(4) * t93;
t146 = qJD(4) * t95;
t11 = t109 * t177;
t142 = t11 * qJD(1);
t88 = t92 * pkin(2) + qJ(4);
t90 = -t156 * pkin(2) - pkin(3);
t91 = t177 / 0.2e1;
t18 = t91 + (pkin(3) / 0.2e1 - t90 / 0.2e1) * t76 + (qJ(4) / 0.2e1 + t88 / 0.2e1) * t74;
t140 = t18 * qJD(1);
t98 = t92 * t182 - t156 * t76 / 0.2e1;
t33 = (-t94 / 0.2e1 + t98) * pkin(2);
t137 = t33 * qJD(1);
t131 = t73 * qJD(1);
t130 = t74 * qJD(1);
t129 = t74 * qJD(2);
t128 = t76 * qJD(1);
t87 = -t94 ^ 2 + t96 ^ 2;
t127 = t87 * qJD(1);
t126 = t93 * qJD(2);
t125 = t94 * qJD(2);
t124 = t95 * qJD(2);
t123 = t96 * qJD(2);
t121 = pkin(1) * t94 * qJD(1);
t120 = pkin(1) * t150;
t119 = t22 * t195;
t118 = t22 * t196;
t117 = t190 * t196;
t116 = t103 * t195;
t115 = t190 * t128;
t114 = t103 * t128;
t113 = t74 * t128;
t112 = t94 * t150;
t102 = -pkin(4) + t90;
t51 = t181 + t108 / 0.2e1;
t101 = qJD(1) * t51 + qJD(2) * t88;
t84 = t122 * t95;
t83 = t122 * t93;
t59 = t93 * t102 + t95 * t88;
t58 = -t95 * t102 + t93 * t88;
t38 = 0.2e1 * t181 + t176;
t32 = t98 * pkin(2) + t91;
t19 = t88 * t182 + t90 * t76 / 0.2e1 + t91 + t159 / 0.2e1 + t178 / 0.2e1;
t1 = [0, 0, 0, t94 * t123, t87 * qJD(2), 0, 0, 0, -pkin(1) * t125, -pkin(1) * t123, t192, qJD(2) * t11 + t184, qJD(2) * t13 - t74 * t148, t192, qJD(2) * t14 + qJD(4) * t73, qJD(2) * t5 - t148 * t40 + t184, t199 * t103, t122 * t10, 0, 0, 0, qJD(2) * t6 + t144 * t22 + t148 * t190, qJD(2) * t7 + t103 * t148 - t145 * t22; 0, 0, 0, t112, t127, t123, -t125, 0, -pkin(6) * t123 - t121, pkin(6) * t125 - t120, (t156 * t74 - t76 * t92) * t165, t142 + (-t156 * t186 + t54 * t92) * t165 + t32 * qJD(3), -qJD(2) * t186 + t155, (-t74 * t90 - t76 * t88) * qJD(2) - qJD(4) * t74, qJD(2) * t54 + t154, t162 + (t186 * t90 + t54 * t88) * qJD(2) + t19 * qJD(3) + t38 * qJD(4), t116, t202, -t199, -t193, 0, t161 - t208, t160 + t209; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t194, qJD(2) * t32 + t185, 0, t194, 0, qJD(2) * t19 + t185, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t113, -t129, t131, qJD(2) * t38 - t128 * t40, 0, 0, 0, 0, 0, t115, t114; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t117, -t202, t199, t193, 0, t118 + t208, -t119 - t209; 0, 0, 0, -t112, -t127, 0, 0, 0, t121, t120, 0, qJD(3) * t33 - t142, -qJD(3) * t76 - t155, 0, -qJD(3) * t74 - t154, -qJD(3) * t18 + qJD(4) * t51 - t162, -t116, -t202, 0, 0, 0, -t161 - t198, -t160 + t197; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4), t88 * qJD(4), 0, 0, 0, 0, 0, qJD(5) * t59 + t147, -qJD(5) * t58 + t146; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t137, -t128, 0, -t130, -t140, 0, 0, 0, 0, 0, -t196, t195; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2), t101, 0, 0, 0, 0, 0, t126, t124; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t122 * t59, -t122 * t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t194, -qJD(2) * t33 - t185, t76 * qJD(2), -t194, t129, qJD(2) * t18 - t148 - t185, 0, 0, 0, 0, 0, t193, -t199; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t137, t128, 0, t130, t140, 0, 0, 0, 0, 0, t196, -t195; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t128, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t196, t195; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t113, 0, -t131, -t51 * qJD(2) + (qJD(1) * t40 + qJD(3)) * t76, 0, 0, 0, 0, 0, -t115, -t114; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(2), -t101, 0, 0, 0, 0, 0, -t83, -t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t128, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t83, t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t117, t202, 0, 0, 0, -t118 + t198, t119 - t197; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(2) * t59 - t147, qJD(2) * t58 - t146; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t196, -t195; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t126, -t124; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t1;
