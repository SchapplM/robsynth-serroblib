% Calculate minimal parameter regressor of coriolis matrix for
% S5RPRPR7
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
% cmat_reg [(5*%NQJ)%x20]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:20
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RPRPR7_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR7_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR7_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR7_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:19:53
% EndTime: 2019-12-31 18:19:56
% DurationCPUTime: 1.14s
% Computational Cost: add. (1635->124), mult. (3320->216), div. (0->0), fcn. (3502->8), ass. (0->130)
t117 = cos(qJ(5));
t113 = sin(pkin(9));
t116 = sin(qJ(3));
t107 = sin(pkin(8)) * pkin(1) + pkin(6);
t170 = qJ(4) + t107;
t142 = t170 * t116;
t188 = cos(pkin(9));
t118 = cos(qJ(3));
t94 = t170 * t118;
t41 = t113 * t94 + t188 * t142;
t189 = t41 * t117;
t115 = sin(qJ(5));
t97 = t113 * t116 - t188 * t118;
t54 = t115 * t97;
t47 = t54 * qJD(5);
t166 = qJD(3) * t117;
t144 = t188 * t116;
t186 = t113 * t118;
t99 = t144 + t186;
t89 = t99 * t166;
t200 = -t89 + t47;
t95 = t97 ^ 2;
t96 = t99 ^ 2;
t66 = t95 + t96;
t120 = -t113 * t142 + t188 * t94;
t199 = t115 * t120;
t198 = t117 * t120;
t111 = t115 ^ 2;
t112 = t117 ^ 2;
t103 = t112 - t111;
t56 = t115 * t99;
t143 = 0.2e1 * t117 * t56;
t123 = qJD(1) * t143 - t103 * qJD(3);
t195 = t116 * pkin(3);
t59 = t99 * pkin(4) + t97 * pkin(7) + t195;
t197 = t59 / 0.2e1;
t196 = -t115 / 0.2e1;
t194 = qJD(3) * pkin(3);
t131 = t120 * t99;
t109 = -cos(pkin(8)) * pkin(1) - pkin(2);
t133 = -t118 * pkin(3) + t109;
t119 = t97 * pkin(4) - t99 * pkin(7) + t133;
t18 = -t117 * t119 + t199;
t58 = t117 * t97;
t1 = t131 * t115 - t18 * t99 + t59 * t58;
t193 = t1 * qJD(1);
t192 = t113 * t97;
t190 = t41 * t115;
t187 = qJD(1) * t99;
t12 = t18 * t97 - t99 * t190;
t185 = t12 * qJD(1);
t19 = t115 * t119 + t198;
t13 = t99 * t189 - t19 * t97;
t184 = t13 * qJD(1);
t14 = -t120 * t97 + t41 * t99;
t183 = t14 * qJD(1);
t159 = t96 - t95;
t27 = t159 * t115;
t180 = t27 * qJD(1);
t28 = t66 * t115;
t179 = t28 * qJD(1);
t29 = t159 * t117;
t178 = t29 * qJD(1);
t145 = t188 * t99;
t122 = -t192 / 0.2e1 - t145 / 0.2e1;
t33 = (-t116 / 0.2e1 + t122) * pkin(3);
t177 = t33 * qJD(1);
t46 = t54 * qJD(1);
t176 = t56 * qJD(1);
t175 = t58 * qJD(1);
t174 = t58 * qJD(3);
t61 = t66 * t117;
t173 = t61 * qJD(1);
t172 = t66 * qJD(1);
t93 = t144 / 0.2e1 + t186 / 0.2e1;
t171 = t93 * qJD(1);
t169 = qJD(1) * t117;
t168 = qJD(1) * t118;
t167 = qJD(3) * t115;
t165 = qJD(4) * t117;
t164 = qJD(5) * t115;
t163 = qJD(5) * t117;
t104 = -t116 ^ 2 + t118 ^ 2;
t162 = t104 * qJD(1);
t161 = t116 * qJD(3);
t160 = t118 * qJD(3);
t158 = t97 * t187;
t157 = t97 * t99 * qJD(3);
t155 = t112 * t187;
t154 = t99 * t164;
t153 = t99 * t163;
t152 = t99 * t169;
t150 = t109 * t116 * qJD(1);
t149 = t109 * t168;
t148 = t115 * t163;
t147 = t115 * t166;
t146 = t116 * t168;
t141 = -qJD(1) * t97 - qJD(5);
t139 = qJD(3) * t143;
t106 = t113 * pkin(3) + pkin(7);
t108 = -t188 * pkin(3) - pkin(4);
t138 = -t106 * t99 - t108 * t97;
t10 = t133 * t195;
t137 = t10 * qJD(1);
t2 = t131 * t117 - t19 * t99 - t59 * t54;
t136 = t2 * qJD(1);
t134 = t141 * t117;
t132 = t106 * t97 / 0.2e1 - t108 * t99 / 0.2e1;
t130 = t93 * qJD(5) + t158;
t129 = t99 * t134;
t121 = t197 + t132;
t8 = t121 * t117;
t128 = t8 * qJD(1) - t108 * t167;
t6 = t121 * t115;
t127 = -t6 * qJD(1) - t108 * t166;
t52 = (t111 / 0.2e1 - t112 / 0.2e1) * t99;
t126 = -t52 * qJD(1) + t147;
t125 = t96 * t115 * t169 + t52 * qJD(3);
t60 = t103 * t96;
t124 = t60 * qJD(1) + t139;
t90 = t93 * qJD(3);
t88 = t99 * t167;
t51 = t58 * qJD(5);
t45 = t54 * qJD(3);
t44 = t52 * qJD(5);
t32 = t195 / 0.2e1 + t122 * pkin(3);
t31 = -t164 - t46;
t9 = t190 / 0.2e1 - t41 * t196 + (t197 - t132) * t117;
t7 = t132 * t115 + t59 * t196 + t189;
t3 = [0, 0, 0, 0, t116 * t160, t104 * qJD(3), 0, 0, 0, t109 * t161, t109 * t160, t66 * qJD(4), t10 * qJD(3) + t14 * qJD(4), -t112 * t157 - t96 * t148, -t60 * qJD(5) + t97 * t139, t29 * qJD(3) - t97 * t154, -t27 * qJD(3) - t153 * t97, t157, t1 * qJD(3) + t28 * qJD(4) + t13 * qJD(5), t2 * qJD(3) + t61 * qJD(4) + t12 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, t146, t162, t160, -t161, 0, -t107 * t160 + t150, t107 * t161 + t149, (-t113 * t99 + t188 * t97) * t194, (-t113 * t41 - t120 * t188) * t194 + t32 * qJD(4) + t137, -t44 + (-t147 - t155) * t97, t123 * t97 - 0.2e1 * t99 * t148, t88 + t178, t89 - t180, t130, t193 + (t115 * t138 - t198) * qJD(3) + t9 * qJD(5), (t117 * t138 + t199) * qJD(3) + t7 * qJD(5) + t136; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t172, t32 * qJD(3) + t183, 0, 0, 0, 0, 0, t179, t173; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t125, -t124, t141 * t56, t129, t90, t9 * qJD(3) - t19 * qJD(5) + t184, t7 * qJD(3) + t18 * qJD(5) + t185; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t161, -t160, 0, (-t145 - t192) * t194, 0, 0, 0, 0, 0, t200, t51 + t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45 - t153, t154 + t174; 0, 0, 0, 0, -t146, -t162, 0, 0, 0, -t150, -t149, 0, t33 * qJD(4) - t137, t97 * t155 - t44, 0.2e1 * t115 * t129, t51 - t178, -t47 + t180, -t130, -t8 * qJD(5) - t165 * t99 - t193, t56 * qJD(4) + t6 * qJD(5) - t136; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t148, t103 * qJD(5), 0, 0, 0, t108 * t164, t108 * t163; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t177, 0, 0, 0, 0, 0, -t152, t176; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t126, -t123, t163 + t175, t31, -t171, -t106 * t163 - t128, t106 * t164 - t127; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t172, -t33 * qJD(3) - t183, 0, 0, 0, 0, 0, -t179 - t200, -t56 * qJD(3) - t163 * t97 - t173; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t177, 0, 0, 0, 0, 0, t152, -t176; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, t134; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t125, t124, t115 * t158 - t174, t152 * t97 + t45, t90, t8 * qJD(3) + t54 * qJD(4) - t184, -t6 * qJD(3) + t165 * t97 - t185; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t126, t123, -t175, t46, t171, t128, t127; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, t97 * t169; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t3;
