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
% cmat_reg [(5*%NQJ)%x22]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 12:07
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
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
% StartTime: 2021-01-15 12:06:31
% EndTime: 2021-01-15 12:06:36
% DurationCPUTime: 1.54s
% Computational Cost: add. (1689->140), mult. (3428->237), div. (0->0), fcn. (3624->8), ass. (0->140)
t120 = sin(qJ(5));
t211 = 0.2e1 * t120;
t118 = sin(pkin(9));
t123 = cos(qJ(3));
t111 = sin(pkin(8)) * pkin(1) + pkin(6);
t209 = qJ(4) + t111;
t171 = t209 * t123;
t192 = cos(pkin(9));
t121 = sin(qJ(3));
t95 = t209 * t121;
t43 = t118 * t171 + t192 * t95;
t210 = t43 / 0.2e1;
t207 = -t118 * t95 + t192 * t171;
t198 = t120 * t207;
t122 = cos(qJ(5));
t196 = t122 * t207;
t99 = t118 * t121 - t192 * t123;
t56 = t120 * t99;
t49 = t56 * qJD(5);
t145 = t192 * t121;
t180 = t118 * t123;
t101 = t145 + t180;
t163 = t101 * qJD(3);
t91 = t122 * t163;
t208 = -t91 + t49;
t97 = t99 ^ 2;
t98 = t101 ^ 2;
t68 = t97 + t98;
t116 = t120 ^ 2;
t117 = t122 ^ 2;
t107 = t117 - t116;
t58 = t120 * t101;
t144 = 0.2e1 * t122 * t58;
t127 = qJD(1) * t144 - qJD(3) * t107;
t204 = t121 * pkin(3);
t202 = qJD(3) * pkin(3);
t113 = -cos(pkin(8)) * pkin(1) - pkin(2);
t103 = -pkin(3) * t123 + t113;
t124 = t99 * pkin(4) - t101 * pkin(7) + t103;
t18 = -t122 * t124 + t198;
t61 = pkin(4) * t101 + pkin(7) * t99 + t204;
t195 = t122 * t61;
t1 = t195 * t99 + (-t18 + t198) * t101;
t201 = t1 * qJD(1);
t199 = t118 * t99;
t197 = t120 * t61;
t158 = -t43 / 0.2e1 + t210;
t4 = t158 * t101;
t194 = t4 * qJD(1);
t193 = t43 * t122;
t12 = t18 * t99 - t43 * t58;
t191 = qJD(1) * t12;
t19 = t120 * t124 + t196;
t13 = t101 * t193 - t19 * t99;
t190 = qJD(1) * t13;
t14 = t43 * t101 - t207 * t99;
t189 = qJD(1) * t14;
t159 = t98 - t97;
t27 = t159 * t120;
t187 = qJD(1) * t27;
t28 = t68 * t120;
t186 = qJD(1) * t28;
t29 = t159 * t122;
t185 = qJD(1) * t29;
t34 = t101 * t103 + t99 * t204;
t184 = qJD(1) * t34;
t35 = t101 * t204 - t103 * t99;
t183 = qJD(1) * t35;
t63 = t68 * t122;
t182 = qJD(1) * t63;
t60 = t122 * t99;
t181 = qJD(3) * t60;
t146 = t192 * t101;
t126 = -t199 / 0.2e1 - t146 / 0.2e1;
t33 = (-t121 / 0.2e1 + t126) * pkin(3);
t178 = t33 * qJD(1);
t48 = t56 * qJD(1);
t177 = t58 * qJD(1);
t176 = t60 * qJD(1);
t175 = t68 * qJD(1);
t96 = t145 / 0.2e1 + t180 / 0.2e1;
t174 = t96 * qJD(1);
t173 = t99 * qJD(1);
t172 = t99 * qJD(3);
t170 = qJD(1) * t122;
t169 = qJD(1) * t123;
t168 = qJD(3) * t122;
t167 = qJD(4) * t122;
t166 = qJD(5) * t120;
t165 = qJD(5) * t122;
t164 = t101 * qJD(1);
t108 = -t121 ^ 2 + t123 ^ 2;
t162 = t108 * qJD(1);
t161 = t121 * qJD(3);
t160 = t123 * qJD(3);
t157 = t99 * t164;
t156 = t99 * t163;
t155 = t117 * t164;
t154 = t120 * t168;
t153 = t101 * t166;
t152 = t101 * t165;
t151 = t113 * t121 * qJD(1);
t150 = t113 * t169;
t149 = t120 * t165;
t148 = t121 * t169;
t147 = t122 * t164;
t143 = -qJD(5) - t173;
t142 = qJD(3) * t144;
t10 = t103 * t204;
t140 = t10 * qJD(1) + t4 * qJD(2);
t2 = -t197 * t99 + (-t19 + t196) * t101;
t139 = t2 * qJD(1);
t110 = pkin(3) * t118 + pkin(7);
t112 = -t192 * pkin(3) - pkin(4);
t137 = -t101 * t110 - t112 * t99;
t136 = t143 * t122;
t135 = t110 * t99 / 0.2e1 - t112 * t101 / 0.2e1;
t125 = t61 / 0.2e1 + t135;
t8 = t158 * t120 - t125 * t122;
t134 = -qJD(3) * t112 * t120 - qJD(1) * t8;
t6 = t125 * t120 + t158 * t122;
t133 = -qJD(1) * t6 - t112 * t168;
t132 = qJD(5) * t96 + t157;
t131 = t101 * t136;
t54 = (t116 / 0.2e1 - t117 / 0.2e1) * t101;
t130 = -qJD(1) * t54 + t154;
t129 = t120 * t98 * t170 + qJD(3) * t54;
t62 = t107 * t98;
t128 = qJD(1) * t62 + t142;
t92 = t96 * qJD(3);
t90 = t120 * t163;
t53 = t60 * qJD(5);
t47 = t56 * qJD(3);
t46 = t54 * qJD(5);
t32 = t204 / 0.2e1 + t126 * pkin(3);
t31 = -t166 - t48;
t9 = t195 / 0.2e1 - t135 * t122 + t210 * t211;
t7 = t193 / 0.2e1 + t122 * t210 - t197 / 0.2e1 + t135 * t120;
t3 = qJD(3) * t4;
t5 = [0, 0, 0, 0, t121 * t160, t108 * qJD(3), 0, 0, 0, t113 * t161, t113 * t160, t34 * qJD(3), t35 * qJD(3), qJD(4) * t68, qJD(3) * t10 + qJD(4) * t14, -t117 * t156 - t98 * t149, -qJD(5) * t62 + t99 * t142, qJD(3) * t29 - t99 * t153, -qJD(3) * t27 - t99 * t152, t156, qJD(3) * t1 + qJD(4) * t28 + qJD(5) * t13, qJD(3) * t2 + qJD(4) * t63 + qJD(5) * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, t148, t162, t160, -t161, 0, -t111 * t160 + t151, t111 * t161 + t150, -qJD(3) * t207 + t184, qJD(3) * t43 + t183, (-t101 * t118 + t192 * t99) * t202, (-t118 * t43 - t192 * t207) * t202 + t32 * qJD(4) + t140, -t46 + (-t154 - t155) * t99, -0.2e1 * t101 * t149 + t127 * t99, t90 + t185, t91 - t187, t132, t201 + (t137 * t120 - t196) * qJD(3) + t9 * qJD(5), (t122 * t137 + t198) * qJD(3) + t7 * qJD(5) + t139; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t175, qJD(3) * t32 + t189, 0, 0, 0, 0, 0, t186, t182; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t129, -t128, t143 * t58, t131, t92, qJD(3) * t9 - qJD(5) * t19 + t190, qJD(3) * t7 + qJD(5) * t18 + t191; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t161, -t160, -t163, t172, 0, t194 + (-t146 - t199) * t202, 0, 0, 0, 0, 0, t208, t53 + t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47 - t152, t153 + t181; 0, 0, 0, 0, -t148, -t162, 0, 0, 0, -t151, -t150, -qJD(4) * t101 - t184, qJD(4) * t99 - t183, 0, qJD(4) * t33 - t140, t99 * t155 - t46, t131 * t211, t53 - t185, -t49 + t187, -t132, qJD(5) * t8 - t101 * t167 - t201, qJD(4) * t58 + qJD(5) * t6 - t139; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t194, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t149, t107 * qJD(5), 0, 0, 0, t112 * t166, t112 * t165; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t164, t173, 0, t178, 0, 0, 0, 0, 0, -t147, t177; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t130, -t127, t165 + t176, t31, -t174, -t110 * t165 - t134, t110 * t166 - t133; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t163, -t172, -t175, -qJD(3) * t33 - t189, 0, 0, 0, 0, 0, -t186 - t208, -qJD(3) * t58 - t165 * t99 - t182; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t164, -t173, 0, -t178, 0, 0, 0, 0, 0, t147, -t177; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, t136; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t129, t128, t120 * t157 - t181, t99 * t147 + t47, t92, -qJD(3) * t8 + qJD(4) * t56 - t190, -qJD(3) * t6 + t167 * t99 - t191; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t130, t127, -t176, t48, t174, t134, t133; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, t99 * t170; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t5;
