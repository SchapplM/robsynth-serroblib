% Calculate minimal parameter regressor of coriolis matrix for
% S5PRRRR6
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
% cmat_reg [(5*%NQJ)%x21]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5PRRRR6_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR6_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR6_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR6_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:10:15
% EndTime: 2019-12-05 17:10:19
% DurationCPUTime: 1.14s
% Computational Cost: add. (1302->159), mult. (2974->202), div. (0->0), fcn. (3108->8), ass. (0->137)
t165 = qJD(4) + qJD(5);
t131 = cos(qJ(3));
t195 = t131 * pkin(2);
t121 = -pkin(3) - t195;
t210 = pkin(3) / 0.2e1 - t121 / 0.2e1;
t166 = qJD(2) + qJD(3);
t127 = sin(qJ(4));
t198 = cos(qJ(5));
t159 = t198 * t127;
t126 = sin(qJ(5));
t130 = cos(qJ(4));
t185 = t126 * t130;
t137 = t159 + t185;
t158 = t198 * t130;
t186 = t126 * t127;
t96 = -t158 + t186;
t36 = -t137 ^ 2 + t96 ^ 2;
t209 = t166 * t36;
t114 = -t127 ^ 2 + t130 ^ 2;
t208 = t166 * t114;
t162 = -t195 / 0.2e1;
t207 = t162 - t210;
t200 = pkin(7) + pkin(8);
t109 = t200 * t127;
t110 = t200 * t130;
t136 = t185 / 0.2e1 + t159 / 0.2e1;
t128 = sin(qJ(3));
t129 = sin(qJ(2));
t132 = cos(qJ(2));
t99 = t128 * t129 - t131 * t132;
t22 = (-t137 / 0.2e1 + t136) * t99;
t184 = t22 * qJD(1);
t206 = -t184 + t165 * (t126 * t109 - t198 * t110);
t135 = t158 / 0.2e1 - t186 / 0.2e1;
t23 = (t96 / 0.2e1 + t135) * t99;
t183 = t23 * qJD(1);
t205 = -t183 + t165 * (t198 * t109 + t126 * t110);
t120 = t128 * pkin(2) + pkin(7);
t194 = pkin(8) + t120;
t90 = t194 * t127;
t91 = t194 * t130;
t204 = -t184 + t165 * (t126 * t90 - t198 * t91);
t203 = -t183 + t165 * (t126 * t91 + t198 * t90);
t201 = -t99 / 0.2e1;
t197 = pkin(4) * t126;
t196 = pkin(4) * t127;
t193 = pkin(2) * qJD(3);
t192 = pkin(3) * qJD(3);
t191 = qJD(2) * pkin(2);
t122 = -t130 * pkin(4) - pkin(3);
t107 = t122 - t195;
t190 = t107 * t96;
t189 = t107 * t137;
t188 = t122 * t96;
t187 = t122 * t137;
t87 = t96 * t196;
t34 = t87 + t189;
t178 = t34 * qJD(2);
t88 = t137 * t196;
t35 = t88 - t190;
t177 = t35 * qJD(2);
t172 = qJD(2) * t107;
t171 = qJD(2) * t121;
t170 = qJD(3) * t122;
t169 = qJD(5) * t107;
t168 = qJD(5) * t122;
t167 = t127 * qJD(4);
t125 = t130 * qJD(4);
t164 = t128 * t193;
t163 = t128 * t191;
t161 = t96 * t172;
t160 = t137 * t172;
t157 = t127 * t171;
t156 = t130 * t171;
t155 = t122 / 0.2e1 + t107 / 0.2e1;
t154 = t198 * qJD(4);
t153 = t198 * qJD(5);
t152 = pkin(2) * t166;
t151 = t166 * t96;
t56 = t165 * t137;
t150 = t96 * t163;
t149 = t137 * t163;
t100 = t128 * t132 + t131 * t129;
t148 = t100 * t165;
t58 = t166 * t100;
t147 = t166 * t127;
t146 = t127 * t163;
t145 = t128 * t152;
t134 = t136 * t195;
t26 = -t137 * t155 - t134;
t12 = -t87 + t26;
t43 = t87 + t187;
t144 = -t12 * qJD(2) + t43 * qJD(3);
t133 = t135 * t195;
t27 = t155 * t96 - t133;
t13 = -t88 + t27;
t44 = t88 - t188;
t143 = -t13 * qJD(2) + t44 * qJD(3);
t142 = t162 + t210;
t65 = t142 * t127;
t141 = t65 * qJD(2) + t127 * t192;
t66 = t142 * t130;
t140 = t66 * qJD(2) + t130 * t192;
t139 = -t26 * qJD(2) + t137 * t170;
t138 = -t27 * qJD(2) - t96 * t170;
t28 = -t188 / 0.2e1 - t190 / 0.2e1 - t133;
t29 = t187 / 0.2e1 + t189 / 0.2e1 - t134;
t115 = t127 * t125;
t113 = t127 * t164;
t108 = t114 * qJD(4);
t89 = t130 * t147;
t82 = t137 * t164;
t81 = t96 * t164;
t68 = t207 * t130;
t67 = t207 * t127;
t57 = t166 * t99;
t55 = t165 * t96;
t54 = t99 * t130;
t53 = t99 * t127;
t31 = t137 * t151;
t30 = t96 * t56;
t25 = t136 * t99 - t137 * t201;
t24 = t135 * t99 + t96 * t201;
t17 = t53 * qJD(4) - t130 * t58;
t16 = t54 * qJD(4) + t100 * t147;
t15 = t88 + t28;
t14 = t87 + t29;
t9 = t165 * t36;
t8 = t165 * t22;
t7 = t165 * t23;
t6 = t166 * t23;
t5 = t166 * t22;
t4 = t137 * t58 + t165 * t24;
t3 = t100 * t151 + t165 * t25;
t2 = t96 * t148 + t166 * t25;
t1 = t137 * t148 + t166 * t24;
t10 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t129 * qJD(2), -t132 * qJD(2), 0, -t58, t57, 0, 0, 0, 0, 0, t17, t16, 0, 0, 0, 0, 0, t3, t4; 0, 0, 0, 0, 0, -t58, t57, 0, 0, 0, 0, 0, t17, t16, 0, 0, 0, 0, 0, t3, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t100 * t125 + t166 * t53, t100 * t167 + t166 * t54, 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, -t7; 0, 0, 0, 0, 0, -t164, -t131 * t193, t115, t108, 0, 0, 0, t121 * t167 - t130 * t164, t121 * t125 + t113, -t30, t9, 0, 0, 0, t34 * qJD(4) + t137 * t169 + t81, t35 * qJD(4) - t96 * t169 + t82; 0, 0, 0, 0, 0, -t145, -t131 * t152, t115, t108, 0, 0, 0, t67 * qJD(4) - t130 * t145, t68 * qJD(4) + t113 + t146, -t30, t9, 0, 0, 0, t14 * qJD(4) + t29 * qJD(5) + t150 + t81, t15 * qJD(4) + t28 * qJD(5) + t149 + t82; 0, 0, 0, 0, 0, 0, 0, t89, t208, t125, -t167, 0, t67 * qJD(3) - t120 * t125 + t157, t68 * qJD(3) + t120 * t167 + t156, -t31, t209, -t55, -t56, 0, t14 * qJD(3) + t178 + t204, t15 * qJD(3) + t177 + t203; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, t209, -t55, -t56, 0, t29 * qJD(3) + t160 + t204, t28 * qJD(3) - t161 + t203; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, -t7; 0, 0, 0, 0, 0, t163, t131 * t191, t115, t108, 0, 0, 0, -t65 * qJD(4) + t130 * t163, -t66 * qJD(4) - t146, -t30, t9, 0, 0, 0, -t12 * qJD(4) - t26 * qJD(5) - t150, -t13 * qJD(4) - t27 * qJD(5) - t149; 0, 0, 0, 0, 0, 0, 0, t115, t108, 0, 0, 0, -pkin(3) * t167, -pkin(3) * t125, -t30, t9, 0, 0, 0, t43 * qJD(4) + t137 * t168, t44 * qJD(4) - t96 * t168; 0, 0, 0, 0, 0, 0, 0, t89, t208, t125, -t167, 0, -pkin(7) * t125 - t141, pkin(7) * t167 - t140, -t31, t209, -t55, -t56, 0, t144 + t206, t143 + t205; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, t209, -t55, -t56, 0, t139 + t206, t138 + t205; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6; 0, 0, 0, 0, 0, 0, 0, -t89, -t208, 0, 0, 0, t65 * qJD(3) - t157, t66 * qJD(3) - t156, t31, -t209, 0, 0, 0, t12 * qJD(3) - t178 + t184, t13 * qJD(3) - t177 + t183; 0, 0, 0, 0, 0, 0, 0, -t89, -t208, 0, 0, 0, t141, t140, t31, -t209, 0, 0, 0, -t144 + t184, -t143 + t183; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(5) * t197, -pkin(4) * t153; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t165 * t197, (-t154 - t153) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, -t209, 0, 0, 0, t26 * qJD(3) - t160 + t184, t27 * qJD(3) + t161 + t183; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, -t209, 0, 0, 0, -t139 + t184, -t138 + t183; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4) * t197, pkin(4) * t154; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t10;
