% Calculate minimal parameter regressor of coriolis matrix for
% S5PPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x19]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5PPRRR4_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR4_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR4_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PPRRR4_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:19:46
% EndTime: 2019-12-05 15:19:51
% DurationCPUTime: 1.44s
% Computational Cost: add. (1238->165), mult. (3872->316), div. (0->0), fcn. (4393->12), ass. (0->140)
t126 = sin(pkin(6));
t206 = sin(pkin(5));
t152 = cos(pkin(11)) * t206;
t207 = cos(pkin(6));
t208 = cos(pkin(5));
t217 = t208 * t126 + t207 * t152;
t127 = sin(qJ(5));
t121 = t127 ^ 2;
t130 = cos(qJ(5));
t123 = t130 ^ 2;
t112 = t123 - t121;
t128 = sin(qJ(4));
t181 = t128 * qJD(3);
t171 = t130 * t181;
t216 = t112 * qJD(4) - 0.2e1 * t127 * t171;
t131 = cos(qJ(4));
t211 = t131 * pkin(9);
t212 = t128 * pkin(4);
t107 = -t211 + t212;
t215 = -t107 / 0.2e1;
t214 = -t128 / 0.2e1;
t213 = -t131 / 0.2e1;
t129 = sin(qJ(3));
t132 = cos(qJ(3));
t158 = sin(pkin(11)) * t206;
t67 = t129 * t158 - t217 * t132;
t28 = t67 * t131;
t68 = t217 * t129 + t132 * t158;
t210 = t68 * t127;
t209 = t68 * t130;
t205 = t126 * t129;
t204 = t126 * t132;
t122 = t128 ^ 2;
t203 = t127 * t122;
t202 = t127 * t128;
t201 = t127 * t129;
t200 = t127 * t131;
t199 = t128 * t130;
t198 = t129 * t130;
t197 = t130 * t107;
t196 = t130 * t122;
t195 = t130 * t131;
t194 = t131 * t132;
t193 = t132 * t122;
t117 = pkin(8) * t202;
t153 = -t131 * pkin(4) - t128 * pkin(9);
t105 = -pkin(3) + t153;
t177 = pkin(8) * t200;
t80 = -t130 * t105 + t177;
t41 = t80 * t128 + (-t117 + t197) * t131;
t192 = t41 * qJD(3);
t176 = pkin(8) * t195;
t81 = t127 * t105 + t176;
t42 = t107 * t200 + (-t81 + t176) * t128;
t191 = t42 * qJD(3);
t124 = t131 ^ 2;
t113 = t124 - t122;
t99 = t113 * t127;
t190 = t99 * qJD(3);
t189 = qJD(3) * t126;
t188 = qJD(4) * t127;
t187 = qJD(4) * t130;
t186 = qJD(5) * t127;
t185 = qJD(5) * t130;
t184 = qJD(5) * t131;
t100 = t130 * t124 - t196;
t183 = t100 * qJD(3);
t182 = t113 * qJD(3);
t180 = t128 * qJD(4);
t179 = t131 * qJD(3);
t178 = t131 * qJD(4);
t175 = pkin(3) * t181;
t174 = pkin(3) * t179;
t173 = t130 * t204;
t172 = t28 / 0.2e1;
t170 = t127 * t184;
t169 = t130 * t184;
t168 = t132 * t189;
t167 = t127 * t185;
t166 = t127 * t187;
t165 = t128 * t178;
t164 = t128 * t179;
t163 = t130 * t180;
t162 = -t202 / 0.2e1;
t161 = t202 / 0.2e1;
t160 = t199 / 0.2e1;
t159 = -t194 / 0.2e1;
t156 = -qJD(5) + t179;
t154 = t127 * t163;
t86 = -t126 * t152 + t208 * t207;
t44 = t86 * t128 + t68 * t131;
t151 = t67 * t127 + t44 * t130;
t150 = t44 * t127 - t67 * t130;
t149 = t156 * t128;
t148 = t211 / 0.2e1 - t212 / 0.2e1;
t140 = t215 + t148;
t70 = t140 * t130;
t147 = pkin(4) * t188 - t70 * qJD(3);
t69 = t140 * t127;
t146 = pkin(4) * t187 + t69 * qJD(3);
t89 = t128 * t207 + t131 * t205;
t145 = t127 * t89 + t173;
t144 = t127 * t204 - t130 * t89;
t90 = (t121 / 0.2e1 - t123 / 0.2e1) * t128;
t143 = -t90 * qJD(3) + t166;
t142 = t130 * t149;
t139 = t127 * qJD(3) * t196 + t90 * qJD(4);
t98 = t112 * t122;
t138 = t98 * qJD(3) + 0.2e1 * t154;
t43 = t68 * t128 - t86 * t131;
t137 = t68 / 0.2e1 + t44 * t213 + t43 * t214;
t88 = t128 * t205 - t131 * t207;
t133 = t205 / 0.2e1 + t89 * t213 + t88 * t214;
t21 = t133 * t130;
t5 = t137 * t130;
t66 = -pkin(8) * t196 - t81 * t131;
t135 = t5 * qJD(1) + t21 * qJD(2) + t66 * qJD(3);
t22 = t133 * t127;
t6 = t137 * t127;
t65 = -pkin(8) * t203 - t80 * t131;
t134 = -t6 * qJD(1) - t22 * qJD(2) - t65 * qJD(3);
t119 = t180 / 0.2e1;
t93 = (t179 - qJD(5) / 0.2e1) * t128;
t87 = t90 * qJD(5);
t60 = t117 + t197 / 0.2e1 + t148 * t130;
t59 = pkin(8) * t199 + (-t148 + t215) * t127;
t48 = t88 * t130;
t47 = t88 * t127;
t26 = t67 * t128;
t24 = t144 * t213 + t88 * t160 + (t127 * t159 + t198 / 0.2e1) * t126;
t23 = t145 * t213 + t88 * t162 + (t130 * t159 - t201 / 0.2e1) * t126;
t16 = t144 * t128 / 0.2e1 + t89 * t160 + t161 * t204;
t15 = t89 * t161 + (t145 + t173) * t214;
t12 = t43 * t130;
t10 = t43 * t127;
t8 = t151 * t131 / 0.2e1 + t43 * t160 + t127 * t172 + t209 / 0.2e1;
t7 = t150 * t213 + t43 * t162 + t130 * t172 - t210 / 0.2e1;
t4 = t151 * t214 + t44 * t160 + t67 * t162;
t3 = t150 * t214 + t67 * t160 + t44 * t161;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, -t68 * qJD(3), t67 * qJD(3), 0, 0, 0, 0, 0, t26 * qJD(4) - t68 * t179, t28 * qJD(4) + t68 * t181, 0, 0, 0, 0, 0, (-(t67 * t200 + t209) * t131 - t67 * t203) * qJD(3) + t3 * qJD(4) + t8 * qJD(5), ((-t195 * t67 + t210) * t131 - t67 * t196) * qJD(3) + t4 * qJD(4) + t7 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26 * qJD(3) - t44 * qJD(4), t28 * qJD(3) + t43 * qJD(4), 0, 0, 0, 0, 0, t3 * qJD(3) + t10 * qJD(5) - t187 * t44, t4 * qJD(3) + t12 * qJD(5) + t188 * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8 * qJD(3) + t10 * qJD(4) - qJD(5) * t151, t7 * qJD(3) + t12 * qJD(4) + qJD(5) * t150; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, -t129 * t189, -t168, 0, 0, 0, 0, 0, (-t129 * t179 - t132 * t180) * t126, (t129 * t181 - t132 * t178) * t126, 0, 0, 0, 0, 0, (-(-t127 * t194 + t198) * t131 + t127 * t193) * t189 + t15 * qJD(4) + t24 * qJD(5), ((t130 * t194 + t201) * t131 + t130 * t193) * t189 + t16 * qJD(4) + t23 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t89 * qJD(4) - t128 * t168, t88 * qJD(4) - t131 * t168, 0, 0, 0, 0, 0, t15 * qJD(3) + t47 * qJD(5) - t187 * t89, t16 * qJD(3) + t48 * qJD(5) + t188 * t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24 * qJD(3) + t47 * qJD(4) + qJD(5) * t144, t23 * qJD(3) + t48 * qJD(4) + qJD(5) * t145; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5 * qJD(5), t6 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21 * qJD(5), t22 * qJD(5); 0, 0, 0, 0, 0, t165, t113 * qJD(4), 0, 0, 0, -pkin(3) * t180, -pkin(3) * t178, -t122 * t167 + t123 * t165, -t98 * qJD(5) - 0.2e1 * t131 * t154, -t100 * qJD(4) + t128 * t170, t99 * qJD(4) + t128 * t169, -t165, -t41 * qJD(4) - t66 * qJD(5), t42 * qJD(4) + t65 * qJD(5); 0, 0, 0, 0, 0, t164, t182, t178, -t180, 0, -pkin(8) * t178 - t175, pkin(8) * t180 - t174, -t87 + (t123 * t181 + t166) * t131, -0.2e1 * t128 * t167 + t216 * t131, t127 * t180 - t183, t163 + t190, -t93, -t192 + (t127 * t153 - t176) * qJD(4) + t60 * qJD(5), t191 + (t130 * t153 + t177) * qJD(4) + t59 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t139, -t138, t127 * t149, t142, t119, t60 * qJD(4) - t81 * qJD(5) - t135, t59 * qJD(4) + t80 * qJD(5) - t134; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -t164, -t182, 0, 0, 0, t175, t174, -t123 * t164 - t87, 0.2e1 * t127 * t142, -t169 + t183, t170 - t190, t93, t70 * qJD(5) + t192, -t69 * qJD(5) - t191; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t167, t112 * qJD(5), 0, 0, 0, -pkin(4) * t186, -pkin(4) * t185; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t143, t216, -t156 * t130, t156 * t127, -t181 / 0.2e1, -pkin(9) * t185 - t147, pkin(9) * t186 - t146; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5 * qJD(3), -t6 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21 * qJD(3), -t22 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t139, t138, (-t127 * t181 + t187) * t131, (-t171 - t188) * t131, t119, -t70 * qJD(4) + t135, t69 * qJD(4) + t134; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t143, -t216, t130 * t179, -t127 * t179, t181 / 0.2e1, t147, t146; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t1;
