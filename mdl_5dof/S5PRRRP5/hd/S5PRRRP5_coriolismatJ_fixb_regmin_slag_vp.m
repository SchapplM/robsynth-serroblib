% Calculate minimal parameter regressor of coriolis matrix for
% S5PRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x20]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:49
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5PRRRP5_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP5_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP5_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP5_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:49:20
% EndTime: 2019-12-05 16:49:25
% DurationCPUTime: 1.31s
% Computational Cost: add. (1679->148), mult. (3829->241), div. (0->0), fcn. (4014->6), ass. (0->133)
t130 = cos(qJ(4));
t199 = t130 * pkin(3);
t123 = pkin(4) + t199;
t205 = -t123 / 0.2e1;
t148 = t199 / 0.2e1 + t205;
t161 = qJD(3) + qJD(4);
t131 = cos(qJ(3));
t184 = t130 * t131;
t127 = sin(qJ(4));
t128 = sin(qJ(3));
t187 = t127 * t128;
t218 = t184 - t187;
t210 = pkin(6) + pkin(7);
t117 = t210 * t131;
t111 = t130 * t117;
t116 = t210 * t128;
t188 = t127 * t116;
t221 = -t111 + t188;
t46 = -qJ(5) * t218 + t221;
t225 = -t46 / 0.2e1;
t185 = t130 * t128;
t186 = t127 * t131;
t108 = t185 + t186;
t132 = cos(qJ(2));
t92 = t132 * t108;
t223 = -t92 / 0.2e1;
t52 = t130 * t116 + t127 * t117;
t222 = -t108 * qJ(5) - t52;
t124 = -t131 * pkin(3) - pkin(2);
t170 = qJD(2) * t124;
t138 = -t186 / 0.2e1 - t185 / 0.2e1;
t206 = t108 / 0.2e1;
t62 = (t206 + t138) * t132;
t175 = t62 * qJD(1);
t220 = -t108 * t170 + t175;
t200 = t128 * pkin(3);
t60 = t124 * t108 - t200 * t218;
t219 = -t60 * qJD(2) + t175;
t137 = -t184 / 0.2e1 + t187 / 0.2e1;
t208 = t218 / 0.2e1;
t63 = (t208 + t137) * t132;
t174 = t63 * qJD(1);
t216 = t161 * t52 - t174;
t129 = sin(qJ(2));
t91 = t108 * t129;
t198 = t91 * t46;
t93 = t218 * t129;
t212 = -t93 / 0.2e1;
t215 = t198 / 0.2e1 + (t93 / 0.2e1 + t212) * t222;
t214 = t108 ^ 2;
t209 = t46 * pkin(4);
t207 = -t218 / 0.2e1;
t152 = -t111 / 0.2e1;
t204 = t129 / 0.2e1;
t203 = pkin(3) * t127;
t202 = pkin(4) * t218;
t201 = pkin(4) * t108;
t194 = pkin(3) * qJD(4);
t193 = pkin(4) * qJD(4);
t192 = qJD(3) * pkin(3);
t94 = t132 * t218;
t19 = -t132 * t129 + t91 * t92 + t93 * t94;
t183 = t19 * qJD(1);
t153 = t127 * t208;
t36 = (t205 - pkin(4) / 0.2e1) * t108 + (t153 - t128 / 0.2e1) * pkin(3);
t182 = t36 * qJD(2);
t141 = pkin(4) / 0.2e1 + t148;
t38 = t141 * t218;
t181 = t38 * qJD(2);
t105 = t218 ^ 2;
t51 = t105 - t214;
t180 = t51 * qJD(2);
t61 = t108 * t200 + t124 * t218;
t176 = t61 * qJD(2);
t74 = t105 + t214;
t173 = t74 * qJD(2);
t77 = t152 + t111 / 0.2e1;
t172 = t77 * qJD(2);
t171 = qJD(2) * t108;
t169 = qJD(2) * t131;
t168 = t218 * qJD(4);
t167 = t108 * qJD(4);
t121 = -t128 ^ 2 + t131 ^ 2;
t166 = t121 * qJD(2);
t165 = t128 * qJD(3);
t164 = t129 * qJD(2);
t163 = t131 * qJD(3);
t162 = t132 * qJD(2);
t160 = pkin(2) * t128 * qJD(2);
t159 = pkin(2) * t169;
t158 = pkin(4) * t171;
t157 = t127 * t194;
t156 = t218 * t170;
t154 = t128 * t169;
t150 = pkin(3) * t161;
t73 = t161 * t108;
t147 = t200 + t201;
t136 = t123 * t223 + t94 * t203 / 0.2e1;
t139 = t132 * t147;
t1 = t139 / 0.2e1 - (t46 / 0.2e1 + t225) * t91 + t136;
t87 = t124 - t202;
t9 = t87 * t147;
t145 = -t1 * qJD(1) + t9 * qJD(2);
t10 = t87 * t201;
t135 = -t198 / 0.2e1 + t215;
t4 = (t223 + t92 / 0.2e1) * pkin(4) + t135;
t144 = t4 * qJD(1) + t10 * qJD(2);
t18 = -t108 * t222 - t218 * t46;
t140 = -t91 * t206 + t93 * t207;
t27 = t204 + t140;
t143 = -t27 * qJD(1) + t18 * qJD(2);
t26 = t141 * t93;
t133 = t148 * t46;
t6 = -t209 / 0.2e1 - t133;
t95 = (-t123 + t199) * t203;
t134 = -t26 * qJD(1) - t6 * qJD(2) - t95 * qJD(3);
t75 = t218 * t171;
t72 = t161 * t218;
t65 = t138 * t132 + t223;
t64 = (t137 + t207) * t132;
t56 = t63 * qJD(2);
t54 = t62 * qJD(2);
t53 = 0.2e1 * t152 + t188;
t37 = -t202 / 0.2e1 + t148 * t218;
t35 = pkin(3) * t153 + t108 * t205 + t200 / 0.2e1 + t201 / 0.2e1;
t28 = t204 - t140;
t25 = pkin(4) * t212 + t148 * t93;
t21 = t65 * qJD(2) - t161 * t93;
t20 = t64 * qJD(2) + t161 * t91;
t5 = t209 / 0.2e1 - t133;
t3 = 0.2e1 * t223 * pkin(4) + t135;
t2 = t91 * t225 - t139 / 0.2e1 + t136 + t215;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19 * qJD(2); 0, 0, -t164, -t162, 0, 0, 0, 0, 0, -t131 * t164 - t132 * t165, t128 * t164 - t132 * t163, 0, 0, 0, 0, 0, t161 * t65 - t164 * t218, t108 * t164 + t161 * t64, (t92 * t108 + t218 * t94) * qJD(2), t183 + (t129 * t87 - t222 * t92 - t46 * t94) * qJD(2) + t2 * qJD(3) + t3 * qJD(4) + t28 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, -t128 * t162 - t129 * t163, t129 * t165 - t131 * t162, 0, 0, 0, 0, 0, t21, t20, 0, t2 * qJD(2) + (-t93 * t123 - t91 * t203) * qJD(3) + t25 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, t20, 0, t3 * qJD(2) + t25 * qJD(3) - t93 * t193; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28 * qJD(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t161 * t62, -t161 * t63, 0, -t1 * qJD(3) + t4 * qJD(4) - t27 * qJD(5) - t183; 0, 0, 0, 0, t128 * t163, t121 * qJD(3), 0, 0, 0, -pkin(2) * t165, -pkin(2) * t163, t218 * t73, t161 * t51, 0, 0, 0, t60 * qJD(3) + t124 * t167, t61 * qJD(3) + t124 * t168, t74 * qJD(5), t9 * qJD(3) + t10 * qJD(4) + t18 * qJD(5); 0, 0, 0, 0, t154, t166, t163, -t165, 0, -pkin(6) * t163 - t160, pkin(6) * t165 - t159, t75, t180, t72, -t73, 0, qJD(3) * t221 + t53 * qJD(4) - t219, t176 + t216, (-t108 * t203 - t123 * t218) * qJD(3) + t37 * qJD(4), (t46 * t123 + t203 * t222) * qJD(3) + t5 * qJD(4) + t35 * qJD(5) + t145; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, t180, t72, -t73, 0, t53 * qJD(3) + qJD(4) * t221 - t220, t156 + t216, -pkin(4) * t168 + t37 * qJD(3), t5 * qJD(3) + t193 * t46 + t144; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t173, t35 * qJD(3) + t143; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, t56, 0, t1 * qJD(2) + t26 * qJD(4); 0, 0, 0, 0, -t154, -t166, 0, 0, 0, t160, t159, -t75, -t180, 0, 0, 0, t77 * qJD(4) + t219, t174 - t176, t38 * qJD(4), t6 * qJD(4) + t36 * qJD(5) - t145; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t157, -t130 * t194, 0, t95 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t127 * t150 + t172, -t130 * t150, t181, -pkin(4) * t157 - t134; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t182; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, t56, 0, -t4 * qJD(2) - t26 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t75, -t180, 0, 0, 0, -t77 * qJD(3) + t220, t174 - t156, -t38 * qJD(3), -t6 * qJD(3) - qJD(5) * t201 - t144; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t127 * t192 - t172, t130 * t192, -t181, t134; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t158; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27 * qJD(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t173, pkin(4) * t167 - t36 * qJD(3) - t143; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t182; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t158; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t7;
