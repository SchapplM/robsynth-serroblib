% Calculate minimal parameter regressor of coriolis matrix for
% S5PRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x22]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 15:53
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5PRRPR4_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR4_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR4_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR4_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 15:52:46
% EndTime: 2021-01-15 15:52:51
% DurationCPUTime: 1.88s
% Computational Cost: add. (1576->144), mult. (3726->242), div. (0->0), fcn. (4172->8), ass. (0->131)
t245 = qJD(3) + qJD(5);
t135 = sin(qJ(5));
t138 = cos(qJ(5));
t134 = cos(pkin(9));
t136 = sin(qJ(3));
t195 = t134 * t136;
t133 = sin(pkin(9));
t139 = cos(qJ(3));
t197 = t133 * t139;
t113 = t195 + t197;
t218 = qJ(4) + pkin(6);
t123 = t218 * t136;
t124 = t218 * t139;
t160 = -t134 * t123 - t124 * t133;
t150 = -pkin(7) * t113 + t160;
t194 = t134 * t139;
t198 = t133 * t136;
t229 = t194 - t198;
t81 = -t133 * t123 + t134 * t124;
t52 = pkin(7) * t229 + t81;
t260 = t245 * (-t135 * t150 - t138 * t52);
t259 = t245 * (t135 * t52 - t138 * t150);
t140 = cos(qJ(2));
t101 = t140 * t229;
t222 = -t135 / 0.2e1;
t99 = t140 * t113;
t144 = -t99 * t222 - t138 * t101 / 0.2e1;
t67 = t135 * t113 - t138 * t229;
t240 = t140 * t67;
t242 = -t240 / 0.2e1 + t144;
t254 = qJD(1) * t242;
t227 = -t99 / 0.2e1;
t145 = t101 * t222 + t138 * t227;
t110 = t138 * t113;
t193 = t135 * t229;
t232 = t110 + t193;
t238 = t140 * t232;
t243 = t238 / 0.2e1 + t145;
t253 = qJD(1) * t243;
t252 = qJD(2) * t242;
t251 = qJD(2) * t243;
t137 = sin(qJ(2));
t100 = t229 * t137;
t241 = t240 / 0.2e1 + t144;
t98 = t113 * t137;
t250 = qJD(2) * t241 + t245 * (t135 * t100 + t138 * t98);
t244 = -t238 / 0.2e1 + t145;
t249 = qJD(2) * t244 + t245 * (-t138 * t100 + t135 * t98);
t15 = t245 * t67;
t236 = -t232 ^ 2 + t67 ^ 2;
t246 = t236 * qJD(2);
t239 = qJD(4) * t67;
t237 = t67 * qJD(2);
t180 = t232 * qJD(2);
t162 = t110 / 0.2e1;
t226 = t229 / 0.2e1;
t225 = -t229 / 0.2e1;
t224 = -t113 / 0.2e1;
t223 = t113 / 0.2e1;
t221 = t137 / 0.2e1;
t220 = pkin(3) * t133;
t219 = t136 * pkin(3);
t217 = qJD(3) * pkin(3);
t130 = -pkin(3) * t139 - pkin(2);
t92 = -pkin(4) * t229 + t130;
t202 = qJD(2) * t92;
t199 = qJD(5) * t92;
t191 = t140 * t136;
t25 = t100 * t101 - t137 * t140 + t98 * t99;
t190 = t25 * qJD(1);
t32 = 0.2e1 * t162 + t193;
t184 = t32 * qJD(2);
t143 = t133 * t226 + t134 * t224;
t48 = (-t136 / 0.2e1 + t143) * pkin(3);
t183 = t48 * qJD(2);
t65 = t162 - t110 / 0.2e1;
t182 = t65 * qJD(2);
t181 = t65 * qJD(5);
t77 = t113 ^ 2 + t229 ^ 2;
t178 = t77 * qJD(2);
t177 = qJD(2) * t139;
t176 = t229 * qJD(2);
t175 = t113 * qJD(2);
t127 = -t136 ^ 2 + t139 ^ 2;
t174 = t127 * qJD(2);
t173 = t136 * qJD(3);
t172 = t137 * qJD(2);
t171 = t139 * qJD(3);
t170 = t140 * qJD(2);
t169 = pkin(2) * t136 * qJD(2);
t168 = pkin(2) * t177;
t167 = t67 * t180;
t166 = t232 * t237;
t163 = t136 * t177;
t146 = t101 * t133 / 0.2e1 + t134 * t227;
t1 = (t191 / 0.2e1 + t146) * pkin(3);
t11 = t130 * t219;
t159 = -t1 * qJD(1) + t11 * qJD(2);
t93 = pkin(4) * t113 + t219;
t12 = t232 * t92 + t67 * t93;
t158 = t12 * qJD(2) - t253;
t13 = t232 * t93 - t67 * t92;
t157 = t13 * qJD(2) - t254;
t24 = -t113 * t160 + t229 * t81;
t147 = t100 * t225 - t98 * t223;
t27 = t221 + t147;
t156 = -qJD(1) * t27 + qJD(2) * t24;
t57 = t113 * t130 - t219 * t229;
t142 = -t197 / 0.2e1 - t195 / 0.2e1;
t60 = (t223 + t142) * t140;
t155 = -qJD(1) * t60 + qJD(2) * t57;
t58 = t113 * t219 + t130 * t229;
t141 = -t194 / 0.2e1 + t198 / 0.2e1;
t61 = (t226 + t141) * t140;
t154 = -qJD(1) * t61 + qJD(2) * t58;
t129 = pkin(3) * t134 + pkin(4);
t107 = -t129 * t138 + t135 * t220;
t153 = qJD(3) * t107;
t108 = t129 * t135 + t138 * t220;
t152 = qJD(3) * t108;
t151 = qJD(3) * t232 + qJD(5) * t32;
t149 = t202 * t232 - t253;
t148 = -t67 * t202 - t254;
t95 = t108 * qJD(5);
t94 = t107 * qJD(5);
t63 = (t142 + t224) * t140;
t62 = (t141 + t225) * t140;
t47 = t219 / 0.2e1 + t143 * pkin(3);
t28 = t221 - t147;
t2 = (-t191 / 0.2e1 + t146) * pkin(3);
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t25, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t172, -t170, 0, 0, 0, 0, 0, -t139 * t172 - t140 * t173, t136 * t172 - t140 * t171, qJD(3) * t63 - t172 * t229, qJD(3) * t62 + t113 * t172, (t101 * t229 + t113 * t99) * qJD(2), t190 + (t101 * t81 + t130 * t137 - t160 * t99) * qJD(2) + t2 * qJD(3) + t28 * qJD(4), 0, 0, 0, 0, 0, t67 * t172 + t245 * t244, t172 * t232 + t245 * t241; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t136 * t170 - t137 * t171, t137 * t173 - t139 * t170, qJD(2) * t63 - qJD(3) * t100, qJD(2) * t62 + qJD(3) * t98, 0, t2 * qJD(2) + (-t100 * t134 - t133 * t98) * t217, 0, 0, 0, 0, 0, t249, t250; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28 * qJD(2), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t249, t250; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t60 * qJD(3), -t61 * qJD(3), 0, -qJD(3) * t1 - qJD(4) * t27 - t190, 0, 0, 0, 0, 0, -t245 * t243, -t245 * t242; 0, 0, 0, 0, t136 * t171, t127 * qJD(3), 0, 0, 0, -pkin(2) * t173, -pkin(2) * t171, t57 * qJD(3), t58 * qJD(3), qJD(4) * t77, qJD(3) * t11 + qJD(4) * t24, -t15 * t232, t245 * t236, 0, 0, 0, qJD(3) * t12 + t199 * t232, qJD(3) * t13 - t199 * t67; 0, 0, 0, 0, t163, t174, t171, -t173, 0, -pkin(6) * t171 - t169, pkin(6) * t173 - t168, -qJD(3) * t81 + t155, -qJD(3) * t160 + t154, (-t113 * t133 - t134 * t229) * t217, t47 * qJD(4) + (t133 * t160 - t134 * t81) * t217 + t159, -t166, t246, -t15, -t151, 0, t158 + t260, t157 + t259; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t178, qJD(3) * t47 + t156, 0, 0, 0, 0, 0, t181, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t167, t246, -t15, -qJD(3) * t32 - qJD(5) * t232, 0, t65 * qJD(4) + t149 + t260, t148 + t259; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60 * qJD(2), t61 * qJD(2), 0, qJD(2) * t1, 0, 0, 0, 0, 0, t251, t252; 0, 0, 0, 0, -t163, -t174, 0, 0, 0, t169, t168, -qJD(4) * t113 - t155, -qJD(4) * t229 - t154, 0, qJD(4) * t48 - t159, t166, -t246, 0, -t181, 0, -qJD(4) * t232 - t158, -t157 + t239; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t95, t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t175, -t176, 0, t183, 0, 0, 0, 0, 0, -t180, t237; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t182, 0, -t152 - t95, t153 + t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27 * qJD(2), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t113 * qJD(3), t229 * qJD(3), -t178, -qJD(3) * t48 - t156, 0, 0, 0, 0, 0, t151, -t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t175, t176, 0, -t183, 0, 0, 0, 0, 0, t180, -t237; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t184, -t237; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t251, t252; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t167, -t246, 0, t65 * qJD(3), 0, -qJD(4) * t32 - t149, -t148 + t239; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t182, 0, t152, -t153; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t184, t237; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t3;
