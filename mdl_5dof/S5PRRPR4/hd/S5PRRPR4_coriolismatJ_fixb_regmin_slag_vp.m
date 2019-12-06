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
% cmat_reg [(5*%NQJ)%x20]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:24
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 16:23:26
% EndTime: 2019-12-05 16:23:31
% DurationCPUTime: 1.49s
% Computational Cost: add. (1486->121), mult. (3484->204), div. (0->0), fcn. (3914->8), ass. (0->115)
t221 = qJD(3) + qJD(5);
t123 = sin(qJ(5));
t126 = cos(qJ(5));
t121 = sin(pkin(9));
t122 = cos(pkin(9));
t124 = sin(qJ(3));
t127 = cos(qJ(3));
t104 = t121 * t127 + t122 * t124;
t196 = -qJ(4) - pkin(6);
t111 = t196 * t124;
t112 = t196 * t127;
t208 = t122 * t111 + t121 * t112;
t212 = -t104 * pkin(7) + t208;
t103 = -t121 * t124 + t122 * t127;
t75 = t121 * t111 - t122 * t112;
t49 = -t103 * pkin(7) - t75;
t239 = t221 * (-t123 * t212 + t126 * t49);
t238 = t221 * (-t123 * t49 - t126 * t212);
t237 = -t104 / 0.2e1;
t201 = -t123 / 0.2e1;
t128 = cos(qJ(2));
t94 = t104 * t128;
t95 = t103 * t128;
t135 = -t94 * t201 - t126 * t95 / 0.2e1;
t61 = -t126 * t103 + t123 * t104;
t216 = t128 * t61;
t218 = -t216 / 0.2e1 + t135;
t232 = qJD(1) * t218;
t203 = -t94 / 0.2e1;
t134 = t126 * t203 + t95 * t201;
t102 = t126 * t104;
t176 = t123 * t103;
t209 = t102 + t176;
t213 = t128 * t209;
t219 = t213 / 0.2e1 + t134;
t231 = qJD(1) * t219;
t230 = t218 * qJD(2);
t229 = t219 * qJD(2);
t217 = t216 / 0.2e1 + t135;
t125 = sin(qJ(2));
t92 = t103 * t125;
t93 = t104 * t125;
t228 = qJD(2) * t217 + t221 * (t123 * t92 + t126 * t93);
t220 = -t213 / 0.2e1 + t134;
t227 = qJD(2) * t220 + t221 * (t123 * t93 - t126 * t92);
t162 = t61 * qJD(5);
t15 = -t61 * qJD(3) - t162;
t211 = -t209 ^ 2 + t61 ^ 2;
t222 = t211 * qJD(2);
t215 = qJD(2) * t61;
t214 = qJD(4) * t61;
t163 = t209 * qJD(2);
t205 = -t92 / 0.2e1;
t144 = t102 / 0.2e1;
t202 = t121 / 0.2e1;
t200 = t125 / 0.2e1;
t199 = pkin(3) * t121;
t197 = t124 * pkin(3);
t195 = qJD(3) * pkin(3);
t118 = -t127 * pkin(3) - pkin(2);
t86 = -t103 * pkin(4) + t118;
t180 = qJD(2) * t86;
t175 = t128 * t124;
t25 = -t128 * t125 + t92 * t95 + t93 * t94;
t173 = t25 * qJD(1);
t32 = 0.2e1 * t144 + t176;
t167 = t32 * qJD(2);
t129 = t103 * t202 + t122 * t237;
t48 = (-t124 / 0.2e1 + t129) * pkin(3);
t166 = t48 * qJD(2);
t59 = t144 - t102 / 0.2e1;
t165 = t59 * qJD(2);
t164 = t59 * qJD(5);
t159 = t209 * qJD(5);
t71 = t103 ^ 2 + t104 ^ 2;
t158 = t71 * qJD(2);
t157 = qJD(2) * t127;
t115 = -t124 ^ 2 + t127 ^ 2;
t156 = t115 * qJD(2);
t155 = t124 * qJD(3);
t154 = t125 * qJD(2);
t153 = t127 * qJD(3);
t152 = t128 * qJD(2);
t151 = pkin(2) * t124 * qJD(2);
t150 = pkin(2) * t157;
t149 = t61 * t163;
t148 = t209 * t215;
t145 = t124 * t157;
t132 = t122 * t203 + t95 * t202;
t1 = (t175 / 0.2e1 + t132) * pkin(3);
t11 = t118 * t197;
t142 = -t1 * qJD(1) + t11 * qJD(2);
t117 = t122 * pkin(3) + pkin(4);
t99 = -t126 * t117 + t123 * t199;
t141 = t99 * qJD(3);
t87 = t104 * pkin(4) + t197;
t12 = t209 * t86 + t87 * t61;
t140 = t12 * qJD(2) - t231;
t13 = t209 * t87 - t61 * t86;
t139 = t13 * qJD(2) - t232;
t24 = t75 * t103 - t104 * t208;
t133 = t103 * t205 + t93 * t237;
t27 = t200 + t133;
t138 = -t27 * qJD(1) + t24 * qJD(2);
t100 = t123 * t117 + t126 * t199;
t137 = t100 * qJD(3);
t136 = qJD(3) * t209 + t32 * qJD(5);
t131 = t180 * t209 - t231;
t130 = -t61 * t180 - t232;
t89 = t100 * qJD(5);
t88 = t99 * qJD(5);
t47 = t197 / 0.2e1 + t129 * pkin(3);
t28 = t200 - t133;
t2 = (-t175 / 0.2e1 + t132) * pkin(3) + (t92 / 0.2e1 + t205) * t208;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25 * qJD(2), 0, 0, 0, 0, 0, 0, 0; 0, 0, -t154, -t152, 0, 0, 0, 0, 0, -t127 * t154 - t128 * t155, t124 * t154 - t128 * t153, (t95 * t103 + t94 * t104) * qJD(2), t173 + (t125 * t118 - t208 * t94 + t95 * t75) * qJD(2) + t2 * qJD(3) + t28 * qJD(4), 0, 0, 0, 0, 0, t61 * t154 + t221 * t220, t154 * t209 + t221 * t217; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t124 * t152 - t125 * t153, t125 * t155 - t127 * t152, 0, t2 * qJD(2) + (-t121 * t93 - t122 * t92) * t195, 0, 0, 0, 0, 0, t227, t228; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28 * qJD(2), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t227, t228; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1 * qJD(3) - t27 * qJD(4) - t173, 0, 0, 0, 0, 0, -t221 * t219, -t221 * t218; 0, 0, 0, 0, t124 * t153, t115 * qJD(3), 0, 0, 0, -pkin(2) * t155, -pkin(2) * t153, t71 * qJD(4), t11 * qJD(3) + t24 * qJD(4), t15 * t209, t221 * t211, 0, 0, 0, t12 * qJD(3) + t86 * t159, t13 * qJD(3) - t162 * t86; 0, 0, 0, 0, t145, t156, t153, -t155, 0, -pkin(6) * t153 - t151, pkin(6) * t155 - t150, (-t103 * t122 - t104 * t121) * t195, t47 * qJD(4) + (t121 * t208 - t122 * t75) * t195 + t142, -t148, t222, t15, -t136, 0, t140 + t239, t139 + t238; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t158, t47 * qJD(3) + t138, 0, 0, 0, 0, 0, t164, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t149, t222, t15, -t32 * qJD(3) - t159, 0, t59 * qJD(4) + t131 + t239, t130 + t238; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1 * qJD(2), 0, 0, 0, 0, 0, t229, t230; 0, 0, 0, 0, -t145, -t156, 0, 0, 0, t151, t150, 0, t48 * qJD(4) - t142, t148, -t222, 0, -t164, 0, -qJD(4) * t209 - t140, -t139 + t214; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t89, t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t166, 0, 0, 0, 0, 0, -t163, t215; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t165, 0, -t137 - t89, t141 + t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27 * qJD(2), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t158, -t48 * qJD(3) - t138, 0, 0, 0, 0, 0, t136, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t166, 0, 0, 0, 0, 0, t163, -t215; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t167, -t215; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t229, t230; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t149, -t222, 0, t59 * qJD(3), 0, -t32 * qJD(4) - t131, -t130 + t214; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t165, 0, t137, -t141; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t167, t215; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t3;
