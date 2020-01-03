% Calculate minimal parameter regressor of coriolis matrix for
% S4RRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% 
% Output:
% cmat_reg [(4*%NQJ)%x24]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:24
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4RRRR3_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR3_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR3_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR3_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:24:38
% EndTime: 2019-12-31 17:24:42
% DurationCPUTime: 1.55s
% Computational Cost: add. (1676->113), mult. (3356->160), div. (0->0), fcn. (3570->6), ass. (0->110)
t233 = qJD(3) + qJD(4);
t240 = qJD(2) + t233;
t115 = sin(qJ(4));
t118 = cos(qJ(4));
t117 = sin(qJ(2));
t201 = pkin(5) + pkin(6);
t102 = t201 * t117;
t119 = cos(qJ(2));
t103 = t201 * t119;
t116 = sin(qJ(3));
t198 = cos(qJ(3));
t48 = t198 * t102 + t116 * t103;
t97 = -t116 * t119 - t198 * t117;
t209 = t97 * pkin(7) - t48;
t175 = t116 * t102;
t99 = t198 * t103;
t208 = -t99 + t175;
t95 = t116 * t117 - t198 * t119;
t43 = t95 * pkin(7) + t208;
t237 = t240 * (-t115 * t43 - t118 * t209);
t236 = t240 * (-t115 * t209 + t118 * t43);
t151 = qJD(2) + qJD(3);
t27 = t115 * t95 + t118 * t97;
t159 = t27 * qJD(4);
t16 = t151 * t27 + t159;
t136 = t115 * t97 - t118 * t95;
t216 = t136 ^ 2 - t27 ^ 2;
t220 = t216 * qJD(1);
t210 = t136 * qJD(1) * t27;
t160 = qJD(4) * t136;
t17 = t151 * t136 + t160;
t205 = t151 * t48;
t204 = -pkin(3) / 0.2e1;
t131 = -t99 / 0.2e1;
t199 = t97 * pkin(3);
t197 = pkin(2) * t116;
t196 = t117 * pkin(2);
t191 = pkin(3) * qJD(3);
t190 = pkin(3) * qJD(4);
t112 = -t119 * pkin(2) - pkin(1);
t74 = t95 * pkin(3) + t112;
t177 = qJD(1) * t74;
t150 = t198 * pkin(2);
t111 = t150 + pkin(3);
t176 = t115 * t111;
t174 = t116 * t118;
t173 = t118 * t111;
t24 = t74 * t27;
t75 = t196 - t199;
t18 = -t136 * t75 - t24;
t171 = t18 * qJD(1);
t25 = t74 * t136;
t19 = -t27 * t75 + t25;
t170 = t19 * qJD(1);
t20 = -t136 * t199 + t24;
t169 = t20 * qJD(1);
t21 = -t199 * t27 - t25;
t168 = t21 * qJD(1);
t42 = t95 ^ 2 - t97 ^ 2;
t166 = t42 * qJD(1);
t50 = -t112 * t97 + t95 * t196;
t163 = t50 * qJD(1);
t51 = -t112 * t95 - t97 * t196;
t162 = t51 * qJD(1);
t71 = t131 + t99 / 0.2e1;
t158 = t71 * qJD(1);
t157 = qJD(1) * t112;
t156 = qJD(1) * t119;
t155 = qJD(3) * t112;
t106 = -t117 ^ 2 + t119 ^ 2;
t154 = t106 * qJD(1);
t153 = t117 * qJD(2);
t152 = t119 * qJD(2);
t149 = pkin(1) * t117 * qJD(1);
t148 = pkin(1) * t156;
t146 = t136 * t177;
t145 = t27 * t177;
t142 = t95 * t157;
t141 = t97 * t157;
t140 = t198 * t115;
t139 = t117 * t156;
t138 = t198 * qJD(2);
t137 = t198 * qJD(3);
t134 = pkin(3) * t233;
t66 = t151 * t97;
t132 = -t150 / 0.2e1;
t88 = (t140 + t174) * pkin(2);
t129 = t88 * qJD(2);
t120 = t132 + pkin(3) / 0.2e1 + t111 / 0.2e1;
t72 = t120 * t115;
t128 = t72 * qJD(2);
t83 = pkin(2) * t174 + t176;
t127 = t83 * qJD(2);
t73 = t120 * t118;
t126 = t73 * qJD(2);
t107 = t115 * t197;
t82 = t107 - t173;
t125 = t82 * qJD(2);
t89 = t118 * t150 - t107;
t124 = t89 * qJD(2);
t85 = t89 * qJD(3);
t84 = t88 * qJD(3);
t77 = t83 * qJD(4);
t76 = t82 * qJD(4);
t69 = t97 * t95 * qJD(1);
t68 = t107 - t173 / 0.2e1 + (t204 + t132) * t118;
t67 = t115 * t204 - t176 / 0.2e1 + (-t174 - t140 / 0.2e1) * pkin(2);
t65 = t151 * t95;
t49 = 0.2e1 * t131 + t175;
t1 = [0, 0, 0, t117 * t152, t106 * qJD(2), 0, 0, 0, -pkin(1) * t153, -pkin(1) * t152, t95 * t66, t151 * t42, 0, 0, 0, t50 * qJD(2) - t97 * t155, t51 * qJD(2) - t95 * t155, -t17 * t27, (qJD(4) + t151) * t216, 0, 0, 0, t18 * qJD(2) - t20 * qJD(3) - t74 * t159, t19 * qJD(2) - t21 * qJD(3) + t160 * t74; 0, 0, 0, t139, t154, t152, -t153, 0, -pkin(5) * t152 - t149, pkin(5) * t153 - t148, t69, t166, -t65, t66, 0, qJD(2) * t208 + t49 * qJD(3) + t163, t162 + t205, -t210, t220, t17, t16, 0, t171 + t236, t170 + t237; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, t166, -t65, t66, 0, t49 * qJD(2) + qJD(3) * t208 - t141, -t142 + t205, -t210, t220, t17, t16, 0, -t169 + t236, -t168 + t237; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t210, t220, t17, t16, 0, -t145 + t236, t146 + t237; 0, 0, 0, -t139, -t154, 0, 0, 0, t149, t148, -t69, -t166, 0, 0, 0, t71 * qJD(3) - t163, -t162, t210, -t220, 0, 0, 0, -t171, -t170; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3) * t197, -pkin(2) * t137, 0, 0, 0, 0, 0, -t84 - t77, -t85 + t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t151 * t197 + t158, (-t138 - t137) * pkin(2), 0, 0, 0, 0, 0, t67 * qJD(4) - t129 - t84, t68 * qJD(4) - t124 - t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67 * qJD(3) - t127 - t77, t68 * qJD(3) + t125 + t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t69, -t166, 0, 0, 0, -t71 * qJD(2) + t141, t142, t210, -t220, 0, 0, 0, t169, t168; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t197 - t158, pkin(2) * t138, 0, 0, 0, 0, 0, -t72 * qJD(4) + t129, -t73 * qJD(4) + t124; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t115 * t190, -t118 * t190; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t115 * t134 - t128, -t118 * t134 - t126; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t210, -t220, 0, 0, 0, t145, -t146; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72 * qJD(3) + t127, t73 * qJD(3) - t125; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t115 * t191 + t128, t118 * t191 + t126; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t1;
