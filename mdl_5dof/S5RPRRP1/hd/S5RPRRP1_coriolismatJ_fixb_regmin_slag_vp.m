% Calculate minimal parameter regressor of coriolis matrix for
% S5RPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x24]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 12:27
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RPRRP1_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP1_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP1_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP1_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 12:27:24
% EndTime: 2021-01-15 12:27:29
% DurationCPUTime: 1.40s
% Computational Cost: add. (2063->161), mult. (3412->183), div. (0->0), fcn. (3597->4), ass. (0->125)
t149 = qJD(3) + qJD(4);
t120 = sin(qJ(4));
t121 = sin(qJ(3));
t122 = cos(qJ(4));
t123 = cos(qJ(3));
t110 = -t120 * t121 + t122 * t123;
t124 = -pkin(1) - pkin(6);
t201 = -pkin(7) + t124;
t112 = t201 * t121;
t113 = t201 * t123;
t198 = t112 * t120 - t113 * t122;
t51 = -qJ(5) * t110 - t198;
t207 = t149 * t51;
t186 = t122 * pkin(3);
t108 = t120 * t123 + t121 * t122;
t157 = t108 * qJD(4);
t117 = pkin(4) + t186;
t190 = -t117 / 0.2e1;
t132 = t186 / 0.2e1 + t190;
t188 = t108 * pkin(4);
t36 = t188 / 0.2e1 - t132 * t108;
t206 = pkin(4) * t157 + qJD(3) * t36;
t189 = pkin(3) * t120;
t205 = t36 * qJD(4) - (-t117 * t108 + t110 * t189) * qJD(3);
t177 = t108 * qJ(5);
t105 = t122 * t112;
t174 = t120 * t113;
t199 = -t105 - t174;
t52 = t199 + t177;
t200 = (pkin(4) / 0.2e1 + t132) * t108;
t204 = qJD(3) * t200;
t203 = qJD(4) * t200;
t202 = t149 * t198;
t71 = t149 * t110;
t107 = t108 ^ 2;
t196 = t110 ^ 2;
t194 = -t52 / 0.2e1;
t6 = (t194 + t52 / 0.2e1) * t110;
t193 = t6 * qJD(3);
t192 = t52 * pkin(4);
t191 = -t105 / 0.2e1;
t187 = t110 * pkin(4);
t185 = t123 * pkin(3);
t183 = pkin(3) * qJD(4);
t182 = qJD(3) * pkin(3);
t180 = t6 * qJD(1);
t115 = pkin(3) * t121 + qJ(2);
t83 = t115 + t188;
t63 = t83 * t110;
t12 = t108 * t52 - t110 * t51;
t176 = t12 * qJD(1);
t84 = t185 + t187;
t18 = t108 * t84 + t63;
t173 = t18 * qJD(1);
t62 = t83 * t108;
t19 = t110 * t84 - t62;
t172 = t19 * qJD(1);
t20 = -t108 * t187 - t63;
t171 = t20 * qJD(1);
t21 = -pkin(4) * t196 + t62;
t170 = t21 * qJD(1);
t138 = -t120 * t108 / 0.2e1;
t31 = (t190 - pkin(4) / 0.2e1) * t110 + (t138 - t123 / 0.2e1) * pkin(3);
t169 = t31 * qJD(1);
t133 = -t107 / 0.2e1 - t196 / 0.2e1;
t33 = -0.1e1 / 0.2e1 + t133;
t168 = t33 * qJD(1);
t167 = t200 * qJD(1);
t166 = t52 * qJD(4);
t57 = t107 - t196;
t165 = t57 * qJD(1);
t60 = t108 * t185 + t110 * t115;
t164 = t60 * qJD(1);
t61 = -t108 * t115 + t110 * t185;
t163 = t61 * qJD(1);
t69 = t191 + t105 / 0.2e1;
t162 = t69 * qJD(1);
t72 = t107 + t196;
t161 = t72 * qJD(1);
t160 = t83 * qJD(1);
t159 = qJD(1) * qJ(2);
t158 = t108 * qJD(1);
t156 = t110 * qJD(1);
t155 = t110 * qJD(4);
t102 = t110 * qJD(5);
t114 = t121 ^ 2 - t123 ^ 2;
t154 = t114 * qJD(1);
t153 = t121 * qJD(1);
t152 = t121 * qJD(3);
t151 = t123 * qJD(1);
t150 = t123 * qJD(3);
t148 = t120 * t183;
t147 = t122 * t183;
t145 = pkin(4) * t156;
t143 = qJ(2) * t153;
t142 = qJ(2) * t151;
t141 = t115 * t158;
t140 = t115 * t156;
t139 = t121 * t151;
t136 = pkin(3) * t149;
t10 = t83 * t84;
t131 = qJD(1) * t10 + qJD(2) * t6;
t11 = pkin(4) * t63;
t130 = t11 * qJD(1);
t47 = 0.2e1 * t191 - t174;
t125 = t186 * t194 - t190 * t52;
t2 = -t192 / 0.2e1 + t125;
t88 = (-t117 + t186) * t189;
t126 = -qJD(1) * t2 - qJD(2) * t200 - qJD(3) * t88;
t100 = t110 * qJD(2);
t97 = t108 * qJD(2);
t96 = t108 * qJD(5);
t79 = t108 * t156;
t70 = t149 * t108;
t67 = t69 * qJD(3);
t66 = t69 * qJD(4);
t59 = t120 * t182 - t162;
t58 = t122 * t182;
t49 = -t120 * t136 + t162;
t48 = t122 * t136;
t32 = 0.1e1 / 0.2e1 + t133;
t30 = pkin(3) * t138 + t110 * t190 + t185 / 0.2e1 + t187 / 0.2e1;
t23 = t47 + t177;
t1 = t192 / 0.2e1 + t125;
t3 = [0, 0, 0, 0, qJD(2), qJ(2) * qJD(2), -t121 * t150, t114 * qJD(3), 0, 0, 0, qJ(2) * t150 + qJD(2) * t121, -qJ(2) * t152 + qJD(2) * t123, -t108 * t71, t149 * t57, 0, 0, 0, qJD(3) * t60 + t115 * t155 + t97, qJD(3) * t61 - t115 * t157 + t100, qJD(3) * t18 - qJD(4) * t20 + t97, qJD(3) * t19 - qJD(4) * t21 + t100, t72 * qJD(5), qJD(2) * t83 + qJD(3) * t10 + qJD(4) * t11 + qJD(5) * t12; 0, 0, 0, 0, qJD(1), t159, 0, 0, 0, 0, 0, t153, t151, 0, 0, 0, 0, 0, t158, t156, t158, t156, 0, qJD(5) * t32 + t160 + t193; 0, 0, 0, 0, 0, 0, -t139, t154, -t152, -t150, 0, -t124 * t152 + t142, -t124 * t150 - t143, -t79, t165, -t70, -t71, 0, qJD(3) * t199 + qJD(4) * t47 + t164, t163 + t202, qJD(3) * t52 + qJD(4) * t23 + t173, t172 - t207, t205, (t117 * t52 + t189 * t51) * qJD(3) + t1 * qJD(4) + t30 * qJD(5) + t131; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t79, t165, -t70, -t71, 0, qJD(3) * t47 + qJD(4) * t199 + t140, -t141 + t202, qJD(3) * t23 + t166 - t171, -t170 - t207, t206, pkin(4) * t166 + qJD(3) * t1 + t130; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t161, qJD(2) * t32 + qJD(3) * t30 + t176; 0, 0, 0, 0, -qJD(1), -t159, 0, 0, 0, 0, 0, -t153, -t151, 0, 0, 0, 0, 0, -t158, -t156, -t158, -t156, 0, qJD(5) * t33 - t160 + t193; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t152, -t150, 0, 0, 0, 0, 0, -t70, -t71, -t70, -t71, 0, t180 - t205; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t70, -t71, -t70, -t71, 0, -t206; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t168; 0, 0, 0, 0, 0, 0, t139, -t154, 0, 0, 0, -t142, t143, t79, -t165, 0, 0, 0, t66 - t164, -t163, -t102 + t66 - t173, t96 - t172, -t203, qJD(4) * t2 + qJD(5) * t31 - t131; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t180 + t203; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t148, -t147, -t148, -t147, 0, t88 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, -t48, t49, -t48, -t167, -pkin(4) * t148 - t126; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t156, t158, 0, t169; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79, -t165, 0, 0, 0, -t67 - t140, t141, -t102 - t67 + t171, t96 + t170, t204, -pkin(4) * t102 - qJD(3) * t2 - t130; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t204; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, t58, t59, t58, t167, t126; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t156, t158, 0, -t145; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71, -t70, -t161, pkin(4) * t155 - qJD(2) * t33 - qJD(3) * t31 - t176; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t168; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t156, -t158, 0, -t169; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t156, -t158, 0, t145; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t3;
