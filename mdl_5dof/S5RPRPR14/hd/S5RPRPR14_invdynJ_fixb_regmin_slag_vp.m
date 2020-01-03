% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPRPR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% 
% Output:
% tau_reg [5x22]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRPR14_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR14_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR14_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR14_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR14_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR14_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:35:19
% EndTime: 2019-12-31 18:35:23
% DurationCPUTime: 1.52s
% Computational Cost: add. (1558->241), mult. (3127->324), div. (0->0), fcn. (2122->10), ass. (0->139)
t103 = sin(qJ(3));
t106 = cos(qJ(3));
t172 = sin(pkin(8));
t173 = cos(pkin(8));
t57 = -t172 * t103 + t173 * t106;
t115 = t173 * t103 + t172 * t106;
t189 = t115 * qJD(1);
t192 = qJD(5) + t189;
t102 = sin(qJ(5));
t105 = cos(qJ(5));
t155 = t105 * qJD(3);
t54 = t57 * qJD(1);
t35 = t102 * t54 - t155;
t194 = t192 * t35;
t91 = t103 * pkin(3);
t175 = qJ(2) + t91;
t135 = t105 * t192;
t131 = qJDD(1) * t172;
t132 = qJDD(1) * t173;
t31 = -t54 * qJD(3) - t103 * t132 - t106 * t131;
t29 = -qJDD(5) + t31;
t176 = t102 * t29;
t193 = t135 * t192 - t176;
t110 = qJD(1) ^ 2;
t104 = sin(qJ(1));
t107 = cos(qJ(1));
t191 = g(1) * t104 - g(2) * t107;
t117 = -t110 * qJ(2) - t191;
t125 = g(1) * t107 + g(2) * t104;
t97 = qJD(1) * qJD(2);
t149 = 0.2e1 * t97;
t96 = qJDD(1) * qJ(2);
t188 = -t125 + 0.2e1 * t96 + t149;
t150 = t106 * qJDD(1);
t153 = qJD(1) * qJD(4);
t154 = qJD(1) * qJD(3);
t159 = qJD(3) * t103;
t108 = -pkin(1) - pkin(6);
t64 = t108 * qJDD(1) + qJDD(2);
t59 = t106 * t64;
t65 = t108 * qJD(1) + qJD(2);
t15 = -t106 * t153 - t65 * t159 + qJDD(3) * pkin(3) + t59 + (t103 * t154 - t150) * qJ(4);
t133 = -qJ(4) * qJD(1) + t65;
t158 = qJD(3) * t106;
t21 = t133 * t158 + (-qJ(4) * qJDD(1) - t153 + t64) * t103;
t5 = t173 * t15 - t172 * t21;
t1 = -qJDD(3) * pkin(4) - t5;
t44 = t133 * t103;
t146 = t172 * t44;
t160 = qJD(1) * t106;
t45 = -qJ(4) * t160 + t106 * t65;
t42 = qJD(3) * pkin(3) + t45;
t16 = t173 * t42 - t146;
t12 = -qJD(3) * pkin(4) - t16;
t61 = t175 * qJD(1) + qJD(4);
t18 = pkin(4) * t189 - t54 * pkin(7) + t61;
t6 = t172 * t15 + t173 * t21;
t143 = qJDD(3) * pkin(7) + qJD(5) * t18 + t6;
t164 = qJ(4) - t108;
t116 = -t106 * qJD(4) + t164 * t159;
t134 = t164 * t106;
t43 = -qJD(3) * t134 - t103 * qJD(4);
t20 = t172 * t116 + t173 * t43;
t30 = pkin(4) * t115 - t57 * pkin(7) + t175;
t62 = t164 * t103;
t34 = -t172 * t134 - t173 * t62;
t136 = qJD(3) * t172;
t137 = qJD(3) * t173;
t53 = -t103 * t137 - t106 * t136;
t187 = t1 * t57 + t12 * t53 - (qJD(5) * t30 + t20) * t192 - t143 * t115 + t34 * t29;
t77 = t172 * pkin(3) + pkin(7);
t95 = qJ(3) + pkin(8);
t83 = sin(t95);
t84 = cos(t95);
t186 = t191 * t84 + (pkin(3) * t160 + t54 * pkin(4) + pkin(7) * t189 + qJD(5) * t77) * t192 - g(3) * t83 + t1;
t184 = g(3) * t103;
t183 = t12 * t57;
t182 = t30 * t29;
t37 = t102 * qJD(3) + t105 * t54;
t181 = t37 * t54;
t180 = t54 * t35;
t157 = qJD(5) * t102;
t32 = qJD(3) * t189 + t103 * t131 - t106 * t132;
t9 = qJD(5) * t155 + t102 * qJDD(3) - t105 * t32 - t54 * t157;
t179 = t9 * t102;
t40 = t173 * t44;
t17 = t172 * t42 + t40;
t178 = t107 * pkin(1) + t104 * qJ(2);
t26 = t105 * t29;
t100 = t106 ^ 2;
t174 = -t103 ^ 2 + t100;
t171 = pkin(1) * qJDD(1);
t170 = t104 * t102;
t169 = t104 * t105;
t168 = t107 * t102;
t167 = t107 * t105;
t165 = t61 * qJD(1);
t163 = pkin(3) * t158 + qJD(2);
t109 = qJD(3) ^ 2;
t161 = -t109 - t110;
t156 = qJD(5) * t105;
t152 = qJDD(3) * t103;
t151 = t103 * qJDD(1);
t147 = t57 * t157;
t145 = t106 * t154;
t13 = qJD(3) * pkin(7) + t17;
t122 = qJDD(4) + t96 + t97 + (t145 + t151) * pkin(3);
t8 = -t31 * pkin(4) + t32 * pkin(7) + t122;
t144 = qJD(5) * t13 - t8;
t142 = -t105 * qJDD(3) - t102 * t32;
t130 = -qJD(5) * t115 - qJD(1);
t129 = qJDD(2) - t171;
t128 = t53 * t37 + t57 * t9;
t123 = -t192 * t53 + t29 * t57;
t121 = -t26 + (-t102 * t189 - t157) * t192;
t120 = g(3) * t84 - t143;
t119 = 0.2e1 * qJ(2) * t154 + qJDD(3) * t108;
t23 = t173 * t45 - t146;
t114 = t77 * t29 + (t12 + t23) * t192;
t52 = t103 * t136 - t106 * t137;
t112 = t115 * t6 + t16 * t53 - t17 * t52 + t5 * t57 - t191;
t111 = -t108 * t109 + t188;
t101 = -qJ(4) - pkin(6);
t90 = t107 * qJ(2);
t87 = qJDD(3) * t106;
t78 = -t173 * pkin(3) - pkin(4);
t50 = t83 * t167 - t170;
t49 = t83 * t168 + t169;
t48 = t83 * t169 + t168;
t47 = -t83 * t170 + t167;
t33 = t173 * t134 - t172 * t62;
t24 = -t52 * pkin(4) - t53 * pkin(7) + t163;
t22 = t172 * t45 + t40;
t19 = -t173 * t116 + t172 * t43;
t10 = t37 * qJD(5) + t142;
t7 = t105 * t8;
t4 = t102 * t18 + t105 * t13;
t3 = -t102 * t13 + t105 * t18;
t2 = [qJDD(1), t191, t125, qJDD(2) - 0.2e1 * t171 - t191, t188, -t129 * pkin(1) - g(1) * (-t104 * pkin(1) + t90) - g(2) * t178 + (t149 + t96) * qJ(2), t100 * qJDD(1) - 0.2e1 * t103 * t145, -0.2e1 * t103 * t150 - 0.2e1 * t174 * t154, -t109 * t103 + t87, -t109 * t106 - t152, 0, t111 * t103 + t119 * t106, -t119 * t103 + t111 * t106, -t189 * t20 + t19 * t54 + t34 * t31 - t33 * t32 - t112, t6 * t34 + t17 * t20 - t5 * t33 - t16 * t19 + t122 * t175 + t61 * t163 - g(1) * (t107 * t91 + t90 + (-pkin(1) + t101) * t104) - g(2) * (-t107 * t101 + t104 * t91 + t178), t128 * t105 - t37 * t147, (-t102 * t37 - t105 * t35) * t53 + (-t10 * t105 - t179 + (t102 * t35 - t105 * t37) * qJD(5)) * t57, -t123 * t105 + t115 * t9 - t147 * t192 - t37 * t52, -t156 * t192 * t57 - t10 * t115 + t123 * t102 + t35 * t52, -t115 * t29 - t192 * t52, -g(1) * t50 - g(2) * t48 + t33 * t10 + t19 * t35 - t3 * t52 + t7 * t115 + (t24 * t192 - t182 + (-t115 * t13 - t192 * t34 + t183) * qJD(5)) * t105 + t187 * t102, g(1) * t49 - g(2) * t47 + t19 * t37 + t33 * t9 + t4 * t52 + (-(-qJD(5) * t34 + t24) * t192 + t182 + t144 * t115 - qJD(5) * t183) * t102 + t187 * t105; 0, 0, 0, qJDD(1), -t110, t129 + t117, 0, 0, 0, 0, 0, t161 * t103 + t87, t161 * t106 - t152, t115 * t31 + t189 * t52 + t57 * t32 - t53 * t54, t112 - t165, 0, 0, 0, 0, 0, t115 * t176 - t57 * t10 - t53 * t35 + (t102 * t52 + t130 * t105) * t192, t115 * t26 + (-t130 * t102 + t105 * t52) * t192 - t128; 0, 0, 0, 0, 0, 0, t106 * t110 * t103, t174 * t110, t150, -t151, qJDD(3), t117 * t106 + t184 + t59, g(3) * t106 + (-t117 - t64) * t103, (t17 - t22) * t54 - (t16 - t23) * t189 + (t172 * t31 + t173 * t32) * pkin(3), t16 * t22 - t17 * t23 + (t173 * t5 + t172 * t6 + t184 + (-t191 - t165) * t106) * pkin(3), t37 * t135 + t179, (t9 - t194) * t105 + (-t192 * t37 - t10) * t102, -t181 + t193, t121 + t180, -t192 * t54, t78 * t10 + t114 * t102 - t186 * t105 - t22 * t35 - t3 * t54, t186 * t102 + t114 * t105 - t22 * t37 + t4 * t54 + t78 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t189 ^ 2 - t54 ^ 2, t16 * t54 + t17 * t189 + t122 - t125, 0, 0, 0, 0, 0, t121 - t180, -t181 - t193; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37 * t35, -t35 ^ 2 + t37 ^ 2, t9 + t194, -t142 + (-qJD(5) + t192) * t37, -t29, -g(1) * t47 - g(2) * t49 + t120 * t102 - t12 * t37 - t13 * t156 + t192 * t4 + t7, g(1) * t48 - g(2) * t50 + t102 * t144 + t105 * t120 + t12 * t35 + t192 * t3;];
tau_reg = t2;
