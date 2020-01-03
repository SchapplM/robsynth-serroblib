% Calculate minimal parameter regressor of coriolis matrix for
% S5RPRPR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x22]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RPRPR14_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR14_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR14_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR14_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:35:19
% EndTime: 2019-12-31 18:35:23
% DurationCPUTime: 1.45s
% Computational Cost: add. (1743->160), mult. (3322->251), div. (0->0), fcn. (3576->6), ass. (0->140)
t126 = sin(qJ(5));
t127 = sin(qJ(3));
t130 = -pkin(1) - pkin(6);
t187 = -qJ(4) + t130;
t110 = t187 * t127;
t129 = cos(qJ(3));
t153 = t187 * t129;
t203 = sin(pkin(8));
t204 = cos(pkin(8));
t131 = t204 * t110 + t203 * t153;
t210 = t126 * t131;
t128 = cos(qJ(5));
t208 = t128 * t131;
t115 = t203 * t127;
t154 = t204 * t129;
t107 = t154 - t115;
t108 = -t204 * t127 - t203 * t129;
t213 = qJD(3) * pkin(3);
t221 = (t203 * t107 + t204 * t108) * t213;
t105 = t107 ^ 2;
t106 = t108 ^ 2;
t76 = t105 + t106;
t149 = t105 / 0.2e1 + t106 / 0.2e1;
t124 = t126 ^ 2;
t125 = t128 ^ 2;
t114 = t125 - t124;
t183 = qJD(1) * t128;
t157 = t107 * t183;
t220 = qJD(3) * t114 - 0.2e1 * t126 * t157;
t219 = t108 / 0.2e1;
t218 = -t126 / 0.2e1;
t217 = -t128 / 0.2e1;
t216 = t128 / 0.2e1;
t215 = t129 * pkin(3);
t120 = pkin(3) * t127 + qJ(2);
t134 = -pkin(4) * t108 - pkin(7) * t107 + t120;
t18 = -t128 * t134 + t210;
t60 = pkin(4) * t107 - pkin(7) * t108 + t215;
t209 = t128 * t60;
t1 = -t209 * t108 + (-t18 + t210) * t107;
t212 = t1 * qJD(1);
t211 = t126 * t60;
t19 = t126 * t134 + t208;
t2 = t211 * t108 + (-t19 + t208) * t107;
t207 = t2 * qJD(1);
t49 = t126 * t107;
t72 = t203 * t110 - t204 * t153;
t10 = -t108 * t18 - t72 * t49;
t202 = qJD(1) * t10;
t55 = t128 * t107;
t11 = t108 * t19 + t72 * t55;
t201 = qJD(1) * t11;
t17 = t107 * t72 + t108 * t131;
t200 = qJD(1) * t17;
t25 = -0.1e1 / 0.2e1 - t149;
t20 = t25 * t126;
t199 = qJD(1) * t20;
t27 = t25 * t128;
t198 = qJD(1) * t27;
t168 = t105 - t106;
t29 = t168 * t126;
t197 = qJD(1) * t29;
t30 = t76 * t126;
t196 = qJD(1) * t30;
t31 = t168 * t128;
t195 = qJD(1) * t31;
t62 = t76 * t128;
t194 = qJD(1) * t62;
t193 = t25 * qJD(1);
t132 = t203 * t219 - t204 * t107 / 0.2e1;
t38 = (-t129 / 0.2e1 + t132) * pkin(3);
t192 = t38 * qJD(1);
t191 = t49 * qJD(1);
t50 = t126 * t108;
t190 = t50 * qJD(1);
t53 = t128 * t108;
t189 = t53 * qJD(1);
t188 = t76 * qJD(1);
t186 = qJD(1) * qJ(2);
t185 = qJD(1) * t108;
t184 = qJD(1) * t120;
t182 = qJD(2) * t108;
t181 = qJD(3) * t126;
t180 = qJD(3) * t128;
t179 = qJD(4) * t128;
t178 = qJD(5) * t126;
t177 = qJD(5) * t128;
t102 = t115 / 0.2e1 - t154 / 0.2e1;
t176 = t102 * qJD(1);
t113 = t127 ^ 2 - t129 ^ 2;
t175 = t113 * qJD(1);
t174 = t127 * qJD(1);
t173 = t127 * qJD(3);
t172 = t129 * qJD(1);
t171 = t129 * qJD(3);
t166 = qJ(2) * t174;
t165 = qJ(2) * t172;
t164 = qJD(1) * t107 * t125;
t163 = t126 * t185;
t162 = t126 * t180;
t161 = t108 * t178;
t160 = t108 * t177;
t159 = t108 * t107 * qJD(3);
t158 = t126 * t177;
t100 = t107 * t180;
t156 = t108 * t183;
t155 = t127 * t172;
t152 = -qJD(5) + t185;
t151 = t126 * t100;
t148 = -t49 * qJD(3) + t160;
t12 = t120 * t215;
t147 = t12 * qJD(1);
t118 = t203 * pkin(3) + pkin(7);
t119 = -t204 * pkin(3) - pkin(4);
t146 = -t107 * t118 + t108 * t119;
t145 = t152 * t128;
t142 = t118 * t219 + t119 * t107 / 0.2e1;
t135 = -t60 / 0.2e1 + t142;
t8 = t135 * t128;
t144 = -qJD(1) * t8 - t119 * t181;
t6 = t135 * t126;
t143 = qJD(1) * t6 - t119 * t180;
t47 = (t124 / 0.2e1 - t125 / 0.2e1) * t107;
t141 = -qJD(1) * t47 + t162;
t140 = t107 * t145;
t139 = qJD(5) * t102 + t107 * t185;
t138 = t105 * t126 * t183 + qJD(3) * t47;
t61 = t114 * t105;
t137 = qJD(1) * t61 + 0.2e1 * t151;
t101 = t102 * qJD(3);
t44 = t50 * qJD(5);
t42 = t47 * qJD(5);
t37 = t215 / 0.2e1 + t132 * pkin(3);
t36 = -t178 + t190;
t28 = t76 * t217 + t216;
t24 = 0.1e1 / 0.2e1 - t149;
t21 = t149 * t126 + t218;
t9 = t209 / 0.2e1 + t142 * t128 + (t126 / 0.2e1 - t218) * t72;
t7 = -t211 / 0.2e1 - t142 * t126 + (t216 - t217) * t72;
t3 = [0, 0, 0, 0, qJD(2), qJ(2) * qJD(2), -t127 * t171, t113 * qJD(3), 0, 0, 0, qJ(2) * t171 + qJD(2) * t127, -qJ(2) * t173 + qJD(2) * t129, qJD(4) * t76, qJD(2) * t120 + qJD(3) * t12 + qJD(4) * t17, -t105 * t158 + t125 * t159, -qJD(5) * t61 - 0.2e1 * t108 * t151, qJD(3) * t31 + t107 * t161, -qJD(3) * t29 + t107 * t160, -t159, qJD(3) * t1 + qJD(4) * t30 + qJD(5) * t11 - t128 * t182, qJD(3) * t2 + qJD(4) * t62 + qJD(5) * t10 + t126 * t182; 0, 0, 0, 0, qJD(1), t186, 0, 0, 0, 0, 0, t174, t172, 0, qJD(4) * t24 + t184, 0, 0, 0, 0, 0, qJD(5) * t28 - t156, qJD(5) * t21 + t163; 0, 0, 0, 0, 0, 0, -t155, t175, -t173, -t171, 0, -t130 * t173 + t165, -t130 * t171 - t166, -t221, (-t131 * t204 - t203 * t72) * t213 + t37 * qJD(4) + t147, -t42 + (t162 + t164) * t108, -0.2e1 * t107 * t158 + t220 * t108, t107 * t181 + t195, t100 - t197, -t139, t212 + (t126 * t146 - t208) * qJD(3) + t9 * qJD(5), t207 + (t128 * t146 + t210) * qJD(3) + t7 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t188, qJD(2) * t24 + qJD(3) * t37 + t200, 0, 0, 0, 0, 0, t196, t194; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t138, -t137, t152 * t49, t140, -t101, qJD(2) * t28 + qJD(3) * t9 - qJD(5) * t19 + t201, qJD(2) * t21 + qJD(3) * t7 + qJD(5) * t18 + t202; 0, 0, 0, 0, -qJD(1), -t186, 0, 0, 0, 0, 0, -t174, -t172, 0, qJD(4) * t25 - t184, 0, 0, 0, 0, 0, qJD(5) * t27 + t156, -qJD(5) * t20 - t163; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t173, -t171, 0, t221, 0, 0, 0, 0, 0, -qJD(5) * t49 + t108 * t180, -qJD(5) * t55 - t108 * t181; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t193, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t148 + t198, -qJD(3) * t55 - t161 - t199; 0, 0, 0, 0, 0, 0, t155, -t175, 0, 0, 0, -t165, t166, 0, qJD(4) * t38 - t147, -t108 * t164 - t42, 0.2e1 * t126 * t140, -qJD(5) * t53 - t195, t44 + t197, t139, qJD(5) * t8 - t107 * t179 - t212, qJD(4) * t49 - qJD(5) * t6 - t207; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t158, t114 * qJD(5), 0, 0, 0, t119 * t178, t119 * t177; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t192, 0, 0, 0, 0, 0, -t157, t191; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t141, t220, t177 - t189, t36, t176, -t118 * t177 - t144, t118 * t178 - t143; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t188, -qJD(2) * t25 - qJD(3) * t38 - t200, 0, 0, 0, 0, 0, t100 + t44 - t196, t148 - t194; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t193, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t192, 0, 0, 0, 0, 0, t157, -t191; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, t145; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t138, t137, qJD(3) * t53 - t107 * t163, -qJD(3) * t50 - t107 * t156, -t101, -qJD(2) * t27 - qJD(3) * t8 - qJD(4) * t50 - t201, qJD(2) * t20 + qJD(3) * t6 - t108 * t179 - t202; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t198, t199; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t141, -t220, t189, -t190, -t176, t144, t143; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t190, -t156; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t3;
