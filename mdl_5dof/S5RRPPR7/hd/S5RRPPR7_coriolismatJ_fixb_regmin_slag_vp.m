% Calculate minimal parameter regressor of coriolis matrix for
% S5RRPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x23]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RRPPR7_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR7_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR7_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR7_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:36:25
% EndTime: 2019-12-31 19:36:30
% DurationCPUTime: 1.64s
% Computational Cost: add. (1661->191), mult. (3375->281), div. (0->0), fcn. (3621->6), ass. (0->164)
t111 = cos(qJ(5));
t108 = sin(pkin(8));
t110 = sin(qJ(2));
t212 = -qJ(3) - pkin(6);
t94 = t212 * t110;
t206 = t108 * t94;
t197 = cos(pkin(8));
t112 = cos(qJ(2));
t95 = t212 * t112;
t93 = t197 * t95;
t225 = -t93 + t206;
t102 = t197 * t112;
t195 = t108 * t110;
t88 = -t102 + t195;
t229 = -t88 * pkin(4) + t225;
t202 = t111 * t229;
t109 = sin(qJ(5));
t57 = -t108 * t95 - t197 * t94;
t90 = t108 * t112 + t110 * t197;
t125 = pkin(4) * t90 + t57;
t231 = t109 * t125;
t160 = qJD(5) * t109;
t46 = t109 * t90;
t182 = t46 * qJD(1);
t230 = t182 + t160;
t219 = t90 ^ 2;
t85 = t88 ^ 2;
t226 = t85 + t219;
t228 = qJD(3) * t226;
t227 = t226 * qJD(1);
t49 = t111 * t88;
t137 = 0.2e1 * t109 * t49;
t106 = t109 ^ 2;
t107 = t111 ^ 2;
t98 = t106 - t107;
t116 = qJD(1) * t137 - qJD(2) * t98;
t130 = -t225 * t88 + t57 * t90;
t223 = qJD(3) * t130;
t221 = t130 * qJD(1);
t220 = t125 * t88;
t218 = -t90 / 0.2e1;
t134 = -t93 / 0.2e1;
t217 = pkin(3) + pkin(7);
t215 = t90 * pkin(3);
t214 = -t109 / 0.2e1;
t105 = t110 * pkin(2);
t210 = qJD(2) * pkin(2);
t141 = -t112 * pkin(2) - pkin(1);
t126 = -t90 * qJ(4) + t141;
t24 = t217 * t88 + t126;
t13 = t109 * t24 - t111 * t125;
t199 = t88 * qJ(4);
t136 = t105 + t199;
t25 = t217 * t90 + t136;
t1 = t111 * t220 + t13 * t88 - t25 * t46;
t209 = t1 * qJD(1);
t101 = pkin(2) * t108 + qJ(4);
t208 = t101 * t88;
t207 = t101 * t90;
t205 = t111 * t25;
t14 = t111 * t24 + t231;
t2 = -t109 * t220 + t14 * t88 - t90 * t205;
t204 = t2 * qJD(1);
t203 = t229 * t109;
t45 = pkin(3) * t88 + t126;
t51 = t136 + t215;
t7 = t45 * t51;
t201 = t7 * qJD(1);
t8 = t13 * t90 + t202 * t88;
t200 = t8 * qJD(1);
t9 = -t14 * t90 + t203 * t88;
t198 = t9 * qJD(1);
t12 = t141 * t105;
t193 = t12 * qJD(1);
t15 = -t45 * t90 - t51 * t88;
t192 = t15 * qJD(1);
t16 = t45 * t88 - t51 * t90;
t191 = t16 * qJD(1);
t103 = -pkin(2) * t197 - pkin(3);
t104 = t105 / 0.2e1;
t19 = t104 + (pkin(3) / 0.2e1 - t103 / 0.2e1) * t90 + (qJ(4) / 0.2e1 + t101 / 0.2e1) * t88;
t188 = t19 * qJD(1);
t21 = t226 * t109;
t187 = t21 * qJD(1);
t154 = t219 - t85;
t26 = t154 * t109;
t186 = t26 * qJD(1);
t27 = t154 * t111;
t185 = t27 * qJD(1);
t28 = t226 * t111;
t184 = t28 * qJD(1);
t113 = (-t108 * t88 / 0.2e1 + t197 * t218) * pkin(2);
t35 = -t105 / 0.2e1 + t113;
t183 = t35 * qJD(1);
t181 = t46 * qJD(5);
t180 = t49 * qJD(1);
t179 = t49 * qJD(2);
t82 = t195 / 0.2e1 - t102 / 0.2e1;
t176 = t82 * qJD(1);
t175 = t219 * qJD(1);
t174 = t88 * qJD(1);
t173 = t88 * qJD(2);
t172 = t88 * qJD(3);
t171 = t88 * qJD(4);
t170 = t90 * qJD(1);
t169 = t90 * qJD(2);
t79 = t90 * qJD(3);
t168 = t90 * qJD(4);
t99 = -t110 ^ 2 + t112 ^ 2;
t166 = t99 * qJD(1);
t165 = qJD(1) * t109;
t164 = qJD(1) * t111;
t163 = qJD(1) * t112;
t162 = qJD(4) * t109;
t161 = qJD(4) * t111;
t159 = qJD(5) * t111;
t158 = t109 * qJD(2);
t157 = t110 * qJD(2);
t156 = t111 * qJD(2);
t155 = t112 * qJD(2);
t153 = pkin(1) * t110 * qJD(1);
t152 = pkin(1) * t163;
t151 = t45 * t170;
t150 = t88 * t170;
t149 = t88 * t169;
t147 = t106 * t174;
t146 = t90 * t160;
t145 = t90 * t159;
t144 = t88 * t165;
t143 = t88 * t158;
t142 = t219 * t165;
t75 = t90 * t164;
t140 = t110 * t163;
t139 = t109 * t159;
t138 = t109 * t156;
t135 = qJD(5) + t170;
t132 = qJD(2) * t137;
t100 = -pkin(7) + t103;
t129 = t100 * t88 + t207;
t128 = t135 * t109;
t55 = t134 + t93 / 0.2e1;
t127 = qJD(1) * t55 + qJD(2) * t101;
t124 = t100 * t218 + t208 / 0.2e1;
t123 = qJD(5) * t82 + t150;
t122 = t88 * t128;
t114 = t25 / 0.2e1 + t124;
t5 = t114 * t111;
t121 = -qJD(1) * t5 + t101 * t158;
t3 = t114 * t109;
t120 = -qJD(1) * t3 - t101 * t156;
t48 = (t107 / 0.2e1 - t106 / 0.2e1) * t88;
t119 = -qJD(1) * t48 + t138;
t118 = t109 * t164 * t85 + qJD(2) * t48;
t52 = t98 * t85;
t117 = qJD(1) * t52 + t132;
t76 = t82 * qJD(2);
t64 = -t75 - t159;
t43 = t48 * qJD(5);
t39 = 0.2e1 * t134 + t206;
t34 = t104 + t113;
t20 = -t208 / 0.2e1 + t103 * t90 / 0.2e1 + t104 + t199 / 0.2e1 + t215 / 0.2e1;
t6 = -t203 / 0.2e1 - t205 / 0.2e1 + t229 * t214 + t124 * t111;
t4 = t124 * t109 + t25 * t214 + t202;
t10 = [0, 0, 0, t110 * t155, t99 * qJD(2), 0, 0, 0, -pkin(1) * t157, -pkin(1) * t155, t228, qJD(2) * t12 + t223, t228, qJD(2) * t15 + t168 * t88, qJD(2) * t16 + qJD(4) * t219, qJD(2) * t7 - t168 * t45 + t223, t106 * t149 + t139 * t85, -qJD(5) * t52 + t132 * t90, qJD(2) * t26 + t145 * t88, qJD(2) * t27 - t146 * t88, -t149, qJD(2) * t1 + qJD(3) * t28 + qJD(5) * t9 + t162 * t219, qJD(2) * t2 - qJD(3) * t21 + qJD(5) * t8 + t161 * t219; 0, 0, 0, t140, t166, t155, -t157, 0, -pkin(6) * t155 - t153, pkin(6) * t157 - t152, (-t108 * t90 + t197 * t88) * t210, t193 + (-t108 * t57 - t197 * t225) * t210 + t34 * qJD(3), (-t103 * t88 - t207) * qJD(2) - t171, qJD(2) * t225 + t192, -qJD(2) * t57 + t191, t201 + (-t101 * t57 + t103 * t225) * qJD(2) + t20 * qJD(3) + t39 * qJD(4), t43 + (t138 + t147) * t90, t116 * t90 - 0.2e1 * t139 * t88, -t156 * t88 + t186, t143 + t185, -t123, t209 + (-t111 * t129 - t231) * qJD(2) - t49 * qJD(4) + t4 * qJD(5), -t125 * t156 + t204 + t6 * qJD(5) + (qJD(2) * t129 + t171) * t109; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t227, qJD(2) * t34 + t221, t227, 0, 0, qJD(2) * t20 + t221, 0, 0, 0, 0, 0, t184, -t187; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t173, t150, t175, qJD(2) * t39 - t151, 0, 0, 0, 0, 0, t142 - t179, t164 * t219 + t143; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t118, -t117, t135 * t49, -t122, -t76, qJD(2) * t4 - qJD(5) * t14 + t198, qJD(2) * t6 + qJD(5) * t13 + t200; 0, 0, 0, -t140, -t166, 0, 0, 0, t153, t152, 0, qJD(3) * t35 - t193, 0, t79 - t192, -t172 - t191, -qJD(3) * t19 + qJD(4) * t55 - t201, -t147 * t90 + t43, -0.2e1 * t111 * t122, -t146 - t186, -t145 - t185, t123, qJD(5) * t3 - t109 * t172 - t209, -qJD(3) * t49 + qJD(5) * t5 - t204; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4), t101 * qJD(4), -t139, t98 * qJD(5), 0, 0, 0, t101 * t159 + t162, -t101 * t160 + t161; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t183, 0, t170, -t174, -t188, 0, 0, 0, 0, 0, -t144, -t180; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2), t127, 0, 0, 0, 0, 0, t158, t156; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t119, -t116, -t128, t64, t176, -t100 * t160 - t120, -t100 * t159 - t121; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t227, -qJD(2) * t35 - t221, -t227, -t169, t173, qJD(2) * t19 - t168 - t221, 0, 0, 0, 0, 0, t143 - t145 - t184, t179 + t181 + t187; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t183, 0, -t170, t174, t188, 0, 0, 0, 0, 0, t144, t180; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t170, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64, t230; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t150, -t175, -qJD(2) * t55 + t151 + t79, 0, 0, 0, 0, 0, -t142 - t181, (-qJD(5) * t90 - t175) * t111; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(2), -t127, 0, 0, 0, 0, 0, -t158, -t156; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t170, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t230, -t135 * t111; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t118, t117, (-t164 * t88 + t158) * t90, (t144 + t156) * t90, -t76, -qJD(2) * t3 + qJD(4) * t46 + t111 * t79 - t198, -qJD(2) * t5 - qJD(3) * t46 + t161 * t90 - t200; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t119, t116, t90 * t165, t75, -t176, t120, t121; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, -t182; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t182, t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t10;
