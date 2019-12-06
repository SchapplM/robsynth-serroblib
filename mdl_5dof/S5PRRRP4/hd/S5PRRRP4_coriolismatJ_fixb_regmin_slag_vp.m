% Calculate minimal parameter regressor of coriolis matrix for
% S5PRRRP4
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
% cmat_reg [(5*%NQJ)%x18]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:46
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5PRRRP4_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP4_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP4_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP4_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:46:22
% EndTime: 2019-12-05 16:46:26
% DurationCPUTime: 1.05s
% Computational Cost: add. (1084->161), mult. (2428->209), div. (0->0), fcn. (2182->6), ass. (0->141)
t116 = cos(qJ(4));
t172 = t116 * qJ(5);
t113 = sin(qJ(4));
t189 = t113 * pkin(4);
t90 = -t172 + t189;
t174 = t90 * t116;
t117 = cos(qJ(3));
t187 = t117 * pkin(2);
t129 = -t116 * pkin(4) - t113 * qJ(5);
t88 = -pkin(3) + t129;
t81 = t88 - t187;
t70 = t81 * t113;
t83 = t88 * t113;
t182 = t70 / 0.2e1 + t83 / 0.2e1;
t200 = t182 - t174;
t111 = t113 ^ 2;
t112 = t116 ^ 2;
t165 = t111 + t112;
t114 = sin(qJ(3));
t115 = sin(qJ(2));
t191 = cos(qJ(2));
t84 = t114 * t115 - t117 * t191;
t138 = t165 * t84;
t199 = pkin(7) * t138;
t106 = -pkin(3) - t187;
t148 = pkin(3) / 0.2e1 - t106 / 0.2e1;
t198 = t148 * t113;
t101 = t112 - t111;
t159 = qJD(2) + qJD(3);
t197 = t159 * t101;
t123 = t172 / 0.2e1 - t189 / 0.2e1;
t118 = t123 * t187;
t193 = t81 / 0.2e1;
t147 = t88 / 0.2e1 + t193;
t12 = -t147 * t90 + t118;
t18 = (-t90 / 0.2e1 - t123) * t84;
t171 = t18 * qJD(1);
t183 = t90 * t88;
t196 = t12 * qJD(2) - qJD(3) * t183 + t171;
t186 = t81 * t90;
t195 = -qJD(2) * t186 + t171;
t110 = t116 * qJD(4);
t37 = t84 * t113;
t85 = t114 * t191 + t117 * t115;
t3 = -t85 * t110 + t159 * t37;
t161 = t113 * qJD(4);
t39 = t84 * t116;
t194 = t159 * t39 + t85 * t161;
t192 = t84 / 0.2e1;
t190 = pkin(3) * t116;
t188 = t114 * pkin(2);
t185 = t85 * t81;
t184 = t85 * t88;
t180 = pkin(2) * qJD(3);
t179 = qJD(2) * pkin(2);
t6 = (0.1e1 - t165) * t85 * t84;
t178 = t6 * qJD(1);
t177 = t81 * t116;
t176 = t88 * t116;
t175 = t90 * t113;
t173 = t106 * t116;
t170 = t39 * qJD(4);
t41 = t175 + t177;
t169 = t41 * qJD(2);
t42 = -t70 + t174;
t168 = t42 * qJD(2);
t136 = t165 * t117;
t82 = pkin(2) * t136;
t167 = t82 * qJD(2);
t155 = t114 * t180;
t100 = t113 * t155;
t109 = t111 * qJD(5);
t166 = t109 - t100;
t164 = qJD(2) * t113;
t163 = qJD(3) * t113;
t162 = qJD(4) * qJ(5);
t160 = t116 * qJD(5);
t19 = -t123 * t84 + t90 * t192;
t158 = t19 * qJD(4) - t37 * qJD(5) + t178;
t157 = -t18 * qJD(4) - t178;
t154 = pkin(7) * t161;
t153 = t114 * t179;
t152 = pkin(7) * t110;
t150 = -t187 / 0.2e1;
t149 = t187 / 0.2e1;
t144 = t81 * t164;
t143 = t106 * t164;
t142 = qJD(2) * t173;
t105 = pkin(7) + t188;
t141 = t105 * t161;
t140 = t105 * t110;
t139 = t112 / 0.2e1 + t111 / 0.2e1;
t137 = pkin(2) * t159;
t44 = t159 * t85;
t43 = t159 * t84;
t135 = t113 * t159;
t134 = t116 * t155;
t131 = t139 * t105;
t130 = t139 * t117;
t2 = (t193 - t88 / 0.2e1 + pkin(2) * t130) * t85 + (t188 / 0.2e1 - t131 + t139 * pkin(7)) * t84;
t24 = (t105 * t136 + t114 * t81) * pkin(2);
t128 = t2 * qJD(1) + t24 * qJD(2);
t95 = t113 * t150;
t20 = t95 - t200;
t47 = -t83 + t174;
t127 = t20 * qJD(2) + t47 * qJD(3);
t96 = t116 * t149;
t21 = t147 * t116 + t175 + t96;
t46 = t175 + t176;
t126 = t21 * qJD(2) + t46 * qJD(3);
t125 = qJD(4) * t90 - qJD(5) * t113;
t102 = t113 * t160;
t124 = t102 - t134;
t53 = t95 + t198;
t122 = pkin(3) * t163 + t53 * qJD(2);
t97 = t116 * t150;
t54 = t148 * t116 + t97;
t121 = t54 * qJD(2) + qJD(3) * t190;
t94 = t113 * t149;
t25 = t94 + t182;
t119 = t25 * qJD(2) + t88 * t163;
t75 = qJD(4) * t129 + t160;
t103 = t113 * t110;
t99 = t116 * t153;
t98 = t113 * t153;
t92 = t101 * qJD(4);
t89 = t159 * t111;
t77 = t116 * t135;
t76 = t82 * qJD(3);
t56 = -t190 / 0.2e1 + t173 / 0.2e1 + t97;
t55 = t95 - t198;
t26 = t94 - t182;
t23 = t95 + t200;
t22 = -t175 - t177 / 0.2e1 - t176 / 0.2e1 + t96;
t13 = t183 / 0.2e1 + t186 / 0.2e1 + t118;
t11 = t165 * t43;
t10 = t37 * qJD(4) - t116 * t44;
t9 = -t113 * t44 - t170;
t8 = t135 * t85 + t170;
t1 = t185 / 0.2e1 + t184 / 0.2e1 - t84 * t131 + (t114 * t192 + t85 * t130) * pkin(2) - t199 / 0.2e1;
t4 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t159 * t6; 0, 0, -t115 * qJD(2), -t191 * qJD(2), 0, -t44, t43, 0, 0, 0, 0, 0, t10, t8, t10, -t11, t9, (-t105 * t138 + t185) * qJD(2) + t1 * qJD(3) + t158; 0, 0, 0, 0, 0, -t44, t43, 0, 0, 0, 0, 0, t10, t8, t10, -t11, t9, t1 * qJD(2) + (t184 - t199) * qJD(3) + t158; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t194, t3, 0, -t194, t159 * t19 + t75 * t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2 * qJD(3) + t157; 0, 0, 0, 0, 0, -t155, -t117 * t180, t103, t92, 0, 0, 0, t106 * t161 - t134, t106 * t110 + t100, -t42 * qJD(4) + t124, t76, -t41 * qJD(4) + t166, t24 * qJD(3) + t125 * t81; 0, 0, 0, 0, 0, -t114 * t137, -t117 * t137, t103, t92, 0, 0, 0, t55 * qJD(4) - t134 - t99, t56 * qJD(4) + t100 + t98, t23 * qJD(4) + t124 - t99, t76 + t167, t22 * qJD(4) + t166 - t98, t13 * qJD(4) + t26 * qJD(5) + (pkin(7) * t136 + t114 * t88) * t180 + t128; 0, 0, 0, 0, 0, 0, 0, t77, t197, t110, -t161, 0, t55 * qJD(3) - t140 + t143, t56 * qJD(3) + t141 + t142, t23 * qJD(3) - t140 - t168, t75, t22 * qJD(3) - t141 - t169, t13 * qJD(3) + t105 * t75 - t195; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77, t110, t89, t26 * qJD(3) + t140 - t144; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2 * qJD(2) + t157; 0, 0, 0, 0, 0, t153, t117 * t179, t103, t92, 0, 0, 0, -t53 * qJD(4) + t99, -t54 * qJD(4) - t98, -t20 * qJD(4) + t102 + t99, -t167, -t21 * qJD(4) + t109 + t98, -t12 * qJD(4) - t25 * qJD(5) - t128; 0, 0, 0, 0, 0, 0, 0, t103, t92, 0, 0, 0, -pkin(3) * t161, -pkin(3) * t110, -t47 * qJD(4) + t102, 0, -t46 * qJD(4) + t109, t125 * t88; 0, 0, 0, 0, 0, 0, 0, t77, t197, t110, -t161, 0, -t122 - t152, -t121 + t154, -t127 - t152, t75, -t126 - t154, pkin(7) * t75 - t196; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77, t110, t89, -t119 + t152; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t159 * t18; 0, 0, 0, 0, 0, 0, 0, -t77, -t197, 0, 0, 0, t53 * qJD(3) - t143, t54 * qJD(3) - t142, t20 * qJD(3) + t168, 0, t21 * qJD(3) + t169, t12 * qJD(3) + t195; 0, 0, 0, 0, 0, 0, 0, -t77, -t197, 0, 0, 0, t122, t121, t127, 0, t126, t196; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(5), qJ(5) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4), t162; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t77, 0, -t89, t25 * qJD(3) + t144; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t77, 0, -t89, t119; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(4), -t162; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t4;
