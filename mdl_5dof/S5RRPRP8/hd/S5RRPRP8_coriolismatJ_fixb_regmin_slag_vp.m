% Calculate minimal parameter regressor of coriolis matrix for
% S5RRPRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x23]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RRPRP8_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP8_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP8_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP8_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:04:29
% EndTime: 2019-12-31 20:04:33
% DurationCPUTime: 1.32s
% Computational Cost: add. (1673->157), mult. (2952->223), div. (0->0), fcn. (2872->4), ass. (0->136)
t141 = qJD(2) - qJD(4);
t114 = cos(qJ(4));
t113 = sin(qJ(2));
t203 = pkin(6) - pkin(7);
t94 = t203 * t113;
t177 = t114 * t94;
t112 = sin(qJ(4));
t115 = cos(qJ(2));
t95 = t203 * t115;
t84 = t112 * t95;
t200 = t84 - t177;
t82 = -t115 * t112 + t113 * t114;
t123 = -qJ(5) * t82 - t200;
t186 = t112 / 0.2e1;
t208 = t123 * t186;
t198 = t82 ^ 2;
t80 = t113 * t112 + t114 * t115;
t79 = t80 ^ 2;
t39 = t79 - t198;
t207 = t39 * qJD(1);
t142 = t115 * qJD(3);
t153 = -t115 * pkin(2) - t113 * qJ(3);
t204 = qJD(2) * t153 + t142;
t193 = t177 / 0.2e1;
t190 = -pkin(2) - pkin(3);
t88 = qJ(3) * t112 - t114 * t190;
t202 = t88 + pkin(4);
t201 = t141 * t80;
t124 = t141 * t82;
t41 = t112 * t94 + t114 * t95;
t199 = t141 * t41;
t197 = -pkin(4) / 0.2e1;
t116 = -qJ(5) * t80 + t41;
t196 = t116 / 0.2e1;
t194 = -t116 / 0.2e1;
t192 = t202 / 0.2e1;
t191 = -t202 / 0.2e1;
t189 = t80 * pkin(4);
t188 = t82 * pkin(4);
t187 = -t112 / 0.2e1;
t185 = t113 / 0.2e1;
t184 = t114 / 0.2e1;
t90 = -pkin(1) + t153;
t72 = t115 * pkin(3) - t90;
t183 = t72 * t80;
t89 = t114 * qJ(3) + t112 * t190;
t182 = t89 * t80;
t175 = t116 * t114;
t181 = t123 * t187 + t175 / 0.2e1;
t179 = pkin(4) * qJD(4);
t48 = t72 + t189;
t106 = t115 * qJ(3);
t77 = t190 * t113 + t106;
t3 = t48 * (t77 - t188);
t176 = t3 * qJD(1);
t173 = t202 * t112;
t13 = -t116 * t80 - t123 * t82;
t172 = qJD(1) * t13;
t53 = t193 - t177 / 0.2e1;
t171 = qJD(1) * t53;
t170 = qJD(1) * t82;
t169 = qJD(4) * t80;
t168 = qJD(4) * t82;
t167 = qJD(4) * t89;
t130 = t88 / 0.2e1 + t191;
t12 = (pkin(4) / 0.2e1 + t130) * t80;
t166 = t12 * qJD(1);
t117 = (-pkin(2) / 0.2e1 - pkin(3) / 0.2e1) * t113 + t106 / 0.2e1;
t14 = t182 / 0.2e1 + (t197 + t191) * t82 + t117;
t165 = t14 * qJD(1);
t18 = -t72 * t82 + t77 * t80;
t164 = t18 * qJD(1);
t19 = t77 * t82 + t183;
t163 = t19 * qJD(1);
t118 = t82 * t184 + t80 * t186;
t42 = t185 + t118;
t158 = t42 * qJD(1);
t49 = t79 + t198;
t157 = t49 * qJD(1);
t93 = t113 * pkin(2) - t106;
t50 = t93 * t113 + t115 * t90;
t156 = t50 * qJD(1);
t51 = -t90 * t113 + t115 * t93;
t155 = t51 * qJD(1);
t111 = t113 ^ 2;
t96 = t115 ^ 2 - t111;
t154 = t96 * qJD(1);
t152 = qJD(1) * t113;
t151 = qJD(1) * t115;
t150 = qJD(2) * qJ(3);
t149 = qJD(3) * t112;
t148 = qJD(3) * t113;
t147 = qJD(3) * t114;
t146 = t111 * qJD(1);
t145 = t112 * qJD(2);
t144 = t113 * qJD(2);
t143 = t114 * qJD(2);
t104 = t115 * qJD(2);
t140 = pkin(4) * t170;
t139 = pkin(4) * t196;
t138 = pkin(1) * t152;
t137 = pkin(1) * t151;
t136 = pkin(6) * t104;
t135 = pkin(6) * t144;
t134 = qJD(1) * t183;
t133 = t72 * t170;
t132 = t80 * t170;
t131 = t90 * t93 * qJD(1);
t129 = t48 * t152;
t128 = t80 * t152;
t127 = t82 * t152;
t126 = t90 * t152;
t6 = t48 * t188;
t8 = t208 - t175 / 0.2e1 + t181;
t122 = qJD(1) * t6 + qJD(3) * t8;
t1 = t130 * t116 + t139;
t25 = (t202 - t88) * t89;
t121 = qJD(1) * t1 - qJD(2) * t25;
t16 = (t197 + t130) * t112;
t120 = -qJD(1) * t8 + qJD(2) * t16;
t47 = t89 * t114 + t173;
t9 = (t196 + t194) * t114;
t119 = qJD(1) * t9 - qJD(2) * t47;
t97 = t113 * t151;
t92 = t141 * t114;
t91 = t141 * t112;
t43 = t185 - t118;
t40 = -t84 + 0.2e1 * t193;
t24 = t112 * t82 - t114 * t80;
t17 = t173 / 0.2e1 + t202 * t187;
t15 = -t182 / 0.2e1 + t82 * t192 - t188 / 0.2e1 + t117;
t11 = -t189 / 0.2e1 + t130 * t80;
t10 = t116 * t184 + t181 - t208;
t7 = t8 * qJD(4);
t2 = t116 * t192 + t88 * t194 + t139;
t4 = [0, 0, 0, t113 * t104, t96 * qJD(2), 0, 0, 0, -pkin(1) * t144, -pkin(1) * t104, -t51 * qJD(2) + t113 * t142, 0, -qJD(2) * t50 + qJD(3) * t111, (qJD(2) * t93 - t148) * t90, t80 * t124, -t141 * t39, 0, 0, 0, qJD(2) * t18 + t148 * t80 + t168 * t72, qJD(2) * t19 + t148 * t82 - t169 * t72, qJD(5) * t49, qJD(2) * t3 + qJD(4) * t6 + qJD(5) * t13 + t148 * t48; 0, 0, 0, t97, t154, t104, -t144, 0, -t136 - t138, t135 - t137, -t136 - t155, t204, -t135 - t156, t204 * pkin(6) + t131, t132, -t207, -t201, -t124, 0, t164 - t199, t200 * qJD(2) + t40 * qJD(4) + t163, (t202 * t80 + t89 * t82) * qJD(2) + t24 * qJD(3) + t11 * qJD(4), t176 + (-t116 * t202 - t123 * t89) * qJD(2) + t10 * qJD(3) + t2 * qJD(4) + t15 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t97, t104, t146, -t126 + t136, 0, 0, 0, 0, 0, t128, t127, qJD(2) * t24, qJD(2) * t10 + qJD(5) * t43 + t129 + t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t132, t207, t201, t124, 0, t133 + t199, t40 * qJD(2) + t200 * qJD(4) - t134, pkin(4) * t169 + qJD(2) * t11, qJD(2) * t2 - t116 * t179 + t122; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t157, qJD(2) * t15 + qJD(3) * t43 + t172; 0, 0, 0, -t97, -t154, 0, 0, 0, t138, t137, t155, 0, t156, -t131, -t132, t207, 0, 0, 0, -t164, qJD(4) * t53 - t163, qJD(4) * t12, -qJD(3) * t9 - qJD(4) * t1 - qJD(5) * t14 - t176; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3), qJ(3) * qJD(3), 0, 0, 0, 0, 0, t149 + t167, -qJD(4) * t88 + t147, 0, qJD(3) * t47 + qJD(4) * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2), t150, 0, 0, 0, 0, 0, t145, t143, 0, qJD(4) * t17 - t119; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t141 * t89, -t141 * t88 + t171, t166, -pkin(4) * t167 + qJD(3) * t17 - t121; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t165; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t97, 0, -t146, t126, 0, 0, 0, 0, 0, -t128, -t127, 0, qJD(2) * t9 - qJD(5) * t42 - t129 + t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(2), -t150, 0, 0, 0, 0, 0, -t91, -t92, 0, -qJD(4) * t16 + t119; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t91, t92, 0, -t112 * t179 - t120; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t158; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t132, -t207, 0, 0, 0, -t133, -qJD(2) * t53 + t134, -qJD(2) * t12, qJD(2) * t1 - qJD(5) * t188 - t122; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(2) * t89 - t149, qJD(2) * t88 - t147 - t171, -t166, qJD(3) * t16 + t121; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t145, -t143, 0, t120; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t140; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t157, pkin(4) * t168 + qJD(2) * t14 + qJD(3) * t42 - t172; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t165; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t158; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t140; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t4;
