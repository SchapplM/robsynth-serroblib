% Calculate minimal parameter regressor of coriolis matrix for
% S5RPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x22]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 12:46
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RPRRP3_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP3_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP3_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP3_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 12:46:09
% EndTime: 2021-01-15 12:46:14
% DurationCPUTime: 1.23s
% Computational Cost: add. (2034->130), mult. (3684->169), div. (0->0), fcn. (3812->6), ass. (0->116)
t146 = qJD(3) + qJD(4);
t185 = cos(qJ(4));
t145 = t185 * pkin(3);
t115 = t145 + pkin(4);
t125 = t145 / 0.2e1 - t115 / 0.2e1;
t118 = sin(qJ(4));
t119 = sin(qJ(3));
t120 = cos(qJ(3));
t108 = t118 * t120 + t185 * t119;
t112 = sin(pkin(8)) * pkin(1) + pkin(6);
t180 = pkin(7) + t112;
t103 = t180 * t119;
t104 = t180 * t120;
t193 = t185 * t103 + t118 * t104;
t50 = -t108 * qJ(5) - t193;
t195 = t146 * t50;
t106 = t118 * t119 - t185 * t120;
t166 = t106 * qJ(5);
t164 = t118 * t103;
t89 = t185 * t104;
t192 = -t89 + t164;
t48 = t192 + t166;
t194 = t146 * t193;
t191 = t146 * t106;
t70 = t146 * t108;
t105 = t106 ^ 2;
t190 = t108 ^ 2;
t188 = -t89 / 0.2e1;
t187 = t48 * pkin(4);
t184 = pkin(3) * t118;
t183 = t106 * pkin(4);
t182 = t108 * pkin(4);
t181 = t119 * pkin(3);
t178 = pkin(4) * qJD(4);
t113 = -cos(pkin(8)) * pkin(1) - pkin(2);
t110 = -pkin(3) * t120 + t113;
t77 = t110 + t183;
t65 = t77 * t108;
t11 = t48 * t106 - t50 * t108;
t173 = qJD(1) * t11;
t80 = t181 + t182;
t17 = t106 * t80 + t65;
t172 = qJD(1) * t17;
t64 = t77 * t106;
t18 = t108 * t80 - t64;
t171 = qJD(1) * t18;
t21 = -t106 * t182 - t65;
t170 = qJD(1) * t21;
t22 = -pkin(4) * t190 + t64;
t169 = qJD(1) * t22;
t67 = t188 + t89 / 0.2e1;
t168 = qJD(1) * t67;
t167 = qJD(4) * t48;
t165 = t115 * t108;
t162 = t118 * t106;
t134 = -t162 / 0.2e1;
t135 = -t165 / 0.2e1;
t92 = -t182 / 0.2e1;
t31 = t135 + t92 + (t134 - t119 / 0.2e1) * pkin(3);
t159 = t31 * qJD(1);
t42 = (-pkin(4) / 0.2e1 - t125) * t106;
t158 = t42 * qJD(1);
t57 = t106 * t181 + t108 * t110;
t157 = t57 * qJD(1);
t58 = -t106 * t110 + t108 * t181;
t156 = t58 * qJD(1);
t59 = t105 - t190;
t155 = t59 * qJD(1);
t71 = t105 + t190;
t154 = t71 * qJD(1);
t153 = qJD(1) * t120;
t152 = qJD(4) * t110;
t151 = t106 * qJD(1);
t150 = t108 * qJD(1);
t100 = t108 * qJD(5);
t111 = -t119 ^ 2 + t120 ^ 2;
t149 = t111 * qJD(1);
t148 = t119 * qJD(3);
t147 = t120 * qJD(3);
t144 = t108 * t178;
t143 = qJD(4) * t184;
t142 = pkin(4) * t150;
t140 = t110 * t151;
t139 = t110 * t150;
t138 = t113 * t119 * qJD(1);
t137 = t113 * t153;
t136 = t119 * t153;
t133 = t185 * qJD(3);
t132 = t185 * qJD(4);
t129 = pkin(3) * t132;
t9 = pkin(4) * t65;
t127 = qJD(1) * t9;
t8 = t77 * t80;
t126 = t8 * qJD(1);
t123 = t125 * t108;
t121 = t125 * t48;
t3 = -t187 / 0.2e1 - t121;
t93 = t182 / 0.2e1;
t44 = t93 + t123;
t90 = (t145 - t115) * t184;
t122 = -qJD(1) * t3 - qJD(2) * t44 - qJD(3) * t90;
t36 = 0.2e1 * t188 + t164;
t96 = t106 * qJD(5);
t76 = t106 * t150;
t63 = t67 * qJD(3);
t62 = t67 * qJD(4);
t54 = qJD(3) * t184 - t168;
t53 = pkin(3) * t133;
t47 = -t146 * t184 + t168;
t46 = (-t133 - t132) * pkin(3);
t43 = t92 + t123;
t41 = t183 / 0.2e1 - t125 * t106;
t30 = pkin(3) * t134 + t135 + t181 / 0.2e1 + t93;
t20 = t36 + t166;
t2 = t187 / 0.2e1 - t121;
t1 = [0, 0, 0, 0, t119 * t147, t111 * qJD(3), 0, 0, 0, t113 * t148, t113 * t147, -t106 * t70, t146 * t59, 0, 0, 0, qJD(3) * t57 + t108 * t152, qJD(3) * t58 - t106 * t152, qJD(3) * t17 - qJD(4) * t21, qJD(3) * t18 - qJD(4) * t22, qJD(5) * t71, qJD(3) * t8 + qJD(4) * t9 + qJD(5) * t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, t136, t149, t147, -t148, 0, -t112 * t147 + t138, t112 * t148 + t137, -t76, t155, -t191, -t70, 0, qJD(3) * t192 + t36 * qJD(4) + t157, t156 + t194, qJD(3) * t48 + qJD(4) * t20 + t172, t171 - t195, (t115 * t106 - t108 * t184) * qJD(3) + t41 * qJD(4), (t115 * t48 + t184 * t50) * qJD(3) + t2 * qJD(4) + t30 * qJD(5) + t126; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t76, t155, -t191, -t70, 0, t36 * qJD(3) + qJD(4) * t192 + t139, -t140 + t194, qJD(3) * t20 + t167 - t170, -t169 - t195, qJD(3) * t41 + t106 * t178, pkin(4) * t167 + qJD(3) * t2 + t127; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t154, qJD(3) * t30 + t173; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t148, -t147, 0, 0, 0, 0, 0, -t70, t191, -t70, t191, 0, (-pkin(3) * t162 - t165) * qJD(3) + t43 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t70, t191, -t70, t191, 0, qJD(3) * t43 - t144; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, -t136, -t149, 0, 0, 0, -t138, -t137, t76, -t155, 0, 0, 0, t62 - t157, -t156, -t100 + t62 - t172, t96 - t171, qJD(4) * t42, qJD(4) * t3 + qJD(5) * t31 - t126; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4) * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t143, -t129, -t143, -t129, 0, t90 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, t46, t47, t46, t158, -pkin(4) * t143 - t122; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t150, t151, 0, t159; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t76, -t155, 0, 0, 0, -t63 - t139, t140, -t100 - t63 + t170, t96 + t169, -qJD(3) * t42, -pkin(4) * t100 - qJD(3) * t3 - t127; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3) * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, t53, t54, t53, -t158, t122; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t150, t151, 0, -t142; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70, -t191, -t154, -qJD(3) * t31 + t144 - t173; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t150, -t151, 0, -t159; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t150, -t151, 0, t142; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t1;
