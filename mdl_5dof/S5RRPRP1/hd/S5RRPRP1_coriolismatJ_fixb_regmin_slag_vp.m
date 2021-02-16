% Calculate minimal parameter regressor of coriolis matrix for
% S5RRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x18]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 20:09
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RRPRP1_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP1_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 20:08:56
% EndTime: 2021-01-15 20:09:00
% DurationCPUTime: 1.10s
% Computational Cost: add. (1240->142), mult. (2357->186), div. (0->0), fcn. (1880->6), ass. (0->131)
t102 = sin(qJ(4));
t101 = cos(pkin(8));
t100 = sin(pkin(8));
t103 = sin(qJ(2));
t91 = t100 * t103 * pkin(1);
t105 = cos(qJ(2));
t157 = t105 * pkin(1);
t94 = pkin(2) + t157;
t119 = t101 * t94 - t91;
t66 = -pkin(3) - t119;
t93 = -t101 * pkin(2) - pkin(3);
t123 = -t93 / 0.2e1 - t66 / 0.2e1;
t170 = t102 * t123;
t104 = cos(qJ(4));
t169 = t104 * t123;
t130 = qJD(1) + qJD(2);
t98 = t102 ^ 2;
t99 = t104 ^ 2;
t88 = t98 + t99;
t166 = t130 * t88;
t89 = t99 - t98;
t165 = t130 * t89;
t71 = t101 * t157 - t91;
t164 = -t71 / 0.2e1;
t163 = t71 / 0.2e1;
t162 = t98 * pkin(4);
t159 = pkin(4) * t102;
t158 = t104 * pkin(4);
t30 = t88 * t71;
t84 = t88 * qJD(5);
t156 = t30 * qJD(2) + t84;
t135 = t102 * qJD(2);
t144 = t101 * t103;
t145 = t100 * t105;
t70 = (t144 + t145) * pkin(1);
t53 = t70 * t135;
t136 = t102 * qJD(1);
t54 = t70 * t136;
t155 = t53 + t54;
t154 = pkin(1) * qJD(1);
t153 = pkin(1) * qJD(2);
t117 = pkin(1) * t144 + t100 * t94;
t67 = pkin(7) + t117;
t147 = qJ(5) + t67;
t44 = t147 * t104;
t150 = t44 * t104;
t43 = t147 * t102;
t151 = t43 * t102;
t18 = t150 + t151;
t48 = t66 - t158;
t4 = t18 * t71 + t48 * t70;
t152 = t4 * qJD(1);
t45 = t48 * t102;
t92 = t100 * pkin(2) + pkin(7);
t146 = qJ(5) + t92;
t74 = t146 * t102;
t149 = t74 * t102;
t75 = t146 * t104;
t148 = t75 * t104;
t80 = t93 - t158;
t76 = t80 * t102;
t12 = t117 * t71 - t119 * t70;
t143 = t12 * qJD(1);
t142 = t18 * qJD(1);
t141 = t30 * qJD(1);
t129 = t102 * t158;
t33 = -t45 + t129;
t140 = t33 * qJD(1);
t39 = t48 * t104 + t162;
t139 = t39 * qJD(1);
t138 = t44 * qJD(4);
t137 = t75 * qJD(4);
t96 = t102 * qJD(4);
t134 = t102 * qJD(5);
t133 = t104 * qJD(1);
t132 = t104 * qJD(2);
t97 = t104 * qJD(4);
t131 = t104 * qJD(5);
t128 = pkin(4) * t97;
t127 = pkin(4) * t134;
t125 = t76 / 0.2e1 + t45 / 0.2e1;
t124 = t80 / 0.2e1 + t48 / 0.2e1;
t122 = t66 * t136;
t121 = t66 * t133;
t120 = t70 * t132;
t51 = t102 * t164;
t118 = pkin(1) * t130;
t81 = t130 * t102;
t82 = t130 * t104;
t116 = t163 + t124;
t55 = t70 * t133;
t115 = -t55 - t120;
t32 = t148 + t149;
t106 = (t145 / 0.2e1 + t144 / 0.2e1) * pkin(1);
t5 = (-t75 / 0.2e1 - t44 / 0.2e1) * t104 + (-t74 / 0.2e1 - t43 / 0.2e1) * t102 + t106;
t114 = -t5 * qJD(1) + t32 * qJD(2);
t36 = t150 / 0.2e1;
t16 = t36 - t150 / 0.2e1;
t3 = pkin(4) * t45;
t113 = t3 * qJD(1) + t16 * qJD(3);
t10 = (-t158 + t163) * t102 + t125;
t49 = -t76 + t129;
t112 = t10 * qJD(1) - t49 * qJD(2);
t13 = t104 * t116 + t162;
t65 = t80 * t104 + t162;
t111 = t13 * qJD(1) + t65 * qJD(2);
t62 = t148 / 0.2e1;
t29 = t62 - t148 / 0.2e1;
t110 = -t16 * qJD(1) - t29 * qJD(2);
t19 = t51 + t170;
t109 = t19 * qJD(1) - t135 * t93;
t52 = t104 * t164;
t20 = t52 + t169;
t108 = t20 * qJD(1) - t132 * t93;
t1 = t116 * t159;
t9 = pkin(4) * t76;
t107 = t1 * qJD(1) + t9 * qJD(2) + t29 * qJD(3);
t95 = pkin(4) * t96;
t90 = t102 * t97;
t85 = t89 * qJD(4);
t77 = pkin(4) * t81;
t69 = t102 * t82;
t27 = t29 * qJD(4);
t22 = t52 - t169;
t21 = t51 - t170;
t15 = t16 * qJD(4);
t14 = t104 * t124 + t162 + t52;
t11 = t125 + t51 - t129;
t6 = t62 + t36 + t149 / 0.2e1 + t151 / 0.2e1 + t106;
t2 = pkin(4) * t51 + (t48 + t80) * t159 / 0.2e1;
t7 = [0, 0, 0, 0, -t103 * t153, -t105 * t153, t12 * qJD(2), t90, t85, 0, 0, 0, t66 * t96 - t120, t66 * t97 + t53, -t33 * qJD(4) - t120, t39 * qJD(4) + t53, t156, t4 * qJD(2) + t3 * qJD(4) + t18 * qJD(5); 0, 0, 0, 0, -t103 * t118, -t105 * t118, t143 + (t100 * t71 - t101 * t70) * qJD(2) * pkin(2), t90, t85, 0, 0, 0, t21 * qJD(4) + t115, t22 * qJD(4) + t155, t11 * qJD(4) + t115, t14 * qJD(4) + t155, t141 + t156, t152 + (t32 * t71 + t70 * t80) * qJD(2) + t2 * qJD(4) + t6 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15; 0, 0, 0, 0, 0, 0, 0, t69, t165, t97, -t96, 0, t21 * qJD(2) - t67 * t97 + t122, t22 * qJD(2) + t67 * t96 + t121, t11 * qJD(2) - t138 - t140, t14 * qJD(2) + t43 * qJD(4) + t139, -t128, -pkin(4) * t138 + t2 * qJD(2) + t113; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t166, t6 * qJD(2) + t142; 0, 0, 0, 0, t103 * t154, t105 * t154, -t143, t90, t85, 0, 0, 0, -t19 * qJD(4) + t55, -t20 * qJD(4) - t54, t10 * qJD(4) + t55, t13 * qJD(4) - t54, t84 - t141, t1 * qJD(4) - t5 * qJD(5) - t152; 0, 0, 0, 0, 0, 0, 0, t90, t85, 0, 0, 0, t93 * t96, t93 * t97, -t49 * qJD(4), t65 * qJD(4), t84, t9 * qJD(4) + t32 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27; 0, 0, 0, 0, 0, 0, 0, t69, t165, t97, -t96, 0, -t92 * t97 - t109, t92 * t96 - t108, t112 - t137, t74 * qJD(4) + t111, -t128, -pkin(4) * t137 + t107; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t166, t114; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t96, -t97, -t96, -t97, 0, -t110 - t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, -t69, -t165, 0, 0, 0, t19 * qJD(2) - t122, t20 * qJD(2) - t121, -t10 * qJD(2) - t134 + t140, -t13 * qJD(2) - t131 - t139, 0, -t1 * qJD(2) - t113 - t127; 0, 0, 0, 0, 0, 0, 0, -t69, -t165, 0, 0, 0, t109, t108, -t112 - t134, -t111 - t131, 0, -t107 - t127; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t110; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t81, -t82, 0, -t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t96, t97, -t166, t5 * qJD(2) - t142 + t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t96, t97, -t166, -t114 + t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81, t82, 0, t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t7;
