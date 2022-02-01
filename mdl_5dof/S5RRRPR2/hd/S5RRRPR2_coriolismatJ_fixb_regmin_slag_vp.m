% Calculate minimal parameter regressor of coriolis matrix for
% S5RRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x17]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:31
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RRRPR2_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR2_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR2_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR2_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:30:53
% EndTime: 2022-01-20 11:30:56
% DurationCPUTime: 0.94s
% Computational Cost: add. (1364->166), mult. (2981->246), div. (0->0), fcn. (2520->8), ass. (0->144)
t128 = -qJD(2) - qJD(3);
t121 = qJD(1) - t128;
t100 = sin(qJ(5));
t103 = cos(qJ(5));
t90 = -t100 ^ 2 + t103 ^ 2;
t173 = t121 * t90;
t101 = sin(qJ(3));
t105 = cos(qJ(2));
t159 = t105 * pkin(1);
t122 = pkin(2) + t159;
t114 = t101 * t122;
t102 = sin(qJ(2));
t104 = cos(qJ(3));
t145 = t104 * t102;
t74 = pkin(1) * t145 + t114;
t98 = sin(pkin(9));
t158 = t98 * t74;
t147 = t101 * t102;
t87 = t104 * t122;
t73 = pkin(1) * t147 - t87;
t70 = pkin(3) - t73;
t99 = cos(pkin(9));
t37 = t99 * t70 - t158;
t33 = -pkin(4) - t37;
t172 = -t33 / 0.2e1;
t58 = t99 * t74;
t44 = -t98 * t73 + t58;
t171 = t44 / 0.2e1;
t45 = -t99 * t73 - t158;
t170 = -t45 / 0.2e1;
t146 = t101 * t105;
t79 = (t145 + t146) * pkin(1);
t144 = t104 * t105;
t80 = (t144 - t147) * pkin(1);
t46 = t99 * t79 + t98 * t80;
t169 = -t46 / 0.2e1;
t152 = t101 * t98;
t160 = t104 * pkin(2);
t94 = pkin(3) + t160;
t71 = -pkin(2) * t152 + t99 * t94;
t67 = -pkin(4) - t71;
t168 = -t67 / 0.2e1;
t151 = t101 * t99;
t77 = (t104 * t98 + t151) * pkin(2);
t167 = t77 / 0.2e1;
t78 = (t104 * t99 - t152) * pkin(2);
t166 = -t78 / 0.2e1;
t93 = -t99 * pkin(3) - pkin(4);
t165 = -t93 / 0.2e1;
t164 = -t100 / 0.2e1;
t163 = t100 / 0.2e1;
t162 = -t103 / 0.2e1;
t161 = t103 / 0.2e1;
t52 = t101 * pkin(2);
t38 = t98 * t70 + t58;
t157 = pkin(1) * qJD(1);
t156 = pkin(1) * qJD(2);
t155 = pkin(2) * qJD(2);
t154 = pkin(2) * qJD(3);
t153 = pkin(3) * qJD(3);
t3 = -t37 * t44 + t38 * t45;
t150 = t3 * qJD(1);
t47 = -t98 * t79 + t99 * t80;
t4 = -t37 * t46 + t38 * t47;
t149 = t4 * qJD(1);
t148 = t46 * t103;
t127 = t167 + t171;
t116 = t169 + t127;
t11 = t116 * t103;
t143 = t11 * qJD(1);
t142 = t52 * qJD(1);
t53 = t87 / 0.2e1 + (-t159 / 0.2e1 + pkin(2) / 0.2e1) * t104;
t141 = t53 * qJD(1);
t140 = t73 * qJD(1);
t139 = t74 * qJD(1);
t138 = t79 * qJD(1);
t137 = t80 * qJD(1);
t136 = qJD(1) * t100;
t135 = qJD(1) * t103;
t134 = qJD(2) * t100;
t133 = qJD(2) * t103;
t132 = qJD(3) * t100;
t131 = qJD(3) * t103;
t130 = t100 * qJD(5);
t97 = t103 * qJD(5);
t129 = -qJD(1) - qJD(2);
t126 = t33 * t136;
t125 = t33 * t135;
t124 = t44 * t136;
t123 = t46 * t136;
t120 = pkin(1) * t129;
t119 = pkin(2) * t128;
t72 = pkin(2) * t151 + t98 * t94;
t118 = t170 + t165 + t172;
t117 = -t47 / 0.2e1 + t168 + t172;
t115 = t166 + t165 + t168;
t106 = t38 * t166 + t37 * t167 + t72 * t170 + t71 * t171;
t109 = (t47 * t98 / 0.2e1 + t99 * t169) * pkin(3);
t1 = t109 + t106;
t17 = -t71 * t77 + t72 * t78;
t113 = t1 * qJD(1) - t17 * qJD(2);
t5 = t117 * t100;
t112 = t5 * qJD(1) - t67 * t134;
t6 = t117 * t103;
t111 = t6 * qJD(1) - t67 * t133;
t10 = t116 * t100;
t110 = -t10 * qJD(1) - t77 * t134;
t13 = t118 * t100;
t22 = t115 * t100;
t108 = t13 * qJD(1) + t22 * qJD(2) - t93 * t132;
t14 = t118 * t103;
t23 = t115 * t103;
t107 = t14 * qJD(1) + t23 * qJD(2) - t93 * t131;
t92 = t98 * pkin(3) + pkin(8);
t91 = t100 * t97;
t86 = t90 * qJD(5);
t82 = t93 * t161;
t81 = t93 * t163;
t76 = t80 * qJD(2);
t75 = t79 * qJD(2);
t68 = pkin(8) + t72;
t66 = t74 * qJD(3);
t65 = t73 * qJD(3);
t64 = t77 * t132;
t56 = t121 * t103 * t100;
t55 = t67 * t161;
t54 = t67 * t163;
t49 = -t160 / 0.2e1 - t87 / 0.2e1 + (t147 - t144 / 0.2e1) * pkin(1);
t48 = -t52 / 0.2e1 - t114 / 0.2e1 + (-t145 - t146 / 0.2e1) * pkin(1);
t43 = t46 * t134;
t34 = pkin(8) + t38;
t32 = t44 * t132;
t27 = t33 * t161;
t26 = t33 * t163;
t25 = t78 * t162 + t55 + t82;
t24 = t78 * t164 + t54 + t81;
t16 = t45 * t162 + t27 + t82;
t15 = t45 * t164 + t26 + t81;
t12 = -t148 / 0.2e1 - t127 * t103;
t9 = t127 * t100 + t46 * t163;
t8 = t47 * t162 + t27 + t55;
t7 = t47 * t164 + t26 + t54;
t2 = t109 - t106;
t18 = [0, 0, 0, 0, -t102 * t156, -t105 * t156, 0, -t75 - t66, -t76 + t65, t4 * qJD(2) + t3 * qJD(3), t91, t86, 0, 0, 0, t33 * t130 + (-qJD(2) * t46 - qJD(3) * t44) * t103, t33 * t97 + t32 + t43; 0, 0, 0, 0, t102 * t120, t105 * t120, 0, t48 * qJD(3) - t138 - t75, t49 * qJD(3) - t137 - t76, t149 + (-t46 * t71 + t47 * t72) * qJD(2) + t2 * qJD(3), t91, t86, 0, 0, 0, t12 * qJD(3) + t7 * qJD(5) + t129 * t148, t9 * qJD(3) + t8 * qJD(5) + t123 + t43; 0, 0, 0, 0, 0, 0, 0, t48 * qJD(2) - t139 - t66, t49 * qJD(2) + t140 + t65, t150 + t2 * qJD(2) + (-t44 * t99 + t45 * t98) * t153, t91, t86, 0, 0, 0, t12 * qJD(2) + t15 * qJD(5) + (-qJD(1) - qJD(3)) * t44 * t103, t9 * qJD(2) + t16 * qJD(5) + t124 + t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, t173, t97, -t130, 0, t7 * qJD(2) + t15 * qJD(3) - t34 * t97 + t126, t8 * qJD(2) + t16 * qJD(3) + t34 * t130 + t125; 0, 0, 0, 0, t102 * t157, t105 * t157, 0, -t52 * qJD(3) + t138, -t53 * qJD(3) + t137, -t1 * qJD(3) - t149, t91, t86, 0, 0, 0, -t11 * qJD(3) - t5 * qJD(5) + t46 * t135, t10 * qJD(3) - t6 * qJD(5) - t123; 0, 0, 0, 0, 0, 0, 0, -t101 * t154, -t104 * t154, t17 * qJD(3), t91, t86, 0, 0, 0, t67 * t130 - t77 * t131, t67 * t97 + t64; 0, 0, 0, 0, 0, 0, 0, t101 * t119 - t142, t104 * t119 - t141, (-t77 * t99 + t78 * t98) * t153 - t113, t91, t86, 0, 0, 0, t128 * t77 * t103 + t24 * qJD(5) - t143, t25 * qJD(5) - t110 + t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, t173, t97, -t130, 0, t24 * qJD(3) - t68 * t97 - t112, t25 * qJD(3) + t68 * t130 - t111; 0, 0, 0, 0, 0, 0, 0, t52 * qJD(2) + t139, t53 * qJD(2) - t140, t1 * qJD(2) - t150, t91, t86, 0, 0, 0, t11 * qJD(2) - t13 * qJD(5) + t44 * t135, -t10 * qJD(2) - t14 * qJD(5) - t124; 0, 0, 0, 0, 0, 0, 0, t101 * t155 + t142, t104 * t155 + t141, t113, t91, t86, 0, 0, 0, -t22 * qJD(5) + t77 * t133 + t143, -t23 * qJD(5) + t110; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t91, t86, 0, 0, 0, t93 * t130, t93 * t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, t173, t97, -t130, 0, -t92 * t97 - t108, t92 * t130 - t107; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t130, -t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t56, -t173, 0, 0, 0, t5 * qJD(2) + t13 * qJD(3) - t126, t6 * qJD(2) + t14 * qJD(3) - t125; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t56, -t173, 0, 0, 0, t22 * qJD(3) + t112, t23 * qJD(3) + t111; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t56, -t173, 0, 0, 0, t108, t107; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t18;
