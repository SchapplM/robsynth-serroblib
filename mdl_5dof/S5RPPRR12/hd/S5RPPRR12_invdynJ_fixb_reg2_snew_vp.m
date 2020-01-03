% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RPPRR12
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RPPRR12_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR12_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR12_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR12_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR12_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR12_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:07:22
% EndTime: 2019-12-31 18:07:28
% DurationCPUTime: 2.16s
% Computational Cost: add. (5174->255), mult. (11789->349), div. (0->0), fcn. (8090->8), ass. (0->164)
t131 = sin(qJ(4));
t127 = sin(pkin(8));
t134 = cos(qJ(4));
t128 = cos(pkin(8));
t167 = t128 * t131;
t150 = t127 * t134 + t167;
t109 = t150 * qJD(1);
t168 = t127 * t131;
t111 = (t128 * t134 - t168) * qJD(1);
t169 = t111 * t109;
t195 = qJDD(4) - t169;
t197 = t131 * t195;
t196 = t134 * t195;
t130 = sin(qJ(5));
t133 = cos(qJ(5));
t94 = -qJD(4) * t133 + t111 * t130;
t96 = qJD(4) * t130 + t111 * t133;
t73 = t96 * t94;
t164 = t111 * qJD(4);
t187 = t150 * qJDD(1);
t88 = -t187 - t164;
t79 = qJDD(5) - t88;
t189 = -t73 + t79;
t194 = t130 * t189;
t193 = t133 * t189;
t137 = qJD(1) ^ 2;
t132 = sin(qJ(1));
t135 = cos(qJ(1));
t151 = t132 * g(1) - t135 * g(2);
t147 = qJDD(2) - t151;
t144 = -t137 * qJ(2) + t147;
t156 = -0.2e1 * qJD(3) * qJD(1);
t183 = pkin(1) + qJ(3);
t192 = -qJDD(1) * t183 + t144 + t156;
t160 = t127 * qJDD(1);
t123 = t127 ^ 2;
t124 = t128 ^ 2;
t166 = t123 + t124;
t191 = pkin(3) * t160 - (pkin(6) * t166 + t183) * t137;
t190 = -pkin(6) - t183;
t188 = t166 * t137;
t104 = qJD(5) + t109;
t159 = t128 * qJDD(1);
t108 = -t131 * t160 + t134 * t159;
t165 = t109 * qJD(4);
t90 = t108 - t165;
t154 = -qJDD(4) * t133 + t130 * t90;
t44 = (qJD(5) - t104) * t96 + t154;
t92 = t94 ^ 2;
t93 = t96 ^ 2;
t103 = t104 ^ 2;
t106 = t109 ^ 2;
t107 = t111 ^ 2;
t184 = t127 * g(3);
t141 = (t156 + (-pkin(3) * t127 - qJ(2)) * t137 + t190 * qJDD(1) + t147) * t128;
t140 = t141 + t184;
t81 = -t128 * g(3) + t127 * t192;
t74 = -pkin(3) * t123 * t137 - pkin(6) * t160 + t81;
t51 = t131 * t74 - t134 * t140;
t175 = t134 * t74;
t52 = g(3) * t168 + t131 * t141 + t175;
t27 = t131 * t52 - t134 * t51;
t182 = t128 * t27;
t136 = qJD(4) ^ 2;
t82 = pkin(4) * t109 - pkin(7) * t111;
t29 = -qJDD(4) * pkin(4) - pkin(7) * t136 + t111 * t82 + t51;
t181 = t130 * t29;
t55 = t73 + t79;
t180 = t130 * t55;
t125 = qJDD(1) * qJ(2);
t152 = t135 * g(1) + t132 * g(2);
t148 = -t125 + t152;
t146 = -qJDD(3) + t148;
t161 = qJD(2) * qJD(1);
t143 = t146 - 0.2e1 * t161;
t78 = t143 - t191;
t179 = t131 * t78;
t85 = qJDD(4) + t169;
t178 = t131 * t85;
t177 = t133 * t29;
t176 = t133 * t55;
t174 = t134 * t78;
t173 = t134 * t85;
t172 = qJDD(1) * pkin(1);
t171 = t104 * t130;
t170 = t104 * t133;
t162 = qJD(5) + t104;
t158 = t131 * t73;
t157 = t134 * t73;
t155 = -pkin(4) * t134 - pkin(3);
t30 = -t136 * pkin(4) + qJDD(4) * pkin(7) - t109 * t82 + t131 * t140 + t175;
t122 = 0.2e1 * t161;
t38 = t122 + (-t90 + t165) * pkin(7) + (-t88 + t164) * pkin(4) - t146 + t191;
t16 = t130 * t30 - t133 * t38;
t17 = t130 * t38 + t133 * t30;
t6 = t130 * t16 + t133 * t17;
t28 = t131 * t51 + t134 * t52;
t98 = t137 * t183 + t143;
t153 = -t98 + t125;
t1 = t127 * (t131 * t29 + t134 * t6) + t128 * (t131 * t6 - t134 * t29);
t57 = t127 * t81 + t128 * (t128 * t192 + t184);
t5 = t130 * t17 - t133 * t16;
t149 = -t130 * qJDD(4) - t133 * t90;
t64 = -qJD(5) * t94 - t149;
t114 = t166 * qJDD(1);
t113 = t127 * t188;
t112 = t128 * t188;
t105 = -t144 + t172;
t101 = -t107 - t136;
t100 = -t107 + t136;
t99 = t106 - t136;
t89 = t108 - 0.2e1 * t165;
t87 = t187 + 0.2e1 * t164;
t83 = -t136 - t106;
t77 = t104 * t94;
t76 = -t93 + t103;
t75 = t92 - t103;
t72 = -t106 - t107;
t71 = t93 - t92;
t68 = -t93 - t103;
t67 = -t101 * t131 - t173;
t66 = t101 * t134 - t178;
t65 = -t103 - t92;
t63 = -qJD(5) * t96 - t154;
t62 = t92 + t93;
t61 = t108 * t131 - t134 * t187;
t60 = -t108 * t134 - t131 * t187;
t59 = t134 * t83 - t197;
t58 = t131 * t83 + t196;
t53 = (t130 * t96 - t133 * t94) * t104;
t49 = t162 * t94 + t149;
t48 = t64 + t77;
t47 = t64 - t77;
t45 = -t162 * t96 - t154;
t43 = t133 * t64 - t171 * t96;
t42 = -t130 * t63 + t170 * t94;
t41 = t127 * t67 + t128 * t66;
t40 = t133 * t75 - t180;
t39 = -t130 * t76 + t193;
t36 = -t130 * t68 - t176;
t35 = t133 * t68 - t180;
t34 = t127 * t61 + t128 * t60;
t33 = t133 * t65 - t194;
t32 = t130 * t65 + t193;
t31 = t127 * t59 + t128 * t58;
t26 = t130 * t48 - t133 * t44;
t25 = -t130 * t47 + t133 * t45;
t24 = -t130 * t44 - t133 * t48;
t23 = -t131 * t49 + t134 * t36;
t22 = t131 * t36 + t134 * t49;
t21 = -t131 * t45 + t134 * t33;
t20 = t131 * t33 + t134 * t45;
t19 = -t131 * t62 + t134 * t26;
t18 = t131 * t26 + t134 * t62;
t15 = -pkin(7) * t35 + t177;
t13 = -pkin(7) * t32 + t181;
t12 = t127 * t28 + t182;
t11 = -pkin(4) * t35 + t17;
t10 = -pkin(4) * t32 + t16;
t9 = t127 * t23 + t128 * t22;
t8 = t127 * t21 + t128 * t20;
t7 = t127 * t19 + t128 * t18;
t2 = -pkin(7) * t24 - t5;
t3 = [0, 0, 0, 0, 0, qJDD(1), t151, t152, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, t147 - 0.2e1 * t172, t122 + 0.2e1 * t125 - t152, pkin(1) * t105 + qJ(2) * (-t137 * pkin(1) + t122 - t148), t124 * qJDD(1), -0.2e1 * t127 * t159, 0, t123 * qJDD(1), 0, 0, t113 * t183 + t127 * t153, t112 * t183 + t128 * t153, -qJ(2) * t188 + t114 * t183 - t57, -qJ(2) * t98 - t183 * t57, t128 * (-t131 * t164 + t134 * t90) - t127 * (t131 * t90 + t134 * t164), t128 * (-t131 * t89 - t134 * t87) - t127 * (-t131 * t87 + t134 * t89), t128 * (-t100 * t131 + t196) - t127 * (t100 * t134 + t197), t128 * (-t131 * t88 + t134 * t165) - t127 * (t131 * t165 + t134 * t88), t128 * (t134 * t99 - t178) - t127 * (t131 * t99 + t173), (t128 * (-t109 * t134 + t111 * t131) - t127 * (-t109 * t131 - t111 * t134)) * qJD(4), t128 * (-pkin(6) * t58 - t179) - t127 * (-pkin(3) * t87 + pkin(6) * t59 + t174) + qJ(2) * t87 - t183 * t31, t128 * (-pkin(6) * t66 - t174) - t127 * (-pkin(3) * t89 + pkin(6) * t67 - t179) + qJ(2) * t89 - t183 * t41, t128 * (-pkin(6) * t60 - t27) - t127 * (-pkin(3) * t72 + pkin(6) * t61 + t28) + qJ(2) * t72 - t183 * t34, -pkin(6) * t182 - t127 * (pkin(3) * t78 + pkin(6) * t28) - qJ(2) * t78 - t183 * t12, t128 * (t134 * t43 + t158) - t127 * (t131 * t43 - t157), t128 * (t131 * t71 + t134 * t25) - t127 * (t131 * t25 - t134 * t71), t128 * (t131 * t48 + t134 * t39) - t127 * (t131 * t39 - t134 * t48), t128 * (t134 * t42 - t158) - t127 * (t131 * t42 + t157), t128 * (-t131 * t44 + t134 * t40) - t127 * (t131 * t40 + t134 * t44), t128 * (t131 * t79 + t134 * t53) - t127 * (t131 * t53 - t134 * t79), t128 * (-pkin(6) * t20 - t10 * t131 + t13 * t134) - t127 * (-pkin(3) * t32 + pkin(6) * t21 + t10 * t134 + t13 * t131) + qJ(2) * t32 - t183 * t8, t128 * (-pkin(6) * t22 - t11 * t131 + t134 * t15) - t127 * (-pkin(3) * t35 + pkin(6) * t23 + t11 * t134 + t131 * t15) + qJ(2) * t35 - t183 * t9, t128 * (-pkin(6) * t18 + t134 * t2) - t127 * (pkin(6) * t19 + t131 * t2) - t183 * t7 + (pkin(4) * t167 - t127 * t155 + qJ(2)) * t24, (t128 * (pkin(4) * t131 - pkin(7) * t134) - t127 * (-pkin(7) * t131 + t155) + qJ(2)) * t5 + t190 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t137, -t105, 0, 0, 0, 0, 0, 0, -t113, -t112, -t114, t57, 0, 0, 0, 0, 0, 0, t31, t41, t34, t12, 0, 0, 0, 0, 0, 0, t8, t9, t7, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t160, t159, -t188, -t98, 0, 0, 0, 0, 0, 0, t87, t89, t72, -t78, 0, 0, 0, 0, 0, 0, t32, t35, t24, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t169, t107 - t106, t108, -t169, -t187, qJDD(4), -t51, -t52, 0, 0, t130 * t64 + t170 * t96, t130 * t45 + t133 * t47, t133 * t76 + t194, t133 * t63 + t171 * t94, t130 * t75 + t176, (-t130 * t94 - t133 * t96) * t104, pkin(4) * t45 + pkin(7) * t33 - t177, pkin(4) * t49 + pkin(7) * t36 + t181, pkin(4) * t62 + pkin(7) * t26 + t6, -pkin(4) * t29 + pkin(7) * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73, t71, t48, -t73, -t44, t79, -t16, -t17, 0, 0;];
tauJ_reg = t3;
