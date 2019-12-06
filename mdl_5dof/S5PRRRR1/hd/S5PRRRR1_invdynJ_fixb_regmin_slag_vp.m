% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5PRRRR1
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% 
% Output:
% tau_reg [5x25]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:03
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRRRR1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR1_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR1_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5PRRRR1_invdynJ_fixb_regmin_slag_vp: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:03:26
% EndTime: 2019-12-05 17:03:32
% DurationCPUTime: 1.74s
% Computational Cost: add. (1259->220), mult. (3180->333), div. (0->0), fcn. (2642->10), ass. (0->132)
t92 = cos(qJ(2));
t168 = g(1) * t92;
t88 = sin(qJ(2));
t173 = g(3) * t88 + t168;
t84 = qJ(3) + qJ(4);
t78 = sin(t84);
t79 = cos(t84);
t172 = g(2) * t79 + t173 * t78;
t151 = qJD(1) * t88;
t87 = sin(qJ(3));
t131 = t87 * t151;
t123 = qJD(3) * t131;
t91 = cos(qJ(3));
t133 = t91 * t151;
t139 = t88 * qJDD(1);
t150 = qJD(2) * t87;
t174 = qJD(4) * t133 + t87 * t139 - qJDD(3) * pkin(2) - (-qJD(3) * t88 * t91 - t92 * t150) * qJD(1);
t86 = sin(qJ(4));
t90 = cos(qJ(4));
t105 = -t86 * t123 + t174 * t90;
t136 = qJD(1) * qJD(2);
t104 = (t92 * t136 + t139) * t91;
t56 = qJD(3) * pkin(2) - t131;
t99 = qJD(4) * t56 + t104;
t178 = t99 * t86;
t7 = t105 + t178;
t179 = -t7 + t172;
t177 = t99 * t90;
t138 = t91 * qJDD(2);
t140 = t87 * qJDD(2);
t50 = t86 * t91 + t90 * t87;
t81 = qJD(3) + qJD(4);
t175 = t81 * t50;
t17 = qJD(2) * t175 - t90 * t138 + t86 * t140;
t147 = qJD(3) * t87;
t15 = qJDD(5) + t17;
t132 = t86 * t150;
t148 = qJD(2) * t91;
t72 = t90 * t148;
t45 = -t72 + t132;
t44 = qJD(5) + t45;
t171 = (t44 * t147 - t15 * t91) * pkin(2) - t44 * t151;
t121 = -t90 * t123 - t174 * t86;
t57 = -pkin(2) * t148 - t92 * qJD(1);
t126 = qJD(5) * t57 + t121 + t177;
t49 = t86 * t87 - t90 * t91;
t22 = t81 * t49;
t31 = t86 * t133 - t90 * t56;
t109 = t49 * t92;
t37 = qJD(1) * t109;
t170 = -t126 * t49 - t31 * t22 + t7 * t50 + (qJD(5) * t91 * pkin(2) - t37) * t44;
t142 = qJDD(1) - g(3);
t169 = g(1) * t88;
t108 = t142 * t92 + t169;
t47 = -t86 * t148 - t90 * t150;
t85 = sin(qJ(5));
t89 = cos(qJ(5));
t115 = t89 * t47 - t85 * t81;
t135 = qJD(2) * qJD(3);
t128 = t91 * t135;
t16 = qJD(4) * t72 - t81 * t132 + t86 * t138 + (t128 + t140) * t90;
t80 = qJDD(3) + qJDD(4);
t9 = -t115 * qJD(5) + t85 * t16 - t89 * t80;
t167 = g(2) * t78;
t143 = qJD(5) * t89;
t144 = qJD(5) * t85;
t8 = t81 * t143 + t47 * t144 + t89 * t16 + t85 * t80;
t165 = t8 * t85;
t28 = -t85 * t47 - t89 * t81;
t164 = t28 * t44;
t163 = t115 * t44;
t162 = t31 * t50;
t161 = t44 * t47;
t160 = t47 * t45;
t159 = t85 * t15;
t158 = t88 * t85;
t157 = t88 * t89;
t156 = t89 * t15;
t155 = t92 * t85;
t154 = t92 * t89;
t82 = t87 ^ 2;
t153 = -t91 ^ 2 + t82;
t93 = qJD(3) ^ 2;
t94 = qJD(2) ^ 2;
t152 = t93 + t94;
t149 = qJD(2) * t88;
t145 = qJD(4) * t90;
t141 = qJDD(3) * t87;
t137 = t92 * qJDD(2);
t130 = t31 * (-t44 + t45);
t129 = t87 * t135;
t125 = t44 * t89;
t58 = t90 * t133;
t32 = t86 * t56 + t58;
t33 = t88 * t136 - t92 * qJDD(1) + (t129 - t138) * pkin(2);
t124 = qJD(5) * t32 - t33;
t120 = -g(3) * t92 + t169;
t36 = t50 * t151;
t118 = t31 * t45 - t36 * t44;
t39 = t49 * t88;
t117 = -t89 * t39 - t155;
t116 = t85 * t39 - t154;
t114 = qJD(5) * t86 + t150;
t113 = -t126 - t167;
t111 = -t50 * t144 - t89 * t22;
t110 = t50 * t92;
t107 = -t142 * t88 + t168;
t19 = -t85 * t32 + t89 * t57;
t103 = t31 * t144 + t179 * t89 + t19 * t47;
t102 = t173 * t79 + t57 * t45 - t121 - t167;
t98 = t57 * t47 - t105 + t172;
t96 = -t104 + (-pkin(2) * t81 - t56) * qJD(4);
t20 = t89 * t32 + t85 * t57;
t95 = t31 * t143 - t179 * t85 - t20 * t47;
t43 = t79 * t154 + t158;
t42 = -t79 * t155 + t157;
t41 = -t79 * t157 + t155;
t40 = t79 * t158 + t154;
t38 = t50 * t88;
t35 = qJD(1) * t110;
t34 = t86 * t131 - t58;
t27 = t89 * t33;
t18 = -t45 ^ 2 + t47 ^ 2;
t13 = qJD(2) * t110 - t22 * t88;
t12 = -qJD(2) * t109 - t175 * t88;
t11 = -t47 * t81 - t17;
t10 = t45 * t81 + t16;
t4 = -t115 * t47 + t44 * t125 + t159;
t3 = -t44 ^ 2 * t85 - t28 * t47 + t156;
t2 = -t115 * t125 + t165;
t1 = (t8 - t164) * t89 + (-t9 + t163) * t85;
t5 = [t142, 0, -t94 * t88 + t137, -qJDD(2) * t88 - t94 * t92, 0, 0, 0, 0, 0, (-0.2e1 * t129 + t138) * t92 + (-t152 * t91 - t141) * t88, (-qJDD(3) * t88 - 0.2e1 * t135 * t92) * t91 + (t152 * t88 - t137) * t87, 0, 0, 0, 0, 0, -t13 * t81 + t149 * t45 - t92 * t17 - t38 * t80, -t12 * t81 - t149 * t47 - t92 * t16 + t39 * t80, 0, 0, 0, 0, 0, (-qJD(5) * t117 - t85 * t12 + t149 * t89) * t44 + t116 * t15 + t13 * t28 + t38 * t9, -(qJD(5) * t116 + t89 * t12 + t149 * t85) * t44 - t117 * t15 - t13 * t115 + t38 * t8; 0, qJDD(2), t108, t107, t82 * qJDD(2) + 0.2e1 * t87 * t128, -0.2e1 * t153 * t135 + 0.2e1 * t87 * t138, t93 * t91 + t141, qJDD(3) * t91 - t93 * t87, 0, t108 * t91, -t108 * t87, t16 * t50 + t47 * t22, -t16 * t49 - t50 * t17 + t175 * t47 + t22 * t45, -t22 * t81 + t50 * t80, -t175 * t81 - t49 * t80, 0, -t45 * t151 + t57 * t175 + t33 * t49 + t35 * t81 + t120 * t79 + (t147 * t45 - t17 * t91) * pkin(2), t47 * t151 - t57 * t22 + t33 * t50 - t37 * t81 - t120 * t78 + (-t147 * t47 - t16 * t91) * pkin(2), t8 * t89 * t50 - t111 * t115, -(t115 * t85 - t28 * t89) * t22 + (-t165 - t89 * t9 + (t115 * t89 + t28 * t85) * qJD(5)) * t50, t111 * t44 - t115 * t175 + t156 * t50 + t8 * t49, -t50 * t159 - t28 * t175 - t9 * t49 + (-t143 * t50 + t85 * t22) * t44, t15 * t49 + t175 * t44, -g(1) * t41 - g(3) * t43 + t19 * t175 + t27 * t49 - t35 * t28 + ((-t32 * t49 + t162) * qJD(5) + t171) * t89 + t170 * t85, -g(1) * t40 - g(3) * t42 - t20 * t175 + t35 * t115 + t170 * t89 + (-qJD(5) * t162 + t124 * t49 - t171) * t85; 0, 0, 0, 0, -t87 * t94 * t91, t153 * t94, t140, t138, qJDD(3), g(2) * t91 + t107 * t87, -g(2) * t87 + t107 * t91, -t160, t18, t10, t11, t80, -t34 * t81 + (-t150 * t45 + t80 * t90) * pkin(2) + t96 * t86 + t98, -t36 * t81 + (t150 * t47 - t80 * t86) * pkin(2) + t96 * t90 + t102, t2, t1, t4, t3, t161, t34 * t28 + t118 * t85 + (-t90 * t9 + (qJD(4) * t28 - t159) * t86 + (-t114 * t89 - t145 * t85) * t44) * pkin(2) + t103, -t34 * t115 + t118 * t89 + (-t90 * t8 + (-qJD(4) * t115 - t156) * t86 + (t114 * t85 - t145 * t89) * t44) * pkin(2) + t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t160, t18, t10, t11, t80, t32 * t81 - t178 + t98, -t31 * t81 + t102 - t177, t2, t1, t4, t3, t161, t130 * t85 - t32 * t28 + t103, t115 * t32 + t130 * t89 + t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t115 * t28, t115 ^ 2 - t28 ^ 2, t8 + t164, -t9 - t163, t15, -g(1) * t42 + g(3) * t40 + t113 * t85 + t115 * t31 - t143 * t32 + t20 * t44 + t27, g(1) * t43 - g(3) * t41 + t113 * t89 + t124 * t85 + t19 * t44 + t31 * t28;];
tau_reg = t5;
