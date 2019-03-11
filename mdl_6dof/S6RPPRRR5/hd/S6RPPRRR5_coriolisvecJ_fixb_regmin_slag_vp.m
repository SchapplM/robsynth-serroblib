% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6]';
% 
% Output:
% tauc_reg [6x30]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPPRRR5_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR5_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR5_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR5_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:29:01
% EndTime: 2019-03-09 02:29:05
% DurationCPUTime: 1.35s
% Computational Cost: add. (1781->224), mult. (3780->310), div. (0->0), fcn. (2548->6), ass. (0->138)
t100 = cos(qJ(5));
t101 = cos(qJ(4));
t97 = sin(qJ(5));
t98 = sin(qJ(4));
t58 = t100 * t98 + t97 * t101;
t54 = t58 * qJD(1);
t177 = qJD(6) + t54;
t140 = qJD(1) * t101;
t132 = t100 * t140;
t150 = qJD(1) * t98;
t134 = t97 * t150;
t53 = -t132 + t134;
t89 = qJD(4) + qJD(5);
t96 = sin(qJ(6));
t99 = cos(qJ(6));
t43 = -t96 * t53 - t99 * t89;
t180 = t177 * t43;
t115 = t99 * t53 - t96 * t89;
t179 = t115 * t177;
t178 = qJD(6) - t177;
t159 = t54 * t89;
t95 = pkin(1) + qJ(3);
t176 = qJD(1) * t95;
t139 = qJD(4) * t101;
t144 = t98 * qJD(2);
t86 = qJD(1) * qJ(2) + qJD(3);
t65 = -qJD(1) * pkin(7) + t86;
t42 = t65 * t139 + (-pkin(8) * t139 + t144) * qJD(1);
t51 = -pkin(8) * t140 + t101 * t65;
t47 = qJD(4) * pkin(4) + t51;
t175 = (qJD(5) * t47 + t42) * t100;
t145 = t100 * t101;
t112 = t89 * t145;
t136 = qJD(1) * qJD(4);
t131 = t98 * t136;
t153 = -qJD(5) * t134 - t97 * t131;
t32 = qJD(1) * t112 + t153;
t154 = t99 * t32;
t174 = -t177 ^ 2 * t96 + t154;
t16 = -qJD(6) * t115 - t159 * t96;
t94 = -pkin(7) + qJ(2);
t171 = pkin(8) - t94;
t61 = t171 * t98;
t62 = t171 * t101;
t113 = -t100 * t62 + t97 * t61;
t149 = qJD(4) * t98;
t48 = t101 * qJD(2) + t171 * t149;
t49 = -qJD(4) * t62 + t144;
t12 = qJD(5) * t113 + t100 * t49 + t97 * t48;
t50 = -pkin(8) * t150 + t98 * t65;
t155 = t97 * t50;
t25 = t100 * t47 - t155;
t21 = -t89 * pkin(5) - t25;
t66 = -qJD(2) + t176;
t56 = pkin(4) * t150 + t66;
t24 = t54 * pkin(5) + t53 * pkin(9) + t56;
t57 = t97 * t98 - t145;
t76 = t98 * pkin(4) + t95;
t33 = t58 * pkin(5) + t57 * pkin(9) + t76;
t138 = qJD(5) * t100;
t148 = qJD(5) * t97;
t37 = -t100 * t149 - t101 * t148 - t138 * t98 - t139 * t97;
t40 = -t100 * t61 - t97 * t62;
t137 = qJD(1) * qJD(2);
t79 = t101 * t137;
t41 = t79 + (pkin(8) * qJD(1) - t65) * t149;
t129 = -t50 * t148 + t97 * t41;
t6 = t129 + t175;
t130 = -t100 * t41 + t97 * t42;
t151 = t100 * t50;
t26 = t97 * t47 + t151;
t7 = qJD(5) * t26 + t130;
t173 = -(qJD(6) * t33 + t12) * t177 - t40 * t32 - (qJD(6) * t24 + t6) * t58 - t7 * t57 + t21 * t37;
t90 = qJD(3) * qJD(1);
t172 = 0.2e1 * t90;
t146 = qJD(6) * t99;
t147 = qJD(6) * t96;
t15 = t89 * t146 + t53 * t147 - t159 * t99;
t170 = t15 * t96;
t169 = t21 * t54;
t168 = t21 * t57;
t167 = t33 * t32;
t38 = -t98 * t148 - t97 * t149 + t112;
t166 = t38 * t89;
t165 = t177 * t53;
t163 = t53 * t43;
t162 = t53 * t115;
t161 = t53 * t54;
t160 = t53 * t89;
t158 = t57 * t15;
t156 = t96 * t32;
t128 = t101 * t136;
t60 = pkin(4) * t128 + t90;
t152 = -t101 ^ 2 + t98 ^ 2;
t143 = qJD(2) - t66;
t68 = pkin(4) * t139 + qJD(3);
t102 = qJD(4) ^ 2;
t103 = qJD(1) ^ 2;
t141 = -t102 - t103;
t135 = pkin(4) * t140;
t87 = 0.2e1 * t137;
t133 = -pkin(4) * t89 - t47;
t126 = t99 * t177;
t34 = -t53 * pkin(5) + t54 * pkin(9);
t81 = t97 * pkin(4) + pkin(9);
t123 = qJD(6) * t81 + t135 + t34;
t122 = qJD(6) * t58 + qJD(1);
t22 = t89 * pkin(9) + t26;
t9 = t99 * t22 + t96 * t24;
t121 = t21 * t146 - t9 * t53 + t7 * t96;
t120 = -0.2e1 * t128;
t27 = t97 * t51 + t151;
t119 = pkin(4) * t148 - t27;
t28 = t100 * t51 - t155;
t118 = -pkin(4) * t138 + t28;
t117 = -t81 * t32 + t169;
t116 = t96 * t22 - t99 * t24;
t114 = qJD(2) + t66 + t176;
t111 = -t116 * t53 + t21 * t147 - t7 * t99;
t110 = -t102 * t94 + t172;
t109 = t56 * t53 - t130;
t108 = t56 * t54 - t129;
t107 = t57 * t147 + t99 * t37;
t105 = t126 * t177 + t156;
t104 = -t132 * t89 - t153;
t82 = -t100 * pkin(4) - pkin(5);
t35 = t37 * t89;
t23 = t53 ^ 2 - t54 ^ 2;
t20 = t104 - t160;
t14 = t38 * pkin(5) - t37 * pkin(9) + t68;
t13 = qJD(5) * t40 - t100 * t48 + t97 * t49;
t11 = t32 * pkin(5) + pkin(9) * t159 + t60;
t10 = t99 * t11;
t4 = t105 - t162;
t3 = -t163 + t174;
t2 = -t115 * t126 + t170;
t1 = (t15 - t180) * t99 + (-t16 + t179) * t96;
t5 = [0, 0, 0, 0, t87, qJ(2) * t87, t87, t172, t86 * qJD(2) + t66 * qJD(3) + (qJ(2) * qJD(2) + qJD(3) * t95) * qJD(1), t98 * t120, 0.2e1 * t152 * t136, -t102 * t98, -t102 * t101, 0, t110 * t98 + t114 * t139, t101 * t110 - t114 * t149, t159 * t57 - t53 * t37, t159 * t58 + t57 * t32 - t37 * t54 + t53 * t38, t35, -t166, 0, -t13 * t89 + t76 * t32 + t56 * t38 + t68 * t54 + t60 * t58, -t12 * t89 - t159 * t76 + t56 * t37 - t68 * t53 - t60 * t57, -t107 * t115 - t99 * t158 (t115 * t96 - t43 * t99) * t37 + (t170 + t16 * t99 + (-t115 * t99 - t43 * t96) * qJD(6)) * t57, t107 * t177 - t115 * t38 + t15 * t58 - t57 * t154, t57 * t156 - t16 * t58 - t43 * t38 + (t57 * t146 - t96 * t37) * t177, t177 * t38 + t32 * t58, t10 * t58 + t13 * t43 - t113 * t16 - t116 * t38 + (t14 * t177 + t167 + (-t177 * t40 - t22 * t58 - t168) * qJD(6)) * t99 + t173 * t96, -t13 * t115 - t113 * t15 - t9 * t38 + (-(-qJD(6) * t40 + t14) * t177 - t167 - (-qJD(6) * t22 + t11) * t58 + qJD(6) * t168) * t96 + t173 * t99; 0, 0, 0, 0, -t103, -t103 * qJ(2), -t103, 0, -t86 * qJD(1) - t90, 0, 0, 0, 0, 0, t120, 0.2e1 * t131, 0, 0, 0, 0, 0, t104 + t160, 0.2e1 * t159, 0, 0, 0, 0, 0, -t163 - t174, t105 + t162; 0, 0, 0, 0, 0, 0, 0, -t103, t143 * qJD(1), 0, 0, 0, 0, 0, t141 * t98, t141 * t101, 0, 0, 0, 0, 0, -qJD(1) * t54 + t35, qJD(1) * t53 - t166, 0, 0, 0, 0, 0, -t58 * t156 + t57 * t16 - t37 * t43 + (-t122 * t99 - t38 * t96) * t177, -t58 * t154 + t158 + t37 * t115 + (t122 * t96 - t38 * t99) * t177; 0, 0, 0, 0, 0, 0, 0, 0, 0, t101 * t103 * t98, -t152 * t103, 0, 0, 0, -t140 * t66 + t79, -t143 * t150, -t161, t23, 0, t20, 0, -t54 * t135 + t27 * t89 + (t133 * t97 - t151) * qJD(5) + t109, t53 * t135 + t28 * t89 + (qJD(5) * t133 - t42) * t100 + t108, t2, t1, t4, t3, t165, t82 * t16 + t117 * t96 + t119 * t43 + (t118 * t96 - t123 * t99) * t177 + t111, t82 * t15 + t117 * t99 - t119 * t115 + (t118 * t99 + t123 * t96) * t177 + t121; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t161, t23, 0, t20, 0, t109 + (-qJD(5) + t89) * t26, t25 * t89 + t108 - t175, t2, t1, t4, t3, t165, -pkin(5) * t16 - (-t96 * t25 + t99 * t34) * t177 - t26 * t43 + t96 * t169 + (-t146 * t177 - t156) * pkin(9) + t111, -pkin(5) * t15 + (t99 * t25 + t96 * t34) * t177 + t26 * t115 + t99 * t169 + (t147 * t177 - t154) * pkin(9) + t121; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t115 * t43, t115 ^ 2 - t43 ^ 2, t15 + t180, -t16 - t179, t32, t21 * t115 - t178 * t9 - t96 * t6 + t10, -t96 * t11 + t116 * t178 + t21 * t43 - t99 * t6;];
tauc_reg  = t5;
