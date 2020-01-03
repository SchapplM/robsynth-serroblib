% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% tauc_reg [5x23]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRPR5_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR5_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:43:02
% EndTime: 2020-01-03 11:43:10
% DurationCPUTime: 1.42s
% Computational Cost: add. (1638->192), mult. (4613->301), div. (0->0), fcn. (3368->8), ass. (0->139)
t114 = sin(qJ(5));
t116 = cos(qJ(5));
t149 = qJD(5) * t114;
t111 = sin(pkin(8));
t110 = sin(pkin(9));
t112 = cos(pkin(9));
t115 = sin(qJ(3));
t117 = cos(qJ(3));
t89 = t110 * t117 + t112 * t115;
t124 = qJD(1) * t89;
t72 = t111 * t124;
t67 = t116 * t72;
t128 = t110 * t115 - t112 * t117;
t153 = qJD(3) * t111;
t74 = t128 * t153;
t70 = qJD(1) * t74;
t78 = t89 * t153;
t71 = qJD(1) * t78;
t156 = qJD(1) * t111;
t142 = t117 * t156;
t143 = t115 * t156;
t75 = -t110 * t143 + t112 * t142;
t11 = -qJD(5) * t67 + t114 * t70 - t116 * t71 - t149 * t75;
t36 = t114 * t75 + t67;
t113 = cos(pkin(8));
t148 = t113 * qJD(1);
t99 = -qJD(3) + t148;
t95 = -qJD(5) + t99;
t174 = t36 * t95;
t190 = t11 - t174;
t129 = -t114 * t72 + t116 * t75;
t189 = t129 * t36;
t12 = qJD(5) * t129 - t114 * t71 - t116 * t70;
t175 = t129 * t95;
t188 = -t12 - t175;
t187 = t129 ^ 2 - t36 ^ 2;
t84 = pkin(3) * t143 + qJ(2) * t156 + qJD(4);
t52 = t72 * pkin(4) + t84;
t178 = t72 * pkin(7);
t165 = qJ(2) * t117;
t146 = t113 * t165;
t90 = -t113 * pkin(2) - t111 * pkin(6) - pkin(1);
t83 = qJD(1) * t90 + qJD(2);
t54 = -qJ(4) * t143 + qJD(1) * t146 + t115 * t83;
t169 = t112 * t54;
t161 = t113 * t115;
t164 = qJ(4) * t111;
t123 = -qJ(2) * t161 - t117 * t164;
t80 = t117 * t83;
t53 = qJD(1) * t123 + t80;
t44 = -t99 * pkin(3) + t53;
t21 = t110 * t44 + t169;
t10 = t21 - t178;
t9 = t10 * t149;
t186 = t52 * t36 + t9;
t106 = t111 ^ 2;
t184 = 0.2e1 * t106;
t183 = t115 * t117;
t182 = qJD(5) + t95;
t107 = t113 ^ 2;
t181 = t184 + t107;
t150 = qJD(4) * t111;
t119 = qJD(3) * t123 - t115 * t150;
t147 = qJD(1) * qJD(2);
t140 = t113 * t147;
t151 = qJD(3) * t117;
t171 = t117 * t140 + t83 * t151;
t30 = qJD(1) * t119 + t171;
t155 = qJD(2) * t115;
t141 = t113 * t155;
t122 = -t117 * t150 - t141;
t152 = qJD(3) * t115;
t162 = t111 * t115;
t31 = -t83 * t152 + ((qJ(4) * t162 - t146) * qJD(3) + t122) * qJD(1);
t4 = -t110 * t30 + t112 * t31;
t2 = t71 * pkin(7) + t4;
t5 = t110 * t31 + t112 * t30;
t3 = t70 * pkin(7) + t5;
t145 = -t114 * t3 + t116 * t2;
t180 = -t52 * t129 + t145;
t137 = -t90 + t164;
t179 = t115 * t137 - t146;
t177 = t75 * pkin(7);
t176 = pkin(3) * t110;
t154 = qJD(2) * t117;
t170 = t113 * t154 + t90 * t151;
t42 = t119 + t170;
t43 = t179 * qJD(3) + t122;
t16 = t110 * t43 + t112 * t42;
t48 = t110 * t54;
t24 = t112 * t53 - t48;
t166 = qJ(2) * t115;
t58 = -t137 * t117 + (-pkin(3) - t166) * t113;
t26 = t110 * t58 - t112 * t179;
t173 = t89 * qJD(3) - t113 * t124;
t172 = t99 * t128;
t134 = pkin(3) * t142;
t168 = qJD(3) * t134 + t111 * t147;
t167 = (pkin(3) * t151 + qJD(2)) * t111;
t118 = qJD(1) ^ 2;
t163 = t106 * t118;
t160 = qJD(3) + t99;
t159 = pkin(3) * t162 + t111 * qJ(2);
t158 = t106 + t107;
t157 = t115 ^ 2 - t117 ^ 2;
t144 = qJ(2) * t152;
t139 = qJD(1) * qJD(3) * t106;
t15 = -t110 * t42 + t112 * t43;
t20 = t112 * t44 - t48;
t23 = -t110 * t53 - t169;
t25 = t110 * t179 + t112 * t58;
t136 = t158 * t118;
t135 = qJD(1) * t160;
t133 = qJD(5) * t89 + t173;
t132 = -qJD(5) * t128 + t172;
t131 = t111 * t135;
t8 = -t99 * pkin(4) - t177 + t20;
t130 = -t116 * t10 - t114 * t8;
t81 = t89 * t111;
t82 = t128 * t111;
t45 = -t114 * t82 + t116 * t81;
t46 = -t114 * t81 - t116 * t82;
t127 = 0.2e1 * t158 * t147;
t126 = (t99 + t148) * t153;
t121 = -t99 ^ 2 - t163;
t103 = t112 * pkin(3) + pkin(4);
t61 = t75 * pkin(4) + t134;
t59 = t81 * pkin(4) + t159;
t55 = -t74 * pkin(4) + t167;
t47 = -t70 * pkin(4) + t168;
t22 = -t81 * pkin(7) + t26;
t19 = -t113 * pkin(4) + t82 * pkin(7) + t25;
t18 = qJD(5) * t46 - t114 * t78 - t116 * t74;
t17 = -qJD(5) * t45 + t114 * t74 - t116 * t78;
t14 = t24 - t177;
t13 = t23 + t178;
t7 = t74 * pkin(7) + t16;
t6 = t78 * pkin(7) + t15;
t1 = [0, 0, 0, 0, 0, t127, qJ(2) * t127, -0.2e1 * t139 * t183, 0.2e1 * t157 * t139, t115 * t126, t117 * t126, 0, t99 * t141 + (-(-t115 * t90 - t146) * t99 + t83 * t161) * qJD(3) + t181 * qJD(1) * (qJ(2) * t151 + t155), (-t113 * t144 + t170) * t99 + t171 * t113 + (-t181 * t144 + t154 * t184) * qJD(1), -t15 * t75 - t16 * t72 + t20 * t78 + t21 * t74 + t25 * t71 + t26 * t70 + t4 * t82 - t5 * t81, t20 * t15 + t159 * t168 + t21 * t16 + t167 * t84 + t4 * t25 + t5 * t26, t11 * t46 + t129 * t17, -t11 * t45 - t46 * t12 - t129 * t18 - t17 * t36, -t11 * t113 - t17 * t95, t12 * t113 + t18 * t95, 0, -(-t114 * t7 + t116 * t6) * t95 - t145 * t113 + t55 * t36 + t59 * t12 + t47 * t45 + t52 * t18 + (-(-t114 * t19 - t116 * t22) * t95 - t130 * t113) * qJD(5), t59 * t11 - t9 * t113 + t52 * t17 + t55 * t129 + t47 * t46 + ((-qJD(5) * t22 + t6) * t95 + t2 * t113) * t114 + ((qJD(5) * t19 + t7) * t95 + (qJD(5) * t8 + t3) * t113) * t116; 0, 0, 0, 0, 0, -t136, -qJ(2) * t136, 0, 0, 0, 0, 0, t121 * t115, t121 * t117, -t128 * t71 - t172 * t72 + t173 * t75 + t89 * t70, -t128 * t4 - t84 * t156 + t172 * t21 - t173 * t20 + t5 * t89, 0, 0, 0, 0, 0, -t36 * t156 + (t114 * t132 + t116 * t133) * t95, -t129 * t156 + (-t114 * t133 + t116 * t132) * t95; 0, 0, 0, 0, 0, 0, 0, t163 * t183, -t157 * t163, -t115 * t131, -t117 * t131, 0, (-t160 * t83 - t140) * t115 + (-t113 * t135 - t163) * t165, -t80 * t99 + (t148 * t160 + t163) * t166 - t171, (t21 + t23) * t75 - (t20 - t24) * t72 + (t110 * t70 + t112 * t71) * pkin(3), -t20 * t23 - t21 * t24 + (t110 * t5 + t112 * t4 - t142 * t84) * pkin(3), t189, t187, t190, t188, 0, (-t114 * t14 + t116 * t13) * t95 - t61 * t36 + (-(-t103 * t114 - t116 * t176) * t95 + t130) * qJD(5) + t180, -t116 * t3 - t114 * t2 - (t114 * t13 + t116 * t14) * t95 - t61 * t129 + ((t103 * t116 - t114 * t176) * t95 - t116 * t8) * qJD(5) + t186; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t72 ^ 2 - t75 ^ 2, t20 * t75 + t21 * t72 + t168, 0, 0, 0, 0, 0, t12 - t175, t11 + t174; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t189, t187, t190, t188, 0, t182 * t130 + t180, (t10 * t95 - t2) * t114 + (-t182 * t8 - t3) * t116 + t186;];
tauc_reg = t1;
