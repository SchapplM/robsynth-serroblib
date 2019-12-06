% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5PRRPR4
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% tau_reg [5x20]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:24
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRRPR4_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR4_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR4_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR4_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR4_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR4_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:23:25
% EndTime: 2019-12-05 16:23:30
% DurationCPUTime: 1.34s
% Computational Cost: add. (1434->231), mult. (3295->328), div. (0->0), fcn. (2524->12), ass. (0->139)
t180 = qJ(4) + pkin(6);
t107 = qJD(3) + qJD(5);
t115 = sin(qJ(5));
t118 = cos(qJ(5));
t110 = sin(pkin(9));
t112 = cos(pkin(9));
t116 = sin(qJ(3));
t119 = cos(qJ(3));
t140 = t110 * t116 - t112 * t119;
t79 = t140 * qJD(2);
t69 = t118 * t79;
t88 = t110 * t119 + t112 * t116;
t81 = t88 * qJD(2);
t36 = -t115 * t81 - t69;
t173 = t36 * t107;
t165 = qJD(5) * t115;
t80 = t88 * qJD(3);
t41 = -qJD(2) * t80 - t140 * qJDD(2);
t161 = qJD(2) * qJD(3);
t154 = t119 * t161;
t155 = t116 * t161;
t42 = t88 * qJDD(2) - t110 * t155 + t112 * t154;
t4 = -qJD(5) * t69 + t115 * t41 + t118 * t42 - t81 * t165;
t201 = t4 - t173;
t142 = -t115 * t79 + t118 * t81;
t200 = t142 * t36;
t174 = t142 * t107;
t5 = t142 * qJD(5) + t115 * t42 - t118 * t41;
t199 = -t5 + t174;
t120 = cos(qJ(2));
t111 = sin(pkin(8));
t113 = cos(pkin(8));
t147 = g(1) * t113 + g(2) * t111;
t133 = t147 * t120;
t117 = sin(qJ(2));
t182 = g(3) * t117;
t126 = t133 + t182;
t181 = g(3) * t120;
t193 = t147 * t117;
t198 = t193 - t181;
t169 = qJDD(1) - g(3);
t197 = t169 * t120 + t193;
t196 = t142 ^ 2 - t36 ^ 2;
t105 = qJ(3) + pkin(9) + qJ(5);
t100 = sin(t105);
t101 = cos(t105);
t187 = t79 * pkin(7);
t164 = t117 * qJD(1);
t151 = t180 * qJD(2) + t164;
t75 = t151 * t119;
t175 = t112 * t75;
t176 = qJD(3) * pkin(3);
t74 = t151 * t116;
t61 = -t74 + t176;
t24 = t110 * t61 + t175;
t13 = t24 - t187;
t170 = t113 * t120;
t171 = t111 * t120;
t104 = t119 * pkin(3) + pkin(2);
t163 = t120 * qJD(1);
t84 = -t104 * qJD(2) + qJD(4) - t163;
t45 = t79 * pkin(4) + t84;
t195 = -t45 * t36 - g(1) * (-t111 * t100 - t101 * t170) - g(2) * (t113 * t100 - t101 * t171) + t13 * t165 + t101 * t182;
t152 = qJD(3) * t180;
t76 = t119 * qJD(4) - t116 * t152;
t77 = -t116 * qJD(4) - t119 * t152;
t179 = -t110 * t76 + t112 * t77 + t88 * t163;
t178 = t110 * t77 + t112 * t76 + t140 * t163;
t191 = t116 * t176 - t164;
t83 = t140 * qJD(3);
t190 = qJD(5) - t107;
t162 = qJD(1) * qJD(2);
t86 = qJDD(2) * pkin(6) + t117 * qJDD(1) + t120 * t162;
t131 = qJ(4) * qJDD(2) + qJD(2) * qJD(4) + t86;
t139 = qJD(3) * t151;
t21 = qJDD(3) * pkin(3) - t131 * t116 - t119 * t139;
t22 = -t116 * t139 + t131 * t119;
t6 = -t110 * t22 + t112 * t21;
t2 = qJDD(3) * pkin(4) - t42 * pkin(7) + t6;
t7 = t110 * t21 + t112 * t22;
t3 = t41 * pkin(7) + t7;
t189 = -t45 * t142 - g(1) * (-t100 * t170 + t111 * t101) - g(2) * (-t100 * t171 - t113 * t101) - t115 * t3 + t118 * t2 + t100 * t182;
t103 = t117 * t162;
t121 = qJD(3) ^ 2;
t158 = t120 * qJDD(1);
t188 = 0.2e1 * qJDD(2) * pkin(2) - pkin(6) * t121 + (t147 + t162) * t117 - t103 + t158 - t181;
t186 = t81 * pkin(7);
t185 = pkin(3) * t110;
t57 = t110 * t75;
t26 = -t112 * t74 - t57;
t94 = t180 * t116;
t95 = t180 * t119;
t48 = -t110 * t94 + t112 * t95;
t177 = qJD(2) * pkin(2);
t108 = t116 ^ 2;
t168 = -t119 ^ 2 + t108;
t122 = qJD(2) ^ 2;
t167 = t121 + t122;
t166 = qJD(2) * t117;
t160 = qJDD(3) * t116;
t159 = t119 * qJDD(2);
t157 = t120 * qJDD(2);
t23 = t112 * t61 - t57;
t25 = t110 * t74 - t175;
t47 = -t110 * t95 - t112 * t94;
t150 = t80 * pkin(4) + t191;
t32 = -pkin(7) * t140 + t48;
t149 = -t83 * pkin(7) + qJD(5) * t32 - t179;
t31 = -t88 * pkin(7) + t47;
t148 = t80 * pkin(7) - qJD(5) * t31 - t178;
t146 = g(1) * t111 - g(2) * t113;
t11 = qJD(3) * pkin(4) - t186 + t23;
t145 = -t115 * t11 - t118 * t13;
t72 = t88 * t117;
t73 = t140 * t117;
t144 = t115 * t73 - t118 * t72;
t143 = -t115 * t72 - t118 * t73;
t43 = t115 * t88 + t118 * t140;
t44 = -t115 * t140 + t118 * t88;
t102 = t112 * pkin(3) + pkin(4);
t138 = t115 * t102 + t118 * t185;
t137 = t118 * t102 - t115 * t185;
t132 = t146 * t119;
t127 = pkin(3) * t155 - t104 * qJDD(2) + qJDD(4) + t103;
t98 = -t163 - t177;
t124 = -pkin(6) * qJDD(3) + (t163 + t98 - t177) * qJD(3);
t46 = t127 - t158;
t123 = -t98 * qJD(2) + t126 - t86;
t106 = qJDD(3) + qJDD(5);
t65 = pkin(4) * t140 - t104;
t51 = t116 * qJD(2) * pkin(3) + t81 * pkin(4);
t28 = -t117 * t80 - t120 * t79;
t27 = t117 * t83 - t120 * t81;
t15 = t26 - t186;
t14 = t25 + t187;
t12 = -t41 * pkin(4) + t46;
t9 = t44 * qJD(5) - t115 * t83 + t118 * t80;
t8 = -t43 * qJD(5) - t115 * t80 - t118 * t83;
t1 = [t169, 0, -t122 * t117 + t157, -qJDD(2) * t117 - t122 * t120, 0, 0, 0, 0, 0, (-0.2e1 * t155 + t159) * t120 + (-t167 * t119 - t160) * t117, (-qJDD(3) * t117 - 0.2e1 * t120 * t161) * t119 + (t167 * t117 - t157) * t116, -t27 * t81 - t28 * t79 - t73 * t41 + t72 * t42, -t46 * t120 + t84 * t166 + t23 * t27 + t24 * t28 - t6 * t72 - t7 * t73 - g(3), 0, 0, 0, 0, 0, (-qJD(5) * t143 - t115 * t28 + t118 * t27) * t107 + t144 * t106 - t36 * t166 - t120 * t5, -(qJD(5) * t144 + t115 * t27 + t118 * t28) * t107 - t143 * t106 + t142 * t166 - t120 * t4; 0, qJDD(2), t197, -t169 * t117 + t133, t108 * qJDD(2) + 0.2e1 * t116 * t154, 0.2e1 * t116 * t159 - 0.2e1 * t168 * t161, t121 * t119 + t160, qJDD(3) * t119 - t121 * t116, 0, t124 * t116 + t188 * t119, -t188 * t116 + t124 * t119, -t140 * t7 - t178 * t79 - t179 * t81 + t23 * t83 - t24 * t80 + t48 * t41 - t47 * t42 - t6 * t88 - t126, t7 * t48 + t6 * t47 - t46 * t104 - g(3) * (t120 * t104 + t117 * t180) + t191 * t84 + t178 * t24 + t179 * t23 + t147 * (t104 * t117 - t120 * t180), t142 * t8 + t4 * t44, -t142 * t9 + t36 * t8 - t4 * t43 - t44 * t5, t44 * t106 + t8 * t107, -t43 * t106 - t9 * t107, 0, (-t115 * t32 + t118 * t31) * t106 + t65 * t5 + t12 * t43 + t45 * t9 - t150 * t36 + (t115 * t148 - t118 * t149) * t107 + t198 * t101, -(t115 * t31 + t118 * t32) * t106 + t65 * t4 + t12 * t44 + t45 * t8 + t150 * t142 + (t115 * t149 + t118 * t148) * t107 - t198 * t100; 0, 0, 0, 0, -t116 * t122 * t119, t168 * t122, t116 * qJDD(2), t159, qJDD(3), t123 * t116 - t132, t146 * t116 + t123 * t119, (t24 + t25) * t81 - (t23 - t26) * t79 + (t110 * t41 - t112 * t42) * pkin(3), -t23 * t25 - t24 * t26 + (t7 * t110 + t6 * t112 - t132 + (-t84 * qJD(2) + t126) * t116) * pkin(3), -t200, t196, t201, t199, t106, t137 * t106 - (-t115 * t15 + t118 * t14) * t107 + t51 * t36 + (-t107 * t138 + t145) * qJD(5) + t189, -t138 * t106 - t118 * t3 - t115 * t2 + (t115 * t14 + t118 * t15) * t107 - t51 * t142 + (-t107 * t137 - t118 * t11) * qJD(5) + t195; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t79 ^ 2 - t81 ^ 2, t23 * t81 + t24 * t79 + t127 - t197, 0, 0, 0, 0, 0, t5 + t174, t4 + t173; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t200, t196, t201, t199, t106, t190 * t145 + t189, (-t13 * t107 - t2) * t115 + (-t190 * t11 - t3) * t118 + t195;];
tau_reg = t1;
