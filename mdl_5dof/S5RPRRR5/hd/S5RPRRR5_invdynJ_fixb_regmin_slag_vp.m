% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPRRR5
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% tau_reg [5x21]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 09:49
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRRR5_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR5_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR5_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR5_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR5_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR5_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 09:49:02
% EndTime: 2022-01-20 09:49:06
% DurationCPUTime: 0.98s
% Computational Cost: add. (1295->176), mult. (2151->241), div. (0->0), fcn. (1414->14), ass. (0->136)
t193 = pkin(7) + pkin(8);
t103 = qJD(1) + qJD(3);
t109 = sin(qJ(5));
t114 = cos(qJ(4));
t110 = sin(qJ(4));
t113 = cos(qJ(5));
t164 = t113 * t110;
t52 = t109 * t114 + t164;
t45 = t52 * t103;
t107 = sin(pkin(9));
t184 = pkin(1) * t107;
t154 = qJD(1) * t184;
t108 = cos(pkin(9));
t89 = t108 * pkin(1) + pkin(2);
t192 = -qJD(3) * t154 + t89 * qJDD(1);
t94 = qJ(1) + pkin(9) + qJ(3);
t88 = cos(t94);
t186 = g(2) * t88;
t111 = sin(qJ(3));
t115 = cos(qJ(3));
t73 = t89 * qJD(1);
t190 = qJD(3) * t73 + qJDD(1) * t184;
t141 = -t190 * t111 + t192 * t115;
t101 = qJDD(1) + qJDD(3);
t183 = t101 * pkin(3);
t20 = -t141 - t183;
t191 = t20 + t186;
t87 = sin(t94);
t189 = g(1) * t88 + g(2) * t87;
t142 = -t111 * t184 + t115 * t89;
t40 = t111 * t73 + t115 * t154;
t145 = t193 * t103 + t40;
t25 = t114 * qJD(2) - t145 * t110;
t102 = qJD(4) + qJD(5);
t26 = t110 * qJD(2) + t145 * t114;
t188 = t192 * t111 + t190 * t115;
t83 = g(1) * t87;
t174 = t111 * t89 + t115 * t184;
t48 = pkin(7) + t174;
t185 = -pkin(8) - t48;
t182 = t103 * pkin(3);
t181 = t114 * pkin(4);
t163 = t113 * t114;
t150 = t103 * t163;
t167 = t109 * t110;
t151 = t103 * t167;
t43 = -t150 + t151;
t180 = t45 * t43;
t106 = qJ(4) + qJ(5);
t96 = sin(t106);
t179 = t87 * t96;
t97 = cos(t106);
t178 = t87 * t97;
t177 = t88 * t96;
t176 = t88 * t97;
t157 = t110 * qJD(4);
t39 = -t111 * t154 + t115 * t73;
t33 = -t39 - t182;
t175 = t114 * t83 + t33 * t157;
t173 = t113 * t26;
t171 = t40 * t103;
t42 = t174 * qJD(3);
t170 = t42 * t103;
t168 = t103 * t110;
t165 = t110 * t101;
t162 = t114 * t101;
t161 = qJDD(2) - g(3);
t104 = t110 ^ 2;
t160 = -t114 ^ 2 + t104;
t159 = qJD(5) * t109;
t156 = t114 * qJD(4);
t155 = t191 * t110 + t33 * t156;
t153 = pkin(4) * t157;
t92 = -pkin(3) - t181;
t148 = qJD(4) * t193;
t147 = t103 * t156;
t19 = t101 * pkin(7) + t188;
t146 = pkin(8) * t101 + t19;
t144 = qJD(4) * t185;
t47 = -pkin(3) - t142;
t139 = -t40 + t153;
t138 = t109 * t165 - t113 * t162;
t112 = sin(qJ(1));
t116 = cos(qJ(1));
t137 = g(1) * t112 - g(2) * t116;
t100 = qJDD(4) + qJDD(5);
t51 = -t163 + t167;
t30 = t102 * t51;
t21 = t52 * t100 - t30 * t102;
t24 = qJD(4) * pkin(4) + t25;
t135 = -t109 * t24 - t173;
t35 = t185 * t110;
t98 = t114 * pkin(8);
t36 = t114 * t48 + t98;
t134 = -t109 * t36 + t113 * t35;
t133 = t109 * t35 + t113 * t36;
t10 = (t103 * t157 - t162) * pkin(4) + t20;
t27 = t92 * t103 - t39;
t132 = -g(1) * t179 + g(2) * t177 + t10 * t52 - t27 * t30;
t31 = t102 * t52;
t131 = g(1) * t178 - g(2) * t176 + t10 * t51 + t27 * t31;
t79 = t114 * pkin(7) + t98;
t130 = qJD(5) * t79 - t110 * t39 + t114 * t148;
t78 = t193 * t110;
t129 = qJD(5) * t78 + t110 * t148 + t114 * t39;
t127 = t141 + t83 - t186;
t117 = qJD(4) ^ 2;
t126 = pkin(7) * t117 - t171 - t183;
t125 = t101 * t47 + t117 * t48 + t170;
t124 = -t33 * t103 + t189 - t19;
t123 = -pkin(7) * qJDD(4) + (t39 - t182) * qJD(4);
t41 = t142 * qJD(3);
t122 = -qJDD(4) * t48 + (t103 * t47 - t41) * qJD(4);
t13 = qJD(5) * t150 + t101 * t164 - t102 * t151 + t109 * t162 + t113 * t147;
t93 = t114 * qJDD(2);
t4 = qJDD(4) * pkin(4) - t26 * qJD(4) - t146 * t110 + t93;
t121 = t27 * t43 + t26 * t159 + g(2) * t178 + g(1) * t176 + g(3) * t96 + (-t26 * t102 - t4) * t109;
t120 = -t188 + t189;
t5 = t25 * qJD(4) + t110 * qJDD(2) + t146 * t114;
t119 = g(1) * t177 + g(2) * t179 - g(3) * t97 + t135 * qJD(5) - t109 * t5 + t113 * t4 - t27 * t45;
t99 = t103 ^ 2;
t64 = qJDD(4) * t114 - t117 * t110;
t63 = qJDD(4) * t110 + t117 * t114;
t46 = t104 * t101 + 0.2e1 * t110 * t147;
t38 = t47 - t181;
t37 = t42 + t153;
t32 = -0.2e1 * t160 * t103 * qJD(4) + 0.2e1 * t110 * t162;
t22 = -t51 * t100 - t31 * t102;
t18 = -t110 * t41 + t114 * t144;
t17 = t110 * t144 + t114 * t41;
t15 = -t43 ^ 2 + t45 ^ 2;
t14 = t31 * t103 + t138;
t6 = t43 * t102 + t13;
t2 = t13 * t52 - t45 * t30;
t1 = -t13 * t51 - t52 * t14 + t30 * t43 - t45 * t31;
t3 = [qJDD(1), t137, g(1) * t116 + g(2) * t112, (t137 + (t107 ^ 2 + t108 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), t101, t142 * t101 + t127 - t170, -t174 * t101 - t41 * t103 + t120, t46, t32, t63, t64, 0, t122 * t110 + (-t125 - t191) * t114 + t175, t122 * t114 + (t125 - t83) * t110 + t155, t2, t1, t21, t22, 0, t37 * t43 + t38 * t14 + (-t133 * qJD(5) - t109 * t17 + t113 * t18) * t102 + t134 * t100 + t131, t37 * t45 + t38 * t13 - (t134 * qJD(5) + t109 * t18 + t113 * t17) * t102 - t133 * t100 + t132; 0, 0, 0, t161, 0, 0, 0, 0, 0, 0, 0, 0, t64, -t63, 0, 0, 0, 0, 0, t22, -t21; 0, 0, 0, 0, t101, t127 + t171, t39 * t103 + t120, t46, t32, t63, t64, 0, t123 * t110 + (-t126 - t191) * t114 + t175, t123 * t114 + (t126 - t83) * t110 + t155, t2, t1, t21, t22, 0, t92 * t14 + (-t109 * t79 - t113 * t78) * t100 + t139 * t43 + (t129 * t109 - t130 * t113) * t102 + t131, t92 * t13 - (-t109 * t78 + t113 * t79) * t100 + t139 * t45 + (t130 * t109 + t129 * t113) * t102 + t132; 0, 0, 0, 0, 0, 0, 0, -t110 * t99 * t114, t160 * t99, t165, t162, qJDD(4), -g(3) * t114 + t124 * t110 + t93, -t161 * t110 + t124 * t114, t180, t15, t6, -t138, t100, -(-t109 * t25 - t173) * t102 + (t113 * t100 - t102 * t159 - t43 * t168) * pkin(4) + t119, (-qJD(5) * t24 + t25 * t102 - t5) * t113 + (-qJD(5) * t113 * t102 - t109 * t100 - t45 * t168) * pkin(4) + t121; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t180, t15, t6, -t138, t100, -t135 * t102 + t119, (-t5 + (-qJD(5) + t102) * t24) * t113 + t121;];
tau_reg = t3;
