% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5PRPRR7
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
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRPRR7_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR7_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR7_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR7_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR7_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR7_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:00:52
% EndTime: 2019-12-05 16:00:56
% DurationCPUTime: 1.37s
% Computational Cost: add. (1633->249), mult. (3047->318), div. (0->0), fcn. (1992->10), ass. (0->149)
t138 = qJD(2) * qJD(4);
t90 = cos(qJ(4));
t142 = t90 * qJDD(2);
t87 = sin(qJ(4));
t192 = t87 * t138 - t142;
t84 = sin(pkin(8));
t85 = cos(pkin(8));
t185 = -g(1) * t85 - g(2) * t84;
t91 = cos(qJ(2));
t110 = t185 * t91;
t150 = qJDD(1) - g(3);
t88 = sin(qJ(2));
t191 = -t150 * t88 - t110;
t183 = pkin(2) + pkin(6);
t190 = t183 * qJDD(2);
t81 = t87 ^ 2;
t82 = t90 ^ 2;
t165 = t81 + t82;
t86 = sin(qJ(5));
t89 = cos(qJ(5));
t115 = t86 * t87 - t89 * t90;
t38 = t115 * t91;
t141 = t91 * qJDD(1);
t140 = qJD(1) * qJD(2);
t71 = t88 * t140;
t114 = qJDD(3) + t71 - t141;
t36 = t114 - t190;
t132 = t165 * t36;
t151 = t91 * qJD(1);
t125 = qJD(3) - t151;
t49 = -t183 * qJD(2) + t125;
t188 = t165 * t49;
t50 = t86 * t90 + t87 * t89;
t43 = t50 * qJD(2);
t78 = g(3) * t91;
t100 = t185 * t88 + t78;
t149 = qJD(2) * qJ(3);
t152 = t88 * qJD(1);
t57 = t149 + t152;
t153 = t57 * qJD(2);
t187 = t100 - t153;
t95 = qJD(2) ^ 2;
t186 = qJDD(2) * t91 - t88 * t95;
t80 = qJD(4) + qJD(5);
t158 = qJD(4) * t87;
t27 = t90 * t36;
t10 = qJDD(4) * pkin(4) + t192 * pkin(7) - t49 * t158 + t27;
t130 = t90 * t138;
t144 = t87 * qJDD(2);
t157 = qJD(4) * t90;
t13 = t49 * t157 + t87 * t36 + (-t130 - t144) * pkin(7);
t156 = qJD(5) * t86;
t160 = qJD(2) * t90;
t29 = -pkin(7) * t160 + t90 * t49;
t26 = qJD(4) * pkin(4) + t29;
t161 = qJD(2) * t87;
t28 = -pkin(7) * t161 + t49 * t87;
t1 = (qJD(5) * t26 + t13) * t89 + t86 * t10 - t28 * t156;
t184 = (t149 + t57 - t152) * qJD(4) - qJDD(4) * t183;
t182 = pkin(4) * t87;
t179 = g(3) * t88;
t178 = pkin(7) + t183;
t107 = t50 * t88;
t52 = t178 * t87;
t53 = t178 * t90;
t23 = t52 * t86 - t53 * t89;
t47 = t178 * t158;
t48 = qJD(4) * t53;
t177 = -qJD(1) * t107 + t23 * qJD(5) + t86 * t47 - t89 * t48;
t108 = t88 * t115;
t24 = -t52 * t89 - t53 * t86;
t176 = qJD(1) * t108 - t24 * qJD(5) + t89 * t47 + t86 * t48;
t175 = t28 * t86;
t174 = t28 * t89;
t134 = t86 * t161;
t45 = t89 * t160 - t134;
t173 = t45 * t43;
t172 = t57 * t91;
t171 = t84 * t88;
t170 = t85 * t88;
t21 = t80 * t50;
t79 = qJDD(4) + qJDD(5);
t168 = -t115 * t79 - t21 * t80;
t167 = t91 * pkin(2) + t88 * qJ(3);
t166 = t81 - t82;
t94 = qJD(4) ^ 2;
t164 = t94 + t95;
t163 = qJ(3) * t91;
t137 = qJDD(2) * qJ(3);
t143 = t88 * qJDD(1);
t37 = t137 + t143 + (qJD(3) + t151) * qJD(2);
t162 = t37 * qJ(3);
t159 = qJD(2) * t91;
t155 = qJDD(2) * pkin(2);
t70 = qJ(3) + t182;
t46 = t70 * qJD(2) + t152;
t154 = t46 * qJD(2);
t148 = qJDD(2) * t88;
t146 = qJDD(4) * t87;
t139 = qJD(2) * qJD(3);
t136 = t90 * t95 * t87;
t135 = -g(1) * t170 - g(2) * t171 + t78;
t133 = t37 * t88 - g(3);
t61 = pkin(4) * t157 + qJD(3);
t128 = t61 - t151;
t126 = t80 * t90;
t124 = t87 * t130;
t62 = t84 * t163;
t63 = t85 * t163;
t123 = -g(1) * t63 - g(2) * t62;
t121 = g(1) * t84 - g(2) * t85;
t120 = -t89 * t142 + t86 * t144;
t16 = t21 * qJD(2) + t120;
t119 = t115 * t16 - t21 * t45;
t112 = -qJD(5) * t134 - t192 * t86;
t17 = (qJD(2) * t126 + t144) * t89 + t112;
t22 = t89 * t126 - t87 * t156 - t86 * t158;
t118 = t17 * t50 + t22 * t43;
t117 = -t22 * t80 - t50 * t79;
t12 = t26 * t86 + t174;
t116 = (-qJD(2) * pkin(2) + t125) * t88 + t172;
t113 = t135 - t153;
t111 = -t135 + t141;
t109 = t121 * t87;
t105 = t164 * t91 + t148;
t104 = -qJDD(4) * t91 + 0.2e1 * t88 * t138;
t41 = t114 - t155;
t101 = t110 - t179;
t2 = -t12 * qJD(5) + t89 * t10 - t86 * t13;
t11 = t26 * t89 - t175;
t99 = -t1 * t50 + t11 * t21 + t115 * t2 - t12 * t22 - t135;
t83 = qJ(4) + qJ(5);
t74 = sin(t83);
t75 = cos(t83);
t98 = -g(1) * (t75 * t170 - t74 * t84) - g(2) * (t75 * t171 + t74 * t85) - t45 * t46 + t2 + t75 * t78;
t97 = -g(1) * (-t74 * t170 - t75 * t84) - g(2) * (-t74 * t171 + t75 * t85) + t43 * t46 - t74 * t78 - t1;
t96 = -t179 + t137 + t139 + t183 * t94 + t37 + (t185 - t140) * t91;
t92 = -pkin(7) - pkin(6);
t73 = qJDD(4) * t90;
t54 = t91 * t95 + t148;
t39 = t50 * t91;
t20 = t143 + t70 * qJDD(2) + (t61 + t151) * qJD(2);
t18 = -t43 ^ 2 + t45 ^ 2;
t15 = t29 * t89 - t175;
t14 = -t29 * t86 - t174;
t8 = -qJD(2) * t108 + t21 * t91;
t7 = qJD(2) * t107 + t80 * t38;
t4 = t45 * t80 + (-t80 * t160 - t144) * t89 - t112;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t150, 0, 0, 0, 0, 0, 0, t186, -t54, 0, -g(3) + (t88 ^ 2 + t91 ^ 2) * qJDD(1), 0, 0, 0, 0, 0, 0, 0, -t186, t54, t116 * qJD(2) - t41 * t91 + t133, 0, 0, 0, 0, 0, 0, t104 * t90 + t105 * t87, -t104 * t87 + t105 * t90, t186 * t165, -t91 * t132 + (t88 * t188 + t172) * qJD(2) + t133, 0, 0, 0, 0, 0, 0, t43 * t159 + t17 * t88 + t38 * t79 + t8 * t80, t45 * t159 - t16 * t88 + t39 * t79 - t7 * t80, t16 * t38 + t17 * t39 - t43 * t7 - t45 * t8, -t1 * t39 + t11 * t8 + t12 * t7 + t154 * t91 + t2 * t38 + t20 * t88 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t111, t191, 0, 0, qJDD(2), 0, 0, 0, 0, 0, 0, qJDD(3) - t111 - 0.2e1 * t155, 0.2e1 * t137 + 0.2e1 * t139 - t191, t162 + t57 * qJD(3) - t41 * pkin(2) - g(1) * (-pkin(2) * t170 + t63) - g(2) * (-pkin(2) * t171 + t62) - g(3) * t167 - t116 * qJD(1), qJDD(2) * t82 - 0.2e1 * t124, 0.2e1 * t166 * t138 - 0.2e1 * t87 * t142, -t87 * t94 + t73, qJDD(2) * t81 + 0.2e1 * t124, -t90 * t94 - t146, 0, t184 * t90 + t96 * t87, -t184 * t87 + t96 * t90, -t135 + t165 * (-t36 + t71 + t190), t162 - g(3) * (pkin(6) * t91 + t167) + t125 * t57 - t183 * t132 + (-qJD(1) * t188 - t185 * t183) * t88 + t123, t119, t115 * t17 + t16 * t50 + t21 * t43 - t22 * t45, t168, t118, t117, 0, t101 * t74 + t128 * t43 + t17 * t70 + t176 * t80 + t20 * t50 + t22 * t46 + t23 * t79, t101 * t75 - t115 * t20 + t128 * t45 - t16 * t70 - t177 * t80 - t21 * t46 - t24 * t79, t16 * t23 - t17 * t24 - t176 * t45 - t177 * t43 + t99, t1 * t24 + t2 * t23 + t20 * t70 - g(3) * (t182 * t88 - t91 * t92 + t167) + t128 * t46 + t177 * t12 + t176 * t11 + t123 + t185 * (t91 * t182 + (-pkin(2) + t92) * t88); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t95, t113 + t41, 0, 0, 0, 0, 0, 0, -t164 * t87 + t73, -t164 * t90 - t146, -t165 * qJDD(2), t132 + t113, 0, 0, 0, 0, 0, 0, -qJD(2) * t43 + t168, -qJD(2) * t45 + t117, -t118 - t119, -t99 - t154; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t136, -t166 * t95, t142, -t136, -t144, qJDD(4), t187 * t90 + t109 + t27, t121 * t90 + (-t36 - t187) * t87, 0, 0, t173, t18, -t120, -t173, t4, t79, -t14 * t80 + (-t80 * t156 - t43 * t160 + t79 * t89) * pkin(4) + t98, t15 * t80 + (-qJD(5) * t80 * t89 - t45 * t160 - t79 * t86) * pkin(4) + t97, (t12 + t14) * t45 + (-t11 + t15) * t43 + (t16 * t89 - t17 * t86 + (-t43 * t89 + t45 * t86) * qJD(5)) * pkin(4), -t11 * t14 - t12 * t15 + (t1 * t86 + t2 * t89 + t109 + (t100 - t154) * t90 + (-t11 * t86 + t12 * t89) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t173, t18, -t120, -t173, t4, t79, t12 * t80 + t98, t11 * t80 + t97, 0, 0;];
tau_reg = t3;
