% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% 
% Output:
% tauc_reg [5x25]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRRPP5_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP5_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP5_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP5_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:58:28
% EndTime: 2019-12-31 20:58:32
% DurationCPUTime: 1.11s
% Computational Cost: add. (1744->221), mult. (4427->259), div. (0->0), fcn. (2857->4), ass. (0->135)
t111 = sin(qJ(3));
t112 = sin(qJ(2));
t113 = cos(qJ(2));
t174 = cos(qJ(3));
t80 = t111 * t113 + t174 * t112;
t155 = qJD(1) * t80;
t141 = t113 * pkin(2) + pkin(1);
t89 = t141 * qJD(1);
t183 = -t155 * qJ(4) - t89;
t108 = qJD(2) + qJD(3);
t145 = qJD(1) * qJD(2);
t181 = -0.2e1 * t145;
t105 = t108 * qJD(4);
t136 = t174 * qJD(3);
t95 = pkin(2) * t136 + qJD(4);
t180 = t95 * t108 + t105;
t64 = t155 ^ 2;
t179 = -t108 ^ 2 - t64;
t176 = -pkin(7) - pkin(6);
t91 = t176 * t113;
t86 = qJD(1) * t91;
t162 = t111 * t86;
t164 = qJD(2) * pkin(2);
t90 = t176 * t112;
t84 = qJD(1) * t90;
t75 = t84 + t164;
t40 = t174 * t75 + t162;
t178 = qJD(4) - t40;
t139 = t174 * t113;
t131 = qJD(1) * t139;
t147 = qJD(1) * t112;
t137 = t111 * t147;
t65 = -t131 + t137;
t177 = t65 ^ 2;
t114 = -pkin(3) - pkin(4);
t175 = t155 * pkin(4);
t26 = t65 * pkin(3) + t183;
t173 = t26 * t65;
t172 = t26 * t155;
t171 = t155 * t65;
t170 = t89 * t65;
t169 = t89 * t155;
t43 = t174 * t84 + t162;
t60 = t155 * qJ(5);
t24 = t60 + t43;
t168 = -t24 + t95;
t167 = -t43 + t95;
t38 = pkin(3) * t155 + t65 * qJ(4);
t74 = t174 * t86;
t41 = t111 * t75 - t74;
t47 = t111 * t90 - t174 * t91;
t166 = t108 * t131;
t45 = t108 * t80;
t37 = t45 * qJD(1);
t165 = qJ(4) * t37;
t101 = t111 * pkin(2) + qJ(4);
t163 = t101 * t37;
t146 = qJD(3) * t111;
t140 = qJD(2) * t176;
t85 = t112 * t140;
t87 = t113 * t140;
t13 = t111 * t87 + t90 * t136 + t91 * t146 + t174 * t85;
t161 = t13 * t108;
t14 = t47 * qJD(3) + t111 * t85 - t174 * t87;
t160 = t14 * t108;
t159 = t65 * qJ(5);
t158 = t65 * t108;
t154 = t111 * t112;
t116 = qJD(1) ^ 2;
t153 = t113 * t116;
t115 = qJD(2) ^ 2;
t152 = t115 * t112;
t151 = t115 * t113;
t20 = t60 + t40;
t150 = qJD(4) - t20;
t148 = t112 ^ 2 - t113 ^ 2;
t144 = pkin(2) * t146;
t143 = t112 * t164;
t138 = t108 * t146;
t135 = t112 * t145;
t132 = qJD(1) * t140;
t76 = t112 * t132;
t77 = t113 * t132;
t134 = -t111 * t77 - t75 * t136 - t86 * t146 - t174 * t76;
t10 = t111 * t76 - t86 * t136 + t75 * t146 - t174 * t77;
t133 = pkin(1) * t181;
t103 = -t174 * pkin(2) - pkin(3);
t21 = t159 + t41;
t42 = t111 * t84 - t74;
t23 = t42 + t159;
t130 = -t23 + t144;
t129 = -t42 + t144;
t46 = -t111 * t91 - t174 * t90;
t128 = t108 * t154;
t127 = t37 * qJ(5) + t65 * qJD(5) - t134;
t36 = qJD(1) * t128 - t166;
t3 = t36 * qJ(5) - qJD(5) * t155 + t10;
t28 = pkin(2) * t147 + t38;
t126 = t80 * qJ(4) + t141;
t125 = t40 * t108 + t134;
t124 = t41 * t108 - t10;
t123 = t42 * t108 - t10;
t122 = t43 * t108 + t134;
t7 = pkin(2) * t135 + t37 * pkin(3) + t36 * qJ(4) - qJD(4) * t155;
t12 = t114 * t65 + qJD(5) - t183;
t121 = -t12 * t155 + t3;
t120 = t12 * t65 + t127;
t119 = -t108 * t137 + t166;
t44 = -qJD(2) * t139 - t113 * t136 + t128;
t118 = -t44 * qJ(4) + t80 * qJD(4) - t143;
t1 = -t37 * pkin(4) - t7;
t106 = t108 * qJ(4);
t100 = 0.2e1 * t105;
t99 = -pkin(4) + t103;
t92 = pkin(2) * t138;
t79 = -t139 + t154;
t39 = t79 * pkin(3) - t126;
t31 = t106 + t41;
t30 = t79 * qJ(5) + t47;
t29 = -t80 * qJ(5) + t46;
t27 = -t108 * pkin(3) + t178;
t25 = t64 - t177;
t22 = t114 * t79 + t126;
t19 = -t38 - t175;
t17 = t119 + t158;
t16 = t106 + t21;
t15 = -t28 - t175;
t11 = t114 * t108 + t178 - t60;
t9 = t105 - t134;
t8 = t45 * pkin(3) - t118;
t6 = t44 * qJ(5) - t80 * qJD(5) + t14;
t5 = t45 * qJ(5) + t79 * qJD(5) + t13;
t4 = t114 * t45 + t118;
t2 = t105 + t127;
t18 = [0, 0, 0, 0.2e1 * t113 * t135, t148 * t181, t151, -t152, 0, -pkin(6) * t151 + t112 * t133, pkin(6) * t152 + t113 * t133, -t155 * t44 - t36 * t80, -t155 * t45 + t36 * t79 - t80 * t37 + t44 * t65, -t44 * t108, -t45 * t108, 0, -t141 * t37 - t160 - t89 * t45 + (qJD(1) * t79 + t65) * t143, t141 * t36 + 0.2e1 * t155 * t143 + t89 * t44 - t161, t26 * t45 + t39 * t37 + t8 * t65 + t7 * t79 - t160, t10 * t80 - t13 * t65 + t14 * t155 - t27 * t44 - t31 * t45 - t46 * t36 - t47 * t37 - t9 * t79, -t155 * t8 + t26 * t44 + t39 * t36 - t7 * t80 + t161, t10 * t46 + t31 * t13 + t27 * t14 + t26 * t8 + t7 * t39 + t9 * t47, -t1 * t79 - t6 * t108 - t12 * t45 - t22 * t37 - t4 * t65, t1 * t80 + t5 * t108 - t12 * t44 + t155 * t4 - t22 * t36, t11 * t44 - t155 * t6 + t16 * t45 + t2 * t79 + t29 * t36 - t3 * t80 + t30 * t37 + t5 * t65, t1 * t22 + t11 * t6 + t12 * t4 + t16 * t5 + t2 * t30 + t3 * t29; 0, 0, 0, -t112 * t153, t148 * t116, 0, 0, 0, t116 * pkin(1) * t112, pkin(1) * t153, t171, t25, t17, 0, 0, t169 + (-t65 * t147 - t138) * pkin(2) + t123, -t170 + (-t108 * t136 - t147 * t155) * pkin(2) + t122, -t28 * t65 + t123 - t172 - t92, -t163 - t103 * t36 + (t129 + t31) * t155 + (t27 - t167) * t65, t155 * t28 - t122 - t173 + t180, t10 * t103 + t9 * t101 + t129 * t27 + t167 * t31 - t26 * t28, t23 * t108 + t15 * t65 - t121 - t92, -t24 * t108 - t15 * t155 + t120 + t180, t163 + t99 * t36 + (-t130 - t16) * t155 + (-t11 + t168) * t65, t2 * t101 + t11 * t130 - t12 * t15 + t16 * t168 + t3 * t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t171, t25, t17, 0, 0, t124 + t169, t125 - t170, -t38 * t65 + t124 - t172, pkin(3) * t36 - t165 + (t31 - t41) * t155 + (t27 - t178) * t65, t155 * t38 + t100 - t125 - t173, -t10 * pkin(3) + t9 * qJ(4) + t178 * t31 - t26 * t38 - t27 * t41, t21 * t108 + t19 * t65 - t121, -t20 * t108 - t155 * t19 + t100 + t120, t165 + t114 * t36 + (-t16 + t21) * t155 + (-t11 + t150) * t65, t2 * qJ(4) - t11 * t21 + t3 * t114 - t12 * t19 + t150 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t171, t17, t179, -t31 * t108 + t10 + t172, t171, t179, t36 - t158, -t16 * t108 + t121; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t108 * t155, t119 - t158, -t64 - t177, t11 * t155 - t16 * t65 + t1;];
tauc_reg = t18;
