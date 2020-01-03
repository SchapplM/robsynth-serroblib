% Calculate inertial parameters regressor of coriolis matrix for
% S4RRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% cmat_reg [(4*4)x(4*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:05
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4RRPR6_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR6_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR6_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR6_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:04:57
% EndTime: 2019-12-31 17:05:00
% DurationCPUTime: 1.54s
% Computational Cost: add. (2422->112), mult. (4771->167), div. (0->0), fcn. (5308->6), ass. (0->111)
t191 = qJD(2) + qJD(4);
t104 = sin(qJ(4));
t178 = cos(qJ(4));
t103 = sin(pkin(7));
t165 = cos(pkin(7));
t105 = sin(qJ(2));
t173 = -qJ(3) - pkin(5);
t95 = t173 * t105;
t106 = cos(qJ(2));
t96 = t173 * t106;
t111 = t103 * t96 + t165 * t95;
t90 = t103 * t106 + t165 * t105;
t55 = -t90 * pkin(6) + t111;
t68 = -t103 * t95 + t165 * t96;
t88 = t103 * t105 - t165 * t106;
t56 = -t88 * pkin(6) - t68;
t30 = t104 * t55 + t178 * t56;
t192 = t191 * t30;
t29 = -t104 * t56 + t178 * t55;
t190 = t191 * t29;
t169 = t104 * t88;
t84 = t178 * t90;
t181 = t84 - t169;
t174 = t181 ^ 2;
t63 = t104 * t90 + t178 * t88;
t175 = t63 ^ 2;
t186 = -t174 + t175;
t189 = t186 * qJD(1);
t144 = t63 * qJD(4);
t31 = -t63 * qJD(2) - t144;
t146 = t181 * qJD(1);
t184 = t63 * t146;
t183 = qJD(3) * t63;
t182 = t63 * qJD(1);
t180 = t90 ^ 2;
t117 = t84 / 0.2e1;
t179 = t90 * pkin(3);
t177 = pkin(2) * t103;
t176 = t105 * pkin(2);
t172 = qJD(2) * pkin(2);
t101 = -t106 * pkin(2) - pkin(1);
t72 = t88 * pkin(3) + t101;
t73 = t176 + t179;
t5 = t72 * t73;
t167 = t5 * qJD(1);
t7 = -t181 * t29 - t30 * t63;
t166 = t7 * qJD(1);
t164 = qJD(1) * t72;
t15 = t174 + t175;
t162 = t15 * qJD(1);
t116 = t165 * pkin(2) + pkin(3);
t81 = t104 * t177 - t178 * t116;
t82 = t104 * t116 + t178 * t177;
t112 = -t181 * t81 / 0.2e1 + t63 * t82 / 0.2e1;
t102 = t176 / 0.2e1;
t118 = t102 + t179 / 0.2e1;
t17 = t112 + t118;
t160 = t17 * qJD(1);
t19 = t101 * t176;
t159 = t19 * qJD(1);
t20 = t181 * t72 + t73 * t63;
t158 = t20 * qJD(1);
t21 = t181 * t73 - t63 * t72;
t157 = t21 * qJD(1);
t32 = -t111 * t90 + t68 * t88;
t155 = t32 * qJD(1);
t41 = 0.2e1 * t117 - t169;
t153 = t41 * qJD(1);
t109 = -t103 * t88 / 0.2e1 - t165 * t90 / 0.2e1;
t52 = (-t105 / 0.2e1 + t109) * pkin(2);
t152 = t52 * qJD(1);
t87 = t88 ^ 2;
t53 = t87 - t180;
t151 = t53 * qJD(1);
t57 = t101 * t90 + t88 * t176;
t150 = t57 * qJD(1);
t58 = -t101 * t88 + t90 * t176;
t149 = t58 * qJD(1);
t61 = t117 - t84 / 0.2e1;
t148 = t61 * qJD(1);
t147 = t61 * qJD(4);
t145 = t181 * qJD(2);
t141 = t181 * qJD(4);
t67 = t87 + t180;
t140 = t67 * qJD(1);
t139 = t88 * qJD(1);
t138 = t90 * qJD(1);
t137 = t90 * qJD(2);
t98 = -t105 ^ 2 + t106 ^ 2;
t136 = t98 * qJD(1);
t135 = qJD(1) * t106;
t134 = t105 * qJD(2);
t133 = t106 * qJD(2);
t132 = pkin(1) * t105 * qJD(1);
t131 = pkin(1) * t135;
t128 = t181 * t182;
t127 = t63 * t164;
t126 = t181 * t164;
t125 = t88 * t138;
t124 = t88 * t137;
t120 = t105 * t133;
t115 = t82 * qJD(2);
t114 = t81 * qJD(2);
t113 = t41 * qJD(4) + t145;
t99 = t105 * t135;
t86 = t88 * qJD(2);
t75 = t82 * qJD(4);
t74 = t81 * qJD(4);
t51 = t109 * pkin(2) + t102;
t18 = -t112 + t118;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t120, t98 * qJD(2), 0, -t120, 0, 0, -pkin(1) * t134, -pkin(1) * t133, 0, 0, -t124, t53 * qJD(2), 0, t124, 0, 0, t57 * qJD(2), t58 * qJD(2), t67 * qJD(3), t19 * qJD(2) + t32 * qJD(3), t31 * t181, t191 * t186, 0, (t141 + t145) * t63, 0, 0, t20 * qJD(2) + t72 * t141, t21 * qJD(2) - t144 * t72, t15 * qJD(3), t5 * qJD(2) + t7 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t99, t136, t133, -t99, -t134, 0, -pkin(5) * t133 - t132, pkin(5) * t134 - t131, 0, 0, -t125, t151, -t86, t125, -t137, 0, t68 * qJD(2) + t150, -qJD(2) * t111 + t149, (-t103 * t90 + t165 * t88) * t172, t159 + (t103 * t111 + t165 * t68) * t172 + t51 * qJD(3), -t128, t189, t31, t184, -t113, 0, t158 - t192, t157 - t190, (-t181 * t82 - t63 * t81) * qJD(2), t167 + (t29 * t82 + t30 * t81) * qJD(2) + t18 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t140, t51 * qJD(2) + t155, 0, 0, 0, 0, 0, 0, t147, 0, t162, t18 * qJD(2) + t166; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t184, t189, t31, t184, -t41 * qJD(2) - t141, 0, t61 * qJD(3) + t126 - t192, -t127 - t190, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t99, -t136, 0, t99, 0, 0, t132, t131, 0, 0, t125, -t151, 0, -t125, 0, 0, -t90 * qJD(3) - t150, t88 * qJD(3) - t149, 0, t52 * qJD(3) - t159, t128, -t189, 0, -t184, -t147, 0, -qJD(3) * t181 - t158, -t157 + t183, 0, -t17 * qJD(3) - t167; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t75, t74, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t138, t139, 0, t152, 0, 0, 0, 0, 0, 0, -t146, t182, 0, -t160; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t148, 0, -t115 - t75, t114 + t74, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t137, -t86, -t140, -t52 * qJD(2) - t155, 0, 0, 0, 0, 0, 0, t113, t31, -t162, t17 * qJD(2) - t166; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t138, -t139, 0, -t152, 0, 0, 0, 0, 0, 0, t146, -t182, 0, t160; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t153, -t182, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t184, -t189, 0, -t184, t61 * qJD(2), 0, -t41 * qJD(3) - t126, t127 + t183, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t148, 0, t115, -t114, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t153, t182, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t1;
