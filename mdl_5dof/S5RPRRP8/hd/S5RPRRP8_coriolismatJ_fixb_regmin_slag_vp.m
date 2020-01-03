% Calculate minimal parameter regressor of coriolis matrix for
% S5RPRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x20]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RPRRP8_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP8_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP8_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP8_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:47:32
% EndTime: 2019-12-31 18:47:36
% DurationCPUTime: 1.46s
% Computational Cost: add. (1176->203), mult. (1937->227), div. (0->0), fcn. (1422->4), ass. (0->142)
t107 = sin(qJ(4));
t108 = sin(qJ(3));
t110 = cos(qJ(3));
t187 = -pkin(1) - pkin(2);
t75 = t108 * qJ(2) - t110 * t187;
t73 = pkin(3) + t75;
t194 = pkin(3) / 0.2e1 + t73 / 0.2e1;
t197 = t107 * t194;
t109 = cos(qJ(4));
t196 = t109 * t194;
t163 = t109 * qJ(5);
t184 = t107 * pkin(4);
t82 = t163 - t184;
t168 = t82 * t109;
t81 = t109 * pkin(4) + t107 * qJ(5);
t77 = -pkin(3) - t81;
t38 = t75 - t77;
t174 = t38 * t107;
t71 = t77 * t107;
t182 = t174 / 0.2e1 - t71 / 0.2e1;
t195 = t182 - t168;
t105 = t107 ^ 2;
t106 = t109 ^ 2;
t154 = t105 + t106;
t139 = qJD(1) - qJD(3);
t89 = t106 - t105;
t192 = t139 * t89;
t125 = t106 / 0.2e1 + t105 / 0.2e1;
t124 = t125 * pkin(7);
t134 = -t77 / 0.2e1 + t38 / 0.2e1;
t76 = t110 * qJ(2) + t108 * t187;
t74 = -pkin(7) + t76;
t1 = (t125 * t75 - t134) * t108 + (t76 / 0.2e1 - t125 * t74 + t124) * t110;
t32 = (-0.1e1 + t154) * t110 * t108;
t27 = t32 * qJD(2);
t191 = -t1 * qJD(1) + t27;
t142 = t109 * qJD(5);
t111 = -t81 * qJD(4) + t142;
t72 = t154 * t110;
t159 = t72 * qJD(1);
t189 = -t72 * qJD(3) + t159;
t188 = -t82 / 0.2e1;
t183 = t82 * t77;
t150 = qJD(3) * t107;
t54 = t76 * t150;
t152 = qJD(1) * t107;
t55 = t76 * t152;
t180 = t54 - t55;
t148 = qJD(3) * t109;
t56 = t76 * t148;
t151 = qJD(1) * t109;
t57 = t76 * t151;
t179 = t57 - t56;
t91 = t107 * t142;
t143 = t108 * qJD(2);
t94 = t109 * t143;
t178 = t91 + t94;
t176 = t107 * t75;
t175 = t109 * t75;
t173 = t38 * t108;
t172 = t38 * t109;
t24 = t154 * t75;
t5 = -t74 * t24 + t38 * t76;
t171 = t5 * qJD(1);
t170 = t77 * t109;
t169 = t82 * t107;
t167 = t82 * t110;
t101 = t105 * qJD(5);
t90 = t107 * t143;
t166 = t101 + t90;
t165 = qJD(1) * t38;
t10 = t74 * t72 + t173;
t164 = t10 * qJD(1);
t13 = t169 + t172;
t162 = t13 * qJD(1);
t14 = t168 - t174;
t161 = t14 * qJD(1);
t118 = t163 / 0.2e1 - t184 / 0.2e1;
t21 = (t82 / 0.2e1 + t118) * t110;
t160 = t21 * qJD(1);
t158 = t72 * qJD(2);
t155 = t89 * qJD(4);
t153 = qJ(2) * qJD(1);
t149 = qJD(3) * t108;
t147 = qJD(4) * qJ(5);
t146 = t107 * qJD(4);
t145 = t107 * qJD(5);
t144 = t108 * qJD(1);
t103 = t109 * qJD(4);
t141 = t110 * qJD(1);
t140 = t110 * qJD(2);
t137 = pkin(7) * t146;
t136 = pkin(7) * t103;
t135 = t82 * t165;
t133 = t74 * t146;
t132 = t74 * t103;
t129 = t107 * t144;
t128 = t109 * t144;
t127 = t107 * t141;
t126 = t109 * t141;
t80 = t139 * t110;
t25 = -t169 + t170;
t51 = -t175 / 0.2e1;
t7 = t134 * t109 + t169 + t51;
t123 = t7 * qJD(1) + t25 * qJD(3);
t26 = -t71 - t168;
t50 = t176 / 0.2e1;
t6 = t50 - t195;
t122 = t6 * qJD(1) + t26 * qJD(3);
t121 = qJD(1) * t73 - t140;
t120 = qJD(4) * t82 + t145;
t15 = t50 - t197;
t117 = pkin(3) * t150 + t15 * qJD(1);
t52 = t175 / 0.2e1;
t16 = t52 - t196;
t116 = pkin(3) * t148 + t16 * qJD(1);
t49 = -t176 / 0.2e1;
t11 = t49 + t182;
t115 = t11 * qJD(1) + t77 * t150;
t114 = t118 * t75;
t113 = t118 * t110;
t20 = (t188 + t118) * t110;
t3 = t134 * t82 - t114;
t112 = t3 * qJD(1) + t20 * qJD(2) + qJD(3) * t183;
t93 = t107 * t103;
t79 = t139 * t108;
t78 = t139 * t105;
t58 = (-t148 + t151) * t107;
t36 = t108 * t146 + t109 * t80;
t35 = -t108 * t103 + t107 * t80;
t34 = -t108 * t148 - t110 * t146 + t128;
t33 = t110 * t103 - t107 * t149 + t129;
t23 = -t167 / 0.2e1 + t113;
t22 = t167 / 0.2e1 + t113;
t18 = t52 + t196;
t17 = t50 + t197;
t12 = t49 - t182;
t9 = t50 + t195;
t8 = -t169 + t170 / 0.2e1 - t172 / 0.2e1 + t51;
t4 = t183 / 0.2e1 + t38 * t188 - t114;
t2 = t173 / 0.2e1 + t108 * t77 / 0.2e1 + (-t76 / 0.2e1 + t124) * t110 + t154 * (-t75 * t108 / 0.2e1 + t74 * t110 / 0.2e1);
t19 = [0, 0, 0, 0, qJD(2), qJ(2) * qJD(2), 0, t76 * qJD(3) + t143, -t75 * qJD(3) + t140, t93, t155, 0, 0, 0, -t73 * t146 + t56 + t94, -t73 * t103 - t54 - t90, t14 * qJD(4) + t178 + t56, t24 * qJD(3) - t158, t13 * qJD(4) + t166 + t54, t10 * qJD(2) + t5 * qJD(3) + t120 * t38; 0, 0, 0, 0, qJD(1), t153, 0, t144, t141, 0, 0, 0, 0, 0, t128, -t129, t128, -t159, t129, t2 * qJD(3) + t23 * qJD(4) + t164 + t27; 0, 0, 0, 0, 0, 0, 0, t139 * t76, -t139 * t75, -t93, -t155, 0, 0, 0, t17 * qJD(4) + t179, t18 * qJD(4) + t180, t9 * qJD(4) + t179 - t91, t139 * t24, t8 * qJD(4) - t101 - t180, t171 + t2 * qJD(2) + (-pkin(7) * t24 + t76 * t77) * qJD(3) + t4 * qJD(4) + t12 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, t192, -t103, t146, 0, t17 * qJD(3) - t73 * t152 - t132, t18 * qJD(3) - t73 * t151 + t133, t9 * qJD(3) - t132 + t161, -t111, t8 * qJD(3) - t133 + t162, t23 * qJD(2) + t4 * qJD(3) + t111 * t74 + t135; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, -t103, t78, t12 * qJD(3) + t38 * t152 + t132; 0, 0, 0, 0, -qJD(1), -t153, 0, -t79, -t80, 0, 0, 0, 0, 0, -t34, t33, -t34, t189, -t33, -t1 * qJD(3) - t21 * qJD(4) - t110 * t145 - t164; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32 * qJD(3); 0, 0, 0, 0, 0, 0, 0, t79, t80, 0, 0, 0, 0, 0, t34, -t33, t34, -t189, t33, t77 * t149 + t22 * qJD(4) + (t154 * qJD(3) * pkin(7) + t145) * t110 + t191; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, t36, t35, 0, -t36, t22 * qJD(3) + t111 * t108 - t160; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35; 0, 0, 0, 0, 0, 0, 0, -t76 * qJD(1) - t143, t75 * qJD(1) - t140, -t93, -t155, 0, 0, 0, -t15 * qJD(4) - t57 - t94, -t16 * qJD(4) + t55 + t90, -t6 * qJD(4) - t178 - t57, -t24 * qJD(1) + t158, -t7 * qJD(4) - t166 - t55, t1 * qJD(2) - t3 * qJD(4) - t11 * qJD(5) - t171; 0, 0, 0, 0, 0, 0, 0, -t144, -t141, 0, 0, 0, 0, 0, -t128, t129, -t128, t159, -t129, -t20 * qJD(4) - t191; 0, 0, 0, 0, 0, 0, 0, 0, 0, t93, t155, 0, 0, 0, -pkin(3) * t146, -pkin(3) * t103, -t26 * qJD(4) + t91, 0, -t25 * qJD(4) + t101, -t120 * t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t58, -t192, t103, -t146, 0, -t117 - t136, -t116 + t137, -t122 - t136, t111, -t123 - t137, pkin(7) * t111 - t112; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t58, t103, -t78, -t115 + t136; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t58, -t192, 0, 0, 0, t15 * qJD(3) + t121 * t107, t16 * qJD(3) + t121 * t109, t6 * qJD(3) - t107 * t140 - t161, 0, t7 * qJD(3) + t109 * t140 - t162, t21 * qJD(2) + t3 * qJD(3) - t135; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t127, -t126, -t127, 0, t126, t20 * qJD(3) + t160; 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, t192, 0, 0, 0, t117, t116, t122, 0, t123, t112; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(5), qJ(5) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4), t147; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t58, 0, -t78, t11 * qJD(3) + (t140 - t165) * t107; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t127; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, 0, t78, t115; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(4), -t147; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t19;
