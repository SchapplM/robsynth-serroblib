% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRRRP4_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP4_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP4_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP4_inertiaDJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:51:13
% EndTime: 2019-12-31 21:51:17
% DurationCPUTime: 1.17s
% Computational Cost: add. (1486->150), mult. (3443->214), div. (0->0), fcn. (2771->6), ass. (0->103)
t151 = sin(qJ(2));
t128 = t151 * pkin(1);
t115 = t128 + pkin(7);
t111 = -pkin(8) - t115;
t95 = sin(qJ(3));
t106 = t111 * t95;
t152 = cos(qJ(4));
t100 = t152 * t106;
t94 = sin(qJ(4));
t135 = qJD(4) * t94;
t96 = cos(qJ(3));
t91 = t96 * pkin(8);
t64 = t96 * t115 + t91;
t104 = t111 * qJD(3);
t153 = cos(qJ(2));
t129 = t153 * pkin(1);
t117 = qJD(2) * t129;
t112 = t96 * t117;
t98 = t95 * t104 + t112;
t113 = t95 * t117;
t99 = -t96 * t104 + t113;
t12 = -qJD(4) * t100 + t64 * t135 - t152 * t98 + t94 * t99;
t103 = t94 * t106;
t120 = t152 * qJD(4);
t13 = qJD(4) * t103 + t64 * t120 + t152 * t99 + t94 * t98;
t44 = t94 * t64 - t100;
t45 = t152 * t64 + t103;
t126 = t152 * t96;
t146 = t94 * t95;
t161 = qJD(3) + qJD(4);
t47 = -qJD(3) * t126 - t96 * t120 + t161 * t146;
t145 = t94 * t96;
t68 = t152 * t95 + t145;
t48 = t161 * t68;
t67 = -t126 + t146;
t159 = t12 * t67 + t13 * t68 - t44 * t47 - t45 * t48;
t166 = 0.2e1 * t159;
t158 = pkin(8) + pkin(7);
t116 = t158 * t152;
t109 = qJD(3) * t116;
t110 = t95 * t116;
t124 = t158 * qJD(3);
t74 = t96 * pkin(7) + t91;
t31 = qJD(4) * t110 + t95 * t109 + t124 * t145 + t74 * t135;
t32 = t74 * t120 + t96 * t109 + (-qJD(4) * t158 - t124) * t146;
t52 = t94 * t74 + t110;
t53 = -t158 * t146 + t152 * t74;
t160 = t31 * t67 + t32 * t68 - t52 * t47 - t53 * t48;
t165 = 0.2e1 * t160;
t164 = t160 + t159;
t163 = 0.2e1 * t47 * t67 - 0.2e1 * t68 * t48;
t92 = t95 ^ 2;
t93 = t96 ^ 2;
t162 = t92 + t93;
t97 = 2 * qJD(5);
t157 = t96 * pkin(3);
t134 = t95 * qJD(3);
t87 = pkin(3) * t134;
t19 = t48 * pkin(4) + t47 * qJ(5) - t68 * qJD(5) + t87;
t89 = qJD(2) * t128;
t15 = t19 + t89;
t86 = -pkin(2) - t157;
t46 = t67 * pkin(4) - t68 * qJ(5) + t86;
t41 = -t129 + t46;
t155 = t15 * t67 + t41 * t48;
t154 = -t15 * t68 + t41 * t47;
t144 = t19 * t67 + t46 * t48;
t143 = -t19 * t68 + t46 * t47;
t70 = t89 + t87;
t85 = -t129 - pkin(2);
t73 = t85 - t157;
t141 = t73 * t48 + t70 * t67;
t140 = -t73 * t47 + t70 * t68;
t139 = t86 * t48 + t67 * t87;
t138 = -t86 * t47 + t68 * t87;
t90 = t96 * qJD(3);
t137 = t85 * t90 + t95 * t89;
t136 = pkin(1) * qJD(2);
t133 = 0.2e1 * t67 * t48;
t132 = pkin(2) * t90;
t131 = pkin(2) * t134;
t130 = pkin(3) * t135;
t127 = t95 * t90;
t125 = -t45 * t12 + t44 * t13;
t121 = -t53 * t31 + t52 * t32;
t108 = t115 * qJD(3);
t107 = -t12 * t53 + t13 * t52 - t45 * t31 + t44 * t32;
t105 = t85 * t134 - t96 * t89;
t102 = t162 * t153;
t88 = pkin(3) * t120;
t84 = -t152 * pkin(3) - pkin(4);
t82 = t94 * pkin(3) + qJ(5);
t81 = -0.2e1 * t130;
t78 = t88 + qJD(5);
t77 = -0.2e1 * t127;
t76 = 0.2e1 * t127;
t66 = 0.2e1 * (-t92 + t93) * qJD(3);
t60 = t102 * t136;
t38 = -0.2e1 * t68 * t47;
t20 = pkin(4) * t47 - t48 * qJ(5) - t67 * qJD(5);
t7 = (t152 * t47 - t48 * t94 + (-t152 * t67 + t68 * t94) * qJD(4)) * pkin(3);
t3 = t68 * t130 - t84 * t47 - t82 * t48 - t78 * t67;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t89, -0.2e1 * t117, 0, 0, t76, t66, 0, t77, 0, 0, 0.2e1 * t105, 0.2e1 * t137, 0.2e1 * t60, 0.2e1 * (t162 * t115 * t129 + t85 * t128) * qJD(2), t38, t163, 0, t133, 0, 0, 0.2e1 * t141, 0.2e1 * t140, t166, 0.2e1 * t73 * t70 + 0.2e1 * t125, t38, 0, -t163, 0, 0, t133, 0.2e1 * t155, t166, 0.2e1 * t154, 0.2e1 * t41 * t15 + 0.2e1 * t125; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t89, -t117, 0, 0, t76, t66, 0, t77, 0, 0, t105 - t131, -t132 + t137, t60, (-t151 * pkin(2) + t102 * pkin(7)) * t136, t38, t163, 0, t133, 0, 0, t139 + t141, t138 + t140, t164, t70 * t86 + t73 * t87 + t107, t38, 0, -t163, 0, 0, t133, t144 + t155, t164, t143 + t154, t15 * t46 + t41 * t19 + t107; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t76, t66, 0, t77, 0, 0, -0.2e1 * t131, -0.2e1 * t132, 0, 0, t38, t163, 0, t133, 0, 0, 0.2e1 * t139, 0.2e1 * t138, t165, 0.2e1 * t86 * t87 + 0.2e1 * t121, t38, 0, -t163, 0, 0, t133, 0.2e1 * t144, t165, 0.2e1 * t143, 0.2e1 * t46 * t19 + 0.2e1 * t121; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t90, 0, -t134, 0, -t96 * t108 - t113, t95 * t108 - t112, 0, 0, 0, 0, -t47, 0, -t48, 0, -t13, t12, t7, (-t152 * t13 - t12 * t94 + (t152 * t45 + t44 * t94) * qJD(4)) * pkin(3), 0, -t47, 0, 0, t48, 0, -t13, t3, -t12, -t12 * t82 + t13 * t84 + t44 * t130 + t45 * t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t90, 0, -t134, 0, -pkin(7) * t90, pkin(7) * t134, 0, 0, 0, 0, -t47, 0, -t48, 0, -t32, t31, t7, (-t152 * t32 - t31 * t94 + (t152 * t53 + t52 * t94) * qJD(4)) * pkin(3), 0, -t47, 0, 0, t48, 0, -t32, t3, -t31, t52 * t130 - t31 * t82 + t32 * t84 + t53 * t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81, -0.2e1 * t88, 0, 0, 0, 0, 0, 0, 0, 0, t81, 0, 0.2e1 * t78, 0.2e1 * t84 * t130 + 0.2e1 * t82 * t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t47, 0, -t48, 0, -t13, t12, 0, 0, 0, -t47, 0, 0, t48, 0, -t13, t20, -t12, -t13 * pkin(4) - t12 * qJ(5) + t45 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t47, 0, -t48, 0, -t32, t31, 0, 0, 0, -t47, 0, 0, t48, 0, -t32, t20, -t31, -t32 * pkin(4) - t31 * qJ(5) + t53 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t130, -t88, 0, 0, 0, 0, 0, 0, 0, 0, -t130, 0, t97 + t88, -pkin(4) * t130 + t78 * qJ(5) + t82 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t97, qJ(5) * t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t47, 0, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t47, 0, t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t130; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t1;
