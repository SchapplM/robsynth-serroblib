% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S6RPPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPPRRR2_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR2_inertiaDJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR2_inertiaDJ_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR2_inertiaDJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:21:22
% EndTime: 2019-03-09 02:21:29
% DurationCPUTime: 2.68s
% Computational Cost: add. (4694->243), mult. (9950->420), div. (0->0), fcn. (9939->10), ass. (0->127)
t179 = sin(qJ(4));
t136 = t179 * qJD(3);
t181 = cos(qJ(4));
t139 = t181 * qJD(3);
t118 = sin(pkin(10)) * pkin(1) + qJ(3);
t115 = pkin(7) + t118;
t95 = cos(pkin(11));
t112 = t95 * t115;
t94 = sin(pkin(11));
t113 = t94 * t115;
t52 = t179 * t112 + t181 * t113;
t103 = t52 * qJD(4) + t94 * t136 - t95 * t139;
t78 = t179 * t94 - t181 * t95;
t79 = t179 * t95 + t181 * t94;
t82 = -cos(pkin(10)) * pkin(1) - pkin(2) - t95 * pkin(3);
t108 = t78 * pkin(4) - t79 * pkin(8) + t82;
t194 = -qJD(5) * t108 + t103;
t105 = t181 * t112 - t179 * t113;
t73 = t79 * qJD(4);
t183 = t73 * pkin(4);
t72 = t78 * qJD(4);
t184 = t72 * pkin(8);
t193 = qJD(5) * t105 - t183 - t184;
t178 = sin(qJ(6));
t96 = sin(qJ(5));
t141 = t178 * t96;
t133 = t79 * t141;
t180 = cos(qJ(6));
t144 = t180 * t96;
t189 = qJD(5) + qJD(6);
t137 = t180 * qJD(6);
t190 = t180 * qJD(5) + t137;
t97 = cos(qJ(5));
t19 = -t72 * t144 + (-t178 * t72 + t190 * t79) * t97 - t189 * t133;
t81 = t178 * t97 + t144;
t46 = t81 * t79;
t140 = qJD(6) * t178;
t58 = (t178 * qJD(5) + t140) * t96 - t190 * t97;
t159 = -t81 * t19 + t46 * t58;
t59 = t189 * t81;
t143 = t180 * t97;
t80 = t141 - t143;
t18 = t59 * t79 - t80 * t72;
t47 = t79 * t143 - t133;
t162 = t80 * t18 - t47 * t59;
t192 = t162 - t159;
t191 = t78 * t72 - t73 * t79;
t153 = qJD(5) * t97;
t114 = t79 * t153 - t96 * t72;
t154 = qJD(5) * t96;
t148 = t79 * t154;
t167 = t97 * t72;
t49 = t148 + t167;
t92 = t96 ^ 2;
t93 = t97 ^ 2;
t135 = qJD(5) * (t92 - t93);
t188 = 0.2e1 * (t94 ^ 2 + t95 ^ 2) * qJD(3);
t23 = -t96 * t105 + t97 * t108;
t24 = t97 * t105 + t96 * t108;
t123 = t23 * t96 - t24 * t97;
t8 = t193 * t96 + t194 * t97;
t9 = -t193 * t97 + t194 * t96;
t187 = t123 * qJD(5) + t8 * t96 - t9 * t97;
t61 = t78 * t73;
t57 = 0.2e1 * t61;
t185 = -pkin(9) - pkin(8);
t182 = t73 * pkin(5);
t177 = t46 * t19;
t176 = t47 * t18;
t35 = t105 * qJD(4) + t95 * t136 + t94 * t139;
t175 = t52 * t35;
t174 = t79 * t72;
t173 = t79 * t96;
t172 = t80 * t59;
t171 = t81 * t58;
t170 = t92 * t72;
t68 = t93 * t72;
t168 = t96 * t73;
t166 = t97 * t79;
t161 = -t18 * t81 - t47 * t58;
t160 = t19 * t80 + t46 * t59;
t157 = t73 * t166 - t78 * t167;
t152 = -0.2e1 * pkin(4) * qJD(5);
t151 = t96 * t167;
t150 = pkin(5) * t154;
t149 = t78 * t154;
t146 = t96 * t153;
t76 = t79 ^ 2;
t134 = t76 * t146;
t132 = pkin(5) * t137;
t131 = pkin(5) * t140;
t130 = t185 * t180;
t129 = t185 * t178;
t127 = pkin(4) * t72 - pkin(8) * t73;
t126 = pkin(4) * t79 + pkin(8) * t78;
t124 = t23 * t97 + t24 * t96;
t122 = t35 * t78 + t52 * t73;
t121 = t58 * t78 - t73 * t81;
t120 = t59 * t78 + t73 * t80;
t117 = t96 * t130;
t116 = t96 * t129;
t50 = t78 * t153 + t168;
t101 = t78 * pkin(5) - pkin(9) * t166 + t23;
t100 = t180 * t101;
t102 = -t114 * pkin(9) - t8;
t21 = -pkin(9) * t173 + t24;
t98 = t49 * pkin(9) + t182 + t9;
t1 = -qJD(6) * t100 - t180 * t102 + t21 * t140 - t178 * t98;
t83 = t185 * t97;
t64 = -t180 * t83 + t116;
t107 = -t124 * qJD(5) - t8 * t97 - t9 * t96;
t99 = t178 * t101;
t2 = -qJD(6) * t99 - t178 * t102 - t21 * t137 + t180 * t98;
t89 = -pkin(5) * t97 - pkin(4);
t63 = t178 * t83 + t117;
t54 = t79 * t68;
t53 = t79 * t170;
t48 = -t97 * t73 + t149;
t44 = t68 + t170;
t38 = pkin(5) * t173 + t52;
t34 = t79 * t135 + t151;
t33 = -t64 * qJD(6) + (t97 * t130 - t116) * qJD(5);
t32 = -t189 * t117 - t129 * t153 - t83 * t140;
t22 = t114 * pkin(5) + t35;
t7 = t180 * t21 + t99;
t6 = -t178 * t21 + t100;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t188, t118 * t188, -0.2e1 * t174, 0.2e1 * t191, 0, t57, 0, 0, 0.2e1 * t82 * t73, -0.2e1 * t82 * t72, 0.2e1 * t103 * t78 - 0.2e1 * t105 * t73 + 0.2e1 * t35 * t79 - 0.2e1 * t52 * t72, -0.2e1 * t105 * t103 + 0.2e1 * t175, -0.2e1 * t54 - 0.2e1 * t134, 0.2e1 * t76 * t135 + 0.4e1 * t79 * t151, -0.2e1 * t78 * t148 + 0.2e1 * t157, -0.2e1 * t53 + 0.2e1 * t134, -0.2e1 * t114 * t78 - 0.2e1 * t79 * t168, t57, 0.2e1 * t114 * t52 + 0.2e1 * t35 * t173 + 0.2e1 * t23 * t73 + 0.2e1 * t9 * t78, 0.2e1 * t35 * t166 - 0.2e1 * t24 * t73 - 0.2e1 * t49 * t52 + 0.2e1 * t8 * t78, 0.2e1 * t124 * t72 + 0.2e1 * t187 * t79, 0.2e1 * t23 * t9 - 0.2e1 * t24 * t8 + 0.2e1 * t175, -0.2e1 * t176, 0.2e1 * t46 * t18 - 0.2e1 * t47 * t19, -0.2e1 * t18 * t78 + 0.2e1 * t47 * t73, 0.2e1 * t177, -0.2e1 * t78 * t19 - 0.2e1 * t73 * t46, t57, 0.2e1 * t19 * t38 + 0.2e1 * t2 * t78 + 0.2e1 * t22 * t46 + 0.2e1 * t6 * t73, 0.2e1 * t1 * t78 - 0.2e1 * t18 * t38 + 0.2e1 * t22 * t47 - 0.2e1 * t7 * t73, 0.2e1 * t1 * t46 + 0.2e1 * t18 * t6 - 0.2e1 * t19 * t7 - 0.2e1 * t2 * t47, -0.2e1 * t1 * t7 + 0.2e1 * t2 * t6 + 0.2e1 * t22 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t103 * t79 - t105 * t72 + t122, 0, 0, 0, 0, 0, 0, 0, t191 * t97 + t157, 0, t107 * t79 + t123 * t72 + t122, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1 * t47 - t18 * t7 - t19 * t6 - t2 * t46 + t22 * t78 + t38 * t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t174 + 0.2e1 * t61, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t54 - 0.2e1 * t53 + 0.2e1 * t61, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t176 + 0.2e1 * t61 + 0.2e1 * t177; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73, -t72, 0, 0, 0, 0, 0, 0, 0, 0, -t48, -t50, t44, -t187, 0, 0, 0, 0, 0, 0, -t120, t121, -t192, -t1 * t81 - t2 * t80 - t58 * t7 - t59 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t160 + t161; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t171 + 0.2e1 * t172; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t72, 0, -t73, 0, -t35, t103, 0, 0, -t34, -0.4e1 * t79 * t146 + t170 - t68, t50, t34, -t48, 0, -t35 * t97 + t127 * t96 + (-t126 * t97 + t52 * t96) * qJD(5), t35 * t96 + t127 * t97 + (t126 * t96 + t52 * t97) * qJD(5), t107, -t35 * pkin(4) + t107 * pkin(8), t161, t159 + t162, -t121, t160, -t120, 0, t46 * t150 + t19 * t89 + t22 * t80 + t33 * t78 + t38 * t59 + t63 * t73, t47 * t150 - t18 * t89 + t22 * t81 + t32 * t78 - t38 * t58 - t64 * t73, t1 * t80 + t18 * t63 - t19 * t64 - t2 * t81 + t32 * t46 - t33 * t47 + t58 * t6 - t59 * t7, -t1 * t64 + t38 * t150 + t2 * t63 + t22 * t89 - t32 * t7 + t33 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t73, t72, 0, 0, 0, 0, 0, 0, 0, 0, t48, t50, -t44, -t183 + (-t92 - t93) * t184, 0, 0, 0, 0, 0, 0, t120, -t121, t192, pkin(5) * t149 - t18 * t64 - t19 * t63 - t32 * t47 - t33 * t46 + t73 * t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32 * t81 - t33 * t80 - t58 * t64 - t59 * t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t146, -0.2e1 * t135, 0, -0.2e1 * t146, 0, 0, t96 * t152, t97 * t152, 0, 0, -0.2e1 * t171, 0.2e1 * t80 * t58 - 0.2e1 * t81 * t59, 0, 0.2e1 * t172, 0, 0, 0.2e1 * t80 * t150 + 0.2e1 * t59 * t89, 0.2e1 * t81 * t150 - 0.2e1 * t58 * t89, 0.2e1 * t32 * t80 - 0.2e1 * t33 * t81 + 0.2e1 * t58 * t63 - 0.2e1 * t59 * t64, 0.2e1 * t89 * t150 - 0.2e1 * t32 * t64 + 0.2e1 * t33 * t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49, 0, -t114, t73, t9, t8, 0, 0, 0, 0, -t18, 0, -t19, t73, -t78 * t131 + t180 * t182 + t2 (-t78 * t137 - t178 * t73) * pkin(5) + t1 (t180 * t18 - t178 * t19 + (t178 * t47 - t180 * t46) * qJD(6)) * pkin(5) (t180 * t2 - t178 * t1 + (-t178 * t6 + t180 * t7) * qJD(6)) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t114, t49, 0, 0, 0, 0, 0, 0, 0, 0, -t19, t18, 0 (-t180 * t19 - t178 * t18 + (t178 * t46 + t180 * t47) * qJD(6)) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t154, -t153, 0, 0, 0, 0, 0, 0, 0, 0, -t59, t58, 0 (-t180 * t59 - t178 * t58 + (t178 * t80 + t180 * t81) * qJD(6)) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t153, 0, -t154, 0, -pkin(8) * t153, pkin(8) * t154, 0, 0, 0, 0, -t58, 0, -t59, 0, t33, t32 (t180 * t58 - t178 * t59 + (t178 * t81 - t180 * t80) * qJD(6)) * pkin(5) (t180 * t33 - t178 * t32 + (-t178 * t63 + t180 * t64) * qJD(6)) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t131, -0.2e1 * t132, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, 0, -t19, t73, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, t18, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t59, t58, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t58, 0, -t59, 0, t33, t32, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t131, -t132, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t3;
