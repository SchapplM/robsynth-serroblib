% Calculate minimal parameter regressor of coriolis matrix for
% S5RRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x25]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 19:37
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RRPPR5_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR5_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR5_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR5_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 19:36:34
% EndTime: 2021-01-15 19:36:39
% DurationCPUTime: 1.50s
% Computational Cost: add. (1481->141), mult. (2908->194), div. (0->0), fcn. (3222->6), ass. (0->117)
t129 = qJD(2) - qJD(5);
t101 = sin(qJ(5));
t103 = cos(qJ(5));
t166 = cos(pkin(8));
t104 = cos(qJ(2));
t186 = -qJ(3) - pkin(6);
t93 = t186 * t104;
t116 = t166 * t93;
t100 = sin(pkin(8));
t102 = sin(qJ(2));
t92 = t186 * t102;
t183 = t100 * t92;
t197 = -t116 + t183;
t82 = t100 * t102 - t166 * t104;
t198 = t82 * pkin(7) + t197;
t56 = -t100 * t93 - t166 * t92;
t84 = t100 * t104 + t166 * t102;
t36 = t84 * pkin(7) - t56;
t218 = t129 * (t101 * t198 + t103 * t36);
t111 = t101 * t82 + t103 * t84;
t199 = t101 * t84 - t103 * t82;
t10 = t111 ^ 2 - t199 ^ 2;
t217 = t10 * qJD(1);
t216 = t129 * (-t101 * t36 + t103 * t198);
t149 = t199 * qJD(5);
t208 = t199 * qJD(2) - t149;
t207 = qJD(3) * t111;
t206 = qJD(3) * t199;
t81 = t84 ^ 2;
t200 = t82 ^ 2 + t81;
t205 = qJD(3) * t200;
t204 = t111 * qJD(1);
t203 = t199 * qJD(1);
t202 = t200 * qJD(1);
t146 = t111 * qJD(5);
t201 = qJD(2) * t111 - t146;
t190 = -t116 / 0.2e1;
t113 = -t197 * t82 + t56 * t84;
t195 = qJD(3) * t113;
t192 = t113 * qJD(1);
t191 = -t82 / 0.2e1;
t189 = -pkin(3) - pkin(4);
t188 = t84 * pkin(3);
t187 = t102 * pkin(2);
t185 = qJD(2) * pkin(2);
t98 = -t104 * pkin(2) - pkin(1);
t107 = t84 * qJ(4) - t98;
t42 = t82 * pkin(3) - t107;
t168 = t82 * qJ(4);
t112 = -t168 - t187;
t45 = -t112 + t188;
t5 = t42 * t45;
t172 = t5 * qJD(1);
t22 = t189 * t82 + t107;
t23 = t189 * t84 + t112;
t6 = -t111 * t22 + t199 * t23;
t171 = t6 * qJD(1);
t7 = t111 * t23 + t199 * t22;
t170 = t7 * qJD(1);
t11 = t98 * t187;
t163 = t11 * qJD(1);
t13 = t42 * t84 + t45 * t82;
t161 = t13 * qJD(1);
t14 = t42 * t82 - t45 * t84;
t160 = t14 * qJD(1);
t95 = t100 * pkin(2) + qJ(4);
t97 = -t166 * pkin(2) - pkin(3);
t99 = t187 / 0.2e1;
t18 = t99 + (pkin(3) / 0.2e1 - t97 / 0.2e1) * t84 + (qJ(4) / 0.2e1 + t95 / 0.2e1) * t82;
t157 = t18 * qJD(1);
t106 = t100 * t191 - t166 * t84 / 0.2e1;
t33 = (-t102 / 0.2e1 + t106) * pkin(2);
t154 = t33 * qJD(1);
t39 = t82 * t187 + t98 * t84;
t153 = t39 * qJD(1);
t40 = t84 * t187 - t98 * t82;
t152 = t40 * qJD(1);
t143 = t56 * qJD(2);
t142 = t81 * qJD(1);
t141 = t82 * qJD(1);
t72 = t82 * qJD(2);
t140 = t82 * qJD(3);
t139 = t84 * qJD(1);
t138 = t84 * qJD(4);
t94 = -t102 ^ 2 + t104 ^ 2;
t137 = t94 * qJD(1);
t136 = qJD(1) * t104;
t135 = t101 * qJD(2);
t134 = t101 * qJD(4);
t133 = t102 * qJD(2);
t132 = t103 * qJD(2);
t131 = t103 * qJD(4);
t130 = t104 * qJD(2);
t128 = pkin(1) * t102 * qJD(1);
t127 = pkin(1) * t136;
t126 = t22 * t203;
t125 = t22 * t204;
t124 = t199 * t204;
t123 = t111 * t203;
t122 = t199 * t139;
t121 = t111 * t139;
t120 = t82 * t139;
t117 = t102 * t136;
t110 = -pkin(4) + t97;
t54 = t190 + t116 / 0.2e1;
t109 = t54 * qJD(1) + t95 * qJD(2);
t91 = t129 * t103;
t90 = t129 * t101;
t76 = t84 * qJD(3);
t74 = t84 * qJD(2);
t62 = t101 * t110 + t103 * t95;
t61 = t101 * t95 - t103 * t110;
t51 = t197 * qJD(2);
t38 = 0.2e1 * t190 + t183;
t32 = t106 * pkin(2) + t99;
t19 = t95 * t191 + t97 * t84 / 0.2e1 + t99 + t168 / 0.2e1 + t188 / 0.2e1;
t1 = [0, 0, 0, t102 * t130, t94 * qJD(2), 0, 0, 0, -pkin(1) * t133, -pkin(1) * t130, t39 * qJD(2), t40 * qJD(2), t205, t11 * qJD(2) + t195, t13 * qJD(2) - t82 * t138, t205, t14 * qJD(2) + t81 * qJD(4), t5 * qJD(2) - t138 * t42 + t195, t208 * t111, t129 * t10, 0, 0, 0, t6 * qJD(2) + t138 * t199 + t146 * t22, t7 * qJD(2) + t111 * t138 - t149 * t22; 0, 0, 0, t117, t137, t130, -t133, 0, -pkin(6) * t130 - t128, pkin(6) * t133 - t127, -t51 + t153, t143 + t152, (-t100 * t84 + t166 * t82) * t185, t163 + (-t100 * t56 - t166 * t197) * t185 + t32 * qJD(3), -t51 + t161, (-t97 * t82 - t95 * t84) * qJD(2) - qJD(4) * t82, -t143 + t160, t172 + (t197 * t97 - t56 * t95) * qJD(2) + t19 * qJD(3) + t38 * qJD(4), t123, t217, -t208, -t201, 0, t171 - t216, t170 + t218; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t202, t32 * qJD(2) + t192, 0, t202, 0, t19 * qJD(2) + t192, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t120, -t72, t142, t38 * qJD(2) - t139 * t42, 0, 0, 0, 0, 0, t122, t121; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t124, -t217, t208, t201, 0, t125 + t216, -t126 - t218; 0, 0, 0, -t117, -t137, 0, 0, 0, t128, t127, -t76 - t153, t140 - t152, 0, t33 * qJD(3) - t163, -t76 - t161, 0, -t140 - t160, -t18 * qJD(3) + t54 * qJD(4) - t172, -t123, -t217, 0, 0, 0, -t171 - t207, -t170 + t206; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4), t95 * qJD(4), 0, 0, 0, 0, 0, t62 * qJD(5) + t134, -t61 * qJD(5) + t131; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t139, t141, 0, t154, -t139, 0, -t141, -t157, 0, 0, 0, 0, 0, -t204, t203; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2), t109, 0, 0, 0, 0, 0, t135, t132; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t129 * t62, -t129 * t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74, -t72, -t202, -t33 * qJD(2) - t192, t74, -t202, t72, t18 * qJD(2) - t138 - t192, 0, 0, 0, 0, 0, t201, -t208; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t139, -t141, 0, -t154, t139, 0, t141, t157, 0, 0, 0, 0, 0, t204, -t203; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t139, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t204, t203; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t120, 0, -t142, -t54 * qJD(2) + (qJD(1) * t42 + qJD(3)) * t84, 0, 0, 0, 0, 0, -t122, -t121; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(2), -t109, 0, 0, 0, 0, 0, -t90, -t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t139, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t90, t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t124, t217, 0, 0, 0, -t125 + t207, t126 - t206; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t62 * qJD(2) - t134, t61 * qJD(2) - t131; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t204, -t203; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t135, -t132; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t1;
