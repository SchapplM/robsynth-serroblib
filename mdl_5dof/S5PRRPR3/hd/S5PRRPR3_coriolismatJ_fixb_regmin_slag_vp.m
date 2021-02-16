% Calculate minimal parameter regressor of coriolis matrix for
% S5PRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x22]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 15:42
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5PRRPR3_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR3_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR3_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR3_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 15:42:13
% EndTime: 2021-01-15 15:42:16
% DurationCPUTime: 0.93s
% Computational Cost: add. (1243->108), mult. (2594->156), div. (0->0), fcn. (2894->6), ass. (0->98)
t160 = cos(qJ(5));
t151 = -qJ(4) - pkin(6);
t96 = sin(qJ(3));
t86 = t151 * t96;
t97 = cos(qJ(3));
t87 = t151 * t97;
t93 = sin(pkin(9));
t94 = cos(pkin(9));
t58 = -t93 * t86 + t94 * t87;
t78 = t93 * t96 - t94 * t97;
t40 = -t78 * pkin(7) - t58;
t172 = t160 * t40;
t107 = -t172 / 0.2e1;
t108 = t94 * t86 + t93 * t87;
t80 = t93 * t97 + t94 * t96;
t100 = -t80 * pkin(7) + t108;
t163 = t160 * t100;
t95 = sin(qJ(5));
t173 = t95 * t40;
t175 = -t163 + t173;
t165 = t95 * t100;
t174 = -t172 - t165;
t152 = t95 * t78;
t75 = t160 * t80;
t164 = t75 - t152;
t109 = t160 * t78 + t95 * t80;
t167 = t109 * qJD(2);
t171 = t164 * t167;
t166 = t109 ^ 2 - t164 ^ 2;
t170 = t166 * qJD(2);
t106 = -t163 / 0.2e1;
t168 = qJD(4) * t109;
t133 = t109 * qJD(5);
t102 = qJD(3) * t109 + t133;
t105 = t75 / 0.2e1;
t162 = pkin(3) * t93;
t161 = t96 * pkin(3);
t159 = t78 * t93;
t158 = t80 * t94;
t150 = qJD(3) * pkin(3);
t91 = -t97 * pkin(3) - pkin(2);
t62 = t78 * pkin(4) + t91;
t146 = qJD(2) * t62;
t145 = qJD(2) * t97;
t63 = t80 * pkin(4) + t161;
t10 = t109 * t63 + t164 * t62;
t144 = t10 * qJD(2);
t11 = -t109 * t62 + t164 * t63;
t143 = t11 * qJD(2);
t18 = -t108 * t80 + t78 * t58;
t141 = t18 * qJD(2);
t22 = 0.2e1 * t105 - t152;
t139 = t22 * qJD(2);
t99 = -t159 / 0.2e1 - t158 / 0.2e1;
t32 = (-t96 / 0.2e1 + t99) * pkin(3);
t138 = t32 * qJD(2);
t41 = t78 * t161 + t91 * t80;
t137 = t41 * qJD(2);
t42 = t80 * t161 - t91 * t78;
t136 = t42 * qJD(2);
t47 = t105 - t75 / 0.2e1;
t135 = t47 * qJD(2);
t45 = t47 * qJD(5);
t134 = t164 * qJD(2);
t130 = t164 * qJD(5);
t53 = t78 ^ 2 + t80 ^ 2;
t129 = t53 * qJD(2);
t128 = t78 * qJD(2);
t127 = t78 * qJD(3);
t126 = t80 * qJD(2);
t125 = t80 * qJD(3);
t88 = -t96 ^ 2 + t97 ^ 2;
t124 = t88 * qJD(2);
t123 = t96 * qJD(3);
t122 = t97 * qJD(3);
t121 = pkin(2) * t96 * qJD(2);
t120 = pkin(2) * t145;
t117 = t109 * t146;
t116 = t164 * t146;
t115 = t96 * t145;
t9 = t91 * t161;
t104 = t9 * qJD(2);
t3 = t106 + t163 / 0.2e1;
t90 = t94 * pkin(3) + pkin(4);
t72 = -t160 * t90 + t95 * t162;
t103 = -t3 * qJD(2) - t72 * qJD(3);
t101 = qJD(3) * t164 + t22 * qJD(5);
t2 = t107 + t172 / 0.2e1;
t73 = t160 * t162 + t95 * t90;
t98 = t47 * qJD(1) + t2 * qJD(2) + t73 * qJD(3);
t65 = t73 * qJD(5);
t64 = t72 * qJD(5);
t44 = t47 * qJD(3);
t31 = t161 / 0.2e1 + t99 * pkin(3);
t13 = -t22 * qJD(3) - t130;
t5 = 0.2e1 * t107 - t165;
t4 = t173 + 0.2e1 * t106;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t123, -t122, -t125, t127, 0, (-t158 - t159) * t150, 0, 0, 0, 0, 0, -t101, t102; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, t102; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, t96 * t122, t88 * qJD(3), 0, 0, 0, -pkin(2) * t123, -pkin(2) * t122, t41 * qJD(3), t42 * qJD(3), t53 * qJD(4), t9 * qJD(3) + t18 * qJD(4), -t102 * t164, (qJD(3) + qJD(5)) * t166, 0, 0, 0, t10 * qJD(3) + t62 * t130, t11 * qJD(3) - t62 * t133; 0, 0, 0, 0, t115, t124, t122, -t123, 0, -pkin(6) * t122 - t121, pkin(6) * t123 - t120, t58 * qJD(3) + t137, -qJD(3) * t108 + t136, (t78 * t94 - t80 * t93) * t150, t31 * qJD(4) + (t108 * t93 + t58 * t94) * t150 + t104, -t171, t170, -t102, -t101, 0, t174 * qJD(3) + t5 * qJD(5) + t144, t175 * qJD(3) + t4 * qJD(5) + t143; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t129, t31 * qJD(3) + t141, 0, 0, 0, 0, 0, t45, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t171, t170, -t102, t13, 0, t5 * qJD(3) + t47 * qJD(4) + t174 * qJD(5) + t116, t4 * qJD(3) + t175 * qJD(5) - t117; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t45, 0; 0, 0, 0, 0, -t115, -t124, 0, 0, 0, t121, t120, -t80 * qJD(4) - t137, t78 * qJD(4) - t136, 0, t32 * qJD(4) - t104, t171, -t170, 0, -t45, 0, -qJD(4) * t164 - t2 * qJD(5) - t144, t3 * qJD(5) - t143 + t168; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t65, t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t126, t128, 0, t138, 0, 0, 0, 0, 0, -t134, t167; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t135, 0, -t65 - t98, -t103 + t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t125, -t127, -t129, -t32 * qJD(3) - t141, 0, 0, 0, 0, 0, t101, -t102; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t126, -t128, 0, -t138, 0, 0, 0, 0, 0, t134, -t167; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t139, -t167; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t171, -t170, 0, t44, 0, t2 * qJD(3) - t22 * qJD(4) - t116, -t3 * qJD(3) + t117 + t168; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t135, 0, t98, t103; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t139, t167; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t1;
