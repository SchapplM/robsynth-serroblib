% Calculate minimal parameter regressor of coriolis matrix for
% S5RPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x22]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 11:45
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RPRPR4_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR4_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR4_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR4_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 11:44:18
% EndTime: 2021-01-15 11:44:21
% DurationCPUTime: 0.98s
% Computational Cost: add. (1500->110), mult. (2851->158), div. (0->0), fcn. (3151->8), ass. (0->100)
t163 = cos(qJ(5));
t91 = sin(pkin(8)) * pkin(1) + pkin(6);
t150 = qJ(4) + t91;
t99 = sin(qJ(3));
t82 = t150 * t99;
t100 = cos(qJ(3));
t83 = t150 * t100;
t95 = sin(pkin(9));
t96 = cos(pkin(9));
t43 = t95 * t82 - t96 * t83;
t84 = -t96 * t100 + t95 * t99;
t36 = -t84 * pkin(7) - t43;
t175 = t163 * t36;
t110 = -t175 / 0.2e1;
t112 = -t96 * t82 - t95 * t83;
t86 = t95 * t100 + t96 * t99;
t103 = -t86 * pkin(7) + t112;
t166 = t163 * t103;
t98 = sin(qJ(5));
t176 = t98 * t36;
t178 = -t166 + t176;
t168 = t98 * t103;
t177 = -t175 - t168;
t155 = t98 * t84;
t79 = t163 * t86;
t167 = t79 - t155;
t111 = t163 * t84 + t98 * t86;
t170 = t111 * qJD(1);
t174 = t167 * t170;
t169 = t111 ^ 2 - t167 ^ 2;
t173 = t169 * qJD(1);
t109 = -t166 / 0.2e1;
t171 = qJD(4) * t111;
t137 = t111 * qJD(5);
t105 = qJD(3) * t111 + t137;
t108 = t79 / 0.2e1;
t165 = pkin(3) * t95;
t164 = t99 * pkin(3);
t162 = t84 * t95;
t161 = t86 * t96;
t154 = qJD(3) * pkin(3);
t93 = -cos(pkin(8)) * pkin(1) - pkin(2);
t88 = -t100 * pkin(3) + t93;
t62 = t84 * pkin(4) + t88;
t149 = qJD(1) * t62;
t63 = t86 * pkin(4) + t164;
t10 = t111 * t63 + t167 * t62;
t148 = t10 * qJD(1);
t11 = -t111 * t62 + t167 * t63;
t147 = t11 * qJD(1);
t13 = -t112 * t86 + t43 * t84;
t145 = t13 * qJD(1);
t38 = 0.2e1 * t108 - t155;
t143 = t38 * qJD(1);
t102 = -t162 / 0.2e1 - t161 / 0.2e1;
t40 = (-t99 / 0.2e1 + t102) * pkin(3);
t142 = t40 * qJD(1);
t41 = t84 * t164 + t88 * t86;
t141 = t41 * qJD(1);
t42 = t86 * t164 - t88 * t84;
t140 = t42 * qJD(1);
t51 = t108 - t79 / 0.2e1;
t139 = t51 * qJD(1);
t49 = t51 * qJD(5);
t138 = t167 * qJD(1);
t134 = t167 * qJD(5);
t57 = t84 ^ 2 + t86 ^ 2;
t133 = t57 * qJD(1);
t132 = t84 * qJD(1);
t131 = t84 * qJD(3);
t130 = t86 * qJD(1);
t129 = t86 * qJD(3);
t89 = t100 ^ 2 - t99 ^ 2;
t128 = t89 * qJD(1);
t127 = t99 * qJD(3);
t126 = qJD(1) * t100;
t125 = t100 * qJD(3);
t122 = t111 * t149;
t121 = t167 * t149;
t120 = t93 * t99 * qJD(1);
t115 = t93 * t126;
t114 = t99 * t126;
t8 = t88 * t164;
t107 = t8 * qJD(1);
t3 = t109 + t166 / 0.2e1;
t92 = t96 * pkin(3) + pkin(4);
t76 = -t163 * t92 + t98 * t165;
t106 = -t3 * qJD(1) - t76 * qJD(3);
t104 = qJD(3) * t167 + t38 * qJD(5);
t2 = t110 + t175 / 0.2e1;
t77 = t163 * t165 + t98 * t92;
t101 = t2 * qJD(1) + t51 * qJD(2) + t77 * qJD(3);
t65 = t77 * qJD(5);
t64 = t76 * qJD(5);
t48 = t51 * qJD(3);
t39 = t164 / 0.2e1 + t102 * pkin(3);
t14 = -t38 * qJD(3) - t134;
t5 = 0.2e1 * t110 - t168;
t4 = t176 + 0.2e1 * t109;
t1 = [0, 0, 0, 0, t99 * t125, t89 * qJD(3), 0, 0, 0, t93 * t127, t93 * t125, t41 * qJD(3), t42 * qJD(3), t57 * qJD(4), t8 * qJD(3) + t13 * qJD(4), -t105 * t167, (qJD(3) + qJD(5)) * t169, 0, 0, 0, t10 * qJD(3) + t134 * t62, t11 * qJD(3) - t137 * t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, t114, t128, t125, -t127, 0, -t91 * t125 + t120, t91 * t127 + t115, t43 * qJD(3) + t141, -qJD(3) * t112 + t140, (t84 * t96 - t86 * t95) * t154, t39 * qJD(4) + (t112 * t95 + t43 * t96) * t154 + t107, -t174, t173, -t105, -t104, 0, t177 * qJD(3) + t5 * qJD(5) + t148, t178 * qJD(3) + t4 * qJD(5) + t147; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t133, t39 * qJD(3) + t145, 0, 0, 0, 0, 0, t49, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t174, t173, -t105, t14, 0, t5 * qJD(3) + t51 * qJD(4) + t177 * qJD(5) + t121, t4 * qJD(3) + t178 * qJD(5) - t122; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t127, -t125, -t129, t131, 0, (-t161 - t162) * t154, 0, 0, 0, 0, 0, -t104, t105; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, t105; 0, 0, 0, 0, -t114, -t128, 0, 0, 0, -t120, -t115, -t86 * qJD(4) - t141, t84 * qJD(4) - t140, 0, t40 * qJD(4) - t107, t174, -t173, 0, -t49, 0, -qJD(4) * t167 - t2 * qJD(5) - t148, t3 * qJD(5) - t147 + t171; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t65, t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t130, t132, 0, t142, 0, 0, 0, 0, 0, -t138, t170; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t139, 0, -t101 - t65, -t106 + t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t129, -t131, -t133, -t40 * qJD(3) - t145, 0, 0, 0, 0, 0, t104, -t105; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t130, -t132, 0, -t142, 0, 0, 0, 0, 0, t138, -t170; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t143, -t170; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t174, -t173, 0, t48, 0, t2 * qJD(3) - t38 * qJD(4) - t121, -t3 * qJD(3) + t122 + t171; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t139, 0, t101, t106; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t143, t170; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t1;
