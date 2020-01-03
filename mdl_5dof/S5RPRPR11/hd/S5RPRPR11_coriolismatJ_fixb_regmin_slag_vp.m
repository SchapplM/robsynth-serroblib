% Calculate minimal parameter regressor of coriolis matrix for
% S5RPRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x25]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RPRPR11_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR11_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR11_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR11_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:27:57
% EndTime: 2019-12-31 18:28:00
% DurationCPUTime: 1.03s
% Computational Cost: add. (1216->110), mult. (2399->139), div. (0->0), fcn. (2760->6), ass. (0->90)
t101 = qJD(3) - qJD(5);
t146 = cos(qJ(3));
t82 = sin(pkin(8));
t83 = cos(pkin(8));
t85 = sin(qJ(3));
t66 = -t146 * t83 + t82 * t85;
t68 = t146 * t82 + t83 * t85;
t84 = sin(qJ(5));
t86 = cos(qJ(5));
t150 = -t66 * t86 + t68 * t84;
t89 = t66 * t84 + t68 * t86;
t10 = t150 ^ 2 - t89 ^ 2;
t165 = t10 * qJD(1);
t138 = pkin(6) + qJ(2);
t74 = t138 * t83;
t91 = t138 * t82;
t45 = t146 * t91 + t74 * t85;
t28 = -pkin(7) * t68 + t45;
t46 = t146 * t74 - t85 * t91;
t90 = pkin(7) * t66 + t46;
t164 = t101 * (t28 * t84 + t86 * t90);
t163 = t101 * (t28 * t86 - t84 * t90);
t158 = qJD(2) * t89;
t157 = t89 * qJD(1);
t118 = t150 * qJD(5);
t156 = qJD(3) * t150 - t118;
t155 = qJD(2) * t150;
t154 = t150 * qJD(1);
t115 = t89 * qJD(5);
t153 = qJD(3) * t89 - t115;
t62 = t68 ^ 2;
t148 = -pkin(3) - pkin(4);
t147 = t66 * pkin(3);
t76 = t82 ^ 2 + t83 ^ 2;
t133 = t68 * qJ(4);
t78 = -t83 * pkin(2) - pkin(1);
t87 = -t78 + t133;
t33 = -t87 + t147;
t134 = t66 * qJ(4);
t42 = pkin(3) * t68 + t134;
t5 = t33 * t42;
t135 = t5 * qJD(1);
t17 = t148 * t66 + t87;
t27 = t148 * t68 - t134;
t7 = t150 * t27 - t17 * t89;
t132 = t7 * qJD(1);
t8 = t150 * t17 + t27 * t89;
t131 = t8 * qJD(1);
t11 = t33 * t68 + t42 * t66;
t127 = t11 * qJD(1);
t12 = t33 * t66 - t42 * t68;
t126 = t12 * qJD(1);
t14 = t45 * t68 - t46 * t66;
t125 = t14 * qJD(1);
t122 = t42 * qJD(1);
t61 = t66 ^ 2;
t31 = t61 - t62;
t121 = t31 * qJD(1);
t44 = t61 + t62;
t114 = t44 * qJD(1);
t113 = t45 * qJD(3);
t43 = t46 * qJD(3);
t112 = t62 * qJD(1);
t111 = t66 * qJD(1);
t54 = t66 * qJD(3);
t110 = t68 * qJD(1);
t56 = t68 * qJD(3);
t109 = t68 * qJD(4);
t72 = t76 * qJ(2);
t108 = t72 * qJD(1);
t107 = t76 * qJD(1);
t106 = t84 * qJD(3);
t105 = t84 * qJD(4);
t104 = t86 * qJD(3);
t103 = t86 * qJD(4);
t102 = qJD(3) * qJ(4);
t100 = t17 * t154;
t99 = t17 * t157;
t98 = t150 * t157;
t97 = t89 * t154;
t96 = t150 * t110;
t95 = t89 * t110;
t94 = t66 * t110;
t93 = t78 * t110;
t75 = t101 * t86;
t73 = t101 * t84;
t71 = qJ(4) * t86 + t148 * t84;
t70 = qJ(4) * t84 - t148 * t86;
t58 = t68 * qJD(2);
t1 = [0, 0, 0, 0, 0, t76 * qJD(2), t72 * qJD(2), -t66 * t56, t31 * qJD(3), 0, 0, 0, t78 * t56, -t78 * t54, qJD(3) * t11 - t109 * t66, t44 * qJD(2), qJD(3) * t12 + qJD(4) * t62, qJD(2) * t14 + qJD(3) * t5 - t109 * t33, t156 * t89, -t101 * t10, 0, 0, 0, qJD(3) * t7 + t109 * t150 + t115 * t17, qJD(3) * t8 + t109 * t89 - t118 * t17; 0, 0, 0, 0, 0, t107, t108, 0, 0, 0, 0, 0, 0, 0, 0, t114, 0, t125, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, -t94, t121, -t54, -t56, 0, -t43 + t93, -t111 * t78 + t113, -t43 + t127, (-t133 + t147) * qJD(3) - t66 * qJD(4), -t113 + t126, t135 + (-pkin(3) * t46 - qJ(4) * t45) * qJD(3) + t46 * qJD(4), t97, -t165, -t156, -t153, 0, t132 - t164, t131 - t163; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t94, -t54, t112, -t110 * t33 + t43, 0, 0, 0, 0, 0, t96, t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t98, t165, t156, t153, 0, t99 + t164, -t100 + t163; 0, 0, 0, 0, 0, -t107, -t108, 0, 0, 0, 0, 0, t56, -t54, t56, -t114, t54, qJD(3) * t42 - t109 - t125, 0, 0, 0, 0, 0, t153, -t156; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t110, -t111, t110, 0, t111, t122, 0, 0, 0, 0, 0, t157, -t154; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t110, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t157, t154; 0, 0, 0, 0, 0, 0, 0, t94, -t121, 0, 0, 0, -t58 - t93, (qJD(1) * t78 + qJD(2)) * t66, -t58 - t127, 0, -qJD(2) * t66 - t126, -qJD(2) * t42 - t135, -t97, t165, 0, 0, 0, -t132 - t158, -t131 + t155; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t110, t111, -t110, 0, -t111, -t122, 0, 0, 0, 0, 0, -t157, t154; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4), qJ(4) * qJD(4), 0, 0, 0, 0, 0, qJD(5) * t71 + t105, -qJD(5) * t70 + t103; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3), t102, 0, 0, 0, 0, 0, t106, t104; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t101 * t71, -t101 * t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t94, 0, -t112, (qJD(1) * t33 + qJD(2)) * t68, 0, 0, 0, 0, 0, -t96, -t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t110, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3), -t102, 0, 0, 0, 0, 0, -t73, -t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73, t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t98, -t165, 0, 0, 0, -t99 + t158, t100 - t155; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t157, -t154; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3) * t71 - t105, qJD(3) * t70 - t103; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t106, -t104; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t1;
