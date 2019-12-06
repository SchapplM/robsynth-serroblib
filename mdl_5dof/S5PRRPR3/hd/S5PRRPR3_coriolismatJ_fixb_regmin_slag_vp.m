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
% cmat_reg [(5*%NQJ)%x20]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 16:19:44
% EndTime: 2019-12-05 16:19:47
% DurationCPUTime: 0.82s
% Computational Cost: add. (1199->98), mult. (2496->140), div. (0->0), fcn. (2782->6), ass. (0->90)
t150 = cos(qJ(5));
t141 = -qJ(4) - pkin(6);
t93 = sin(qJ(3));
t84 = t141 * t93;
t94 = cos(qJ(3));
t85 = t141 * t94;
t90 = sin(pkin(9));
t91 = cos(pkin(9));
t155 = t84 * t91 + t85 * t90;
t77 = t90 * t94 + t91 * t93;
t159 = -pkin(7) * t77 + t155;
t162 = t150 * t159;
t102 = -t162 / 0.2e1;
t56 = -t84 * t90 + t85 * t91;
t76 = -t90 * t93 + t91 * t94;
t40 = pkin(7) * t76 - t56;
t164 = t150 * t40;
t103 = -t164 / 0.2e1;
t92 = sin(qJ(5));
t165 = t92 * t40;
t167 = -t162 + t165;
t163 = t92 * t159;
t166 = -t164 - t163;
t142 = t92 * t76;
t73 = t150 * t77;
t154 = t73 + t142;
t104 = -t150 * t76 + t77 * t92;
t157 = t104 * qJD(2);
t161 = t154 * t157;
t156 = t104 ^ 2 - t154 ^ 2;
t160 = t156 * qJD(2);
t158 = qJD(4) * t104;
t125 = t104 * qJD(5);
t98 = qJD(3) * t104 + t125;
t101 = t73 / 0.2e1;
t153 = pkin(3) * t90;
t151 = t93 * pkin(3);
t149 = t76 * t90;
t148 = t77 * t91;
t140 = qJD(3) * pkin(3);
t110 = -pkin(3) * t94 - pkin(2);
t60 = -pkin(4) * t76 + t110;
t136 = qJD(2) * t60;
t135 = qJD(2) * t94;
t61 = pkin(4) * t77 + t151;
t10 = t104 * t61 + t154 * t60;
t134 = t10 * qJD(2);
t11 = -t104 * t60 + t154 * t61;
t133 = t11 * qJD(2);
t18 = -t155 * t77 - t76 * t56;
t131 = t18 * qJD(2);
t22 = 0.2e1 * t101 + t142;
t129 = t22 * qJD(2);
t96 = t149 / 0.2e1 - t148 / 0.2e1;
t32 = (-t93 / 0.2e1 + t96) * pkin(3);
t128 = t32 * qJD(2);
t45 = t101 - t73 / 0.2e1;
t127 = t45 * qJD(2);
t43 = t45 * qJD(5);
t126 = t154 * qJD(2);
t122 = t154 * qJD(5);
t51 = t76 ^ 2 + t77 ^ 2;
t121 = t51 * qJD(2);
t86 = -t93 ^ 2 + t94 ^ 2;
t120 = t86 * qJD(2);
t119 = t93 * qJD(3);
t118 = t94 * qJD(3);
t117 = pkin(2) * t93 * qJD(2);
t116 = pkin(2) * t135;
t113 = t104 * t136;
t112 = t154 * t136;
t111 = t93 * t135;
t9 = t110 * t151;
t100 = t9 * qJD(2);
t3 = t102 + t162 / 0.2e1;
t88 = pkin(3) * t91 + pkin(4);
t70 = -t150 * t88 + t153 * t92;
t99 = -qJD(2) * t3 - qJD(3) * t70;
t97 = qJD(3) * t154 + qJD(5) * t22;
t2 = t103 + t164 / 0.2e1;
t71 = t150 * t153 + t88 * t92;
t95 = qJD(1) * t45 + qJD(2) * t2 + qJD(3) * t71;
t63 = t71 * qJD(5);
t62 = t70 * qJD(5);
t42 = t45 * qJD(3);
t31 = t151 / 0.2e1 + t96 * pkin(3);
t13 = -qJD(3) * t22 - t122;
t5 = 0.2e1 * t103 - t163;
t4 = t165 + 0.2e1 * t102;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t119, -t118, 0, (-t148 + t149) * t140, 0, 0, 0, 0, 0, -t97, t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, t93 * t118, t86 * qJD(3), 0, 0, 0, -pkin(2) * t119, -pkin(2) * t118, t51 * qJD(4), qJD(3) * t9 + qJD(4) * t18, -t98 * t154, (qJD(3) + qJD(5)) * t156, 0, 0, 0, qJD(3) * t10 + t122 * t60, qJD(3) * t11 - t125 * t60; 0, 0, 0, 0, t111, t120, t118, -t119, 0, -pkin(6) * t118 - t117, pkin(6) * t119 - t116, (-t76 * t91 - t77 * t90) * t140, t31 * qJD(4) + (t155 * t90 + t56 * t91) * t140 + t100, -t161, t160, -t98, -t97, 0, qJD(3) * t166 + t5 * qJD(5) + t134, qJD(3) * t167 + t4 * qJD(5) + t133; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t121, qJD(3) * t31 + t131, 0, 0, 0, 0, 0, t43, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t161, t160, -t98, t13, 0, t5 * qJD(3) + t45 * qJD(4) + qJD(5) * t166 + t112, t4 * qJD(3) + qJD(5) * t167 - t113; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43, 0; 0, 0, 0, 0, -t111, -t120, 0, 0, 0, t117, t116, 0, qJD(4) * t32 - t100, t161, -t160, 0, -t43, 0, -qJD(4) * t154 - qJD(5) * t2 - t134, qJD(5) * t3 - t133 + t158; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t63, t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t128, 0, 0, 0, 0, 0, -t126, t157; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t127, 0, -t63 - t95, t62 - t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t121, -qJD(3) * t32 - t131, 0, 0, 0, 0, 0, t97, -t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t128, 0, 0, 0, 0, 0, t126, -t157; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t129, -t157; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t161, -t160, 0, t42, 0, qJD(3) * t2 - qJD(4) * t22 - t112, -qJD(3) * t3 + t113 + t158; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t127, 0, t95, t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t129, t157; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t1;
