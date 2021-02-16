% Calculate minimal parameter regressor of coriolis matrix for
% S5PRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x22]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 16:24
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5PRRRP3_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP3_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP3_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP3_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 16:23:30
% EndTime: 2021-01-15 16:23:34
% DurationCPUTime: 1.15s
% Computational Cost: add. (1654->128), mult. (3304->167), div. (0->0), fcn. (3432->4), ass. (0->114)
t143 = qJD(3) + qJD(4);
t181 = cos(qJ(4));
t142 = t181 * pkin(3);
t112 = t142 + pkin(4);
t122 = t142 / 0.2e1 - t112 / 0.2e1;
t115 = sin(qJ(4));
t116 = sin(qJ(3));
t117 = cos(qJ(3));
t104 = t115 * t117 + t181 * t116;
t185 = pkin(6) + pkin(7);
t108 = t185 * t116;
t109 = t185 * t117;
t189 = t181 * t108 + t115 * t109;
t192 = t104 * qJ(5) + t189;
t193 = t143 * t192;
t102 = t115 * t116 - t181 * t117;
t164 = t102 * qJ(5);
t107 = t181 * t109;
t161 = t115 * t108;
t188 = -t107 + t161;
t48 = t188 + t164;
t191 = t143 * t189;
t190 = t143 * t102;
t64 = t143 * t104;
t101 = t102 ^ 2;
t187 = t104 ^ 2;
t184 = t48 * pkin(4);
t183 = -t107 / 0.2e1;
t180 = pkin(3) * t115;
t179 = t102 * pkin(4);
t178 = t104 * pkin(4);
t177 = t116 * pkin(3);
t175 = pkin(4) * qJD(4);
t113 = -t117 * pkin(3) - pkin(2);
t79 = t113 + t179;
t61 = t79 * t104;
t11 = t102 * t48 + t104 * t192;
t170 = qJD(2) * t11;
t80 = t177 + t178;
t17 = t102 * t80 + t61;
t169 = qJD(2) * t17;
t60 = t79 * t102;
t18 = t104 * t80 - t60;
t168 = qJD(2) * t18;
t19 = -t102 * t178 - t61;
t167 = qJD(2) * t19;
t20 = -pkin(4) * t187 + t60;
t166 = qJD(2) * t20;
t165 = qJD(4) * t48;
t163 = t112 * t104;
t162 = t115 * t102;
t131 = -t162 / 0.2e1;
t132 = -t163 / 0.2e1;
t86 = -t178 / 0.2e1;
t22 = t132 + t86 + (t131 - t116 / 0.2e1) * pkin(3);
t157 = t22 * qJD(2);
t35 = (-pkin(4) / 0.2e1 - t122) * t102;
t156 = t35 * qJD(2);
t53 = t101 - t187;
t155 = t53 * qJD(2);
t56 = t102 * t177 + t113 * t104;
t154 = t56 * qJD(2);
t57 = -t113 * t102 + t104 * t177;
t153 = t57 * qJD(2);
t69 = t101 + t187;
t152 = t69 * qJD(2);
t76 = t183 + t107 / 0.2e1;
t151 = t76 * qJD(2);
t150 = qJD(2) * t117;
t149 = qJD(4) * t113;
t148 = t102 * qJD(2);
t147 = t104 * qJD(2);
t94 = t104 * qJD(5);
t110 = -t116 ^ 2 + t117 ^ 2;
t146 = t110 * qJD(2);
t145 = t116 * qJD(3);
t144 = t117 * qJD(3);
t141 = pkin(2) * t116 * qJD(2);
t140 = pkin(2) * t150;
t139 = t104 * t175;
t138 = qJD(4) * t180;
t137 = pkin(4) * t147;
t135 = t113 * t148;
t134 = t113 * t147;
t133 = t116 * t150;
t130 = t181 * qJD(3);
t129 = t181 * qJD(4);
t126 = pkin(3) * t129;
t9 = pkin(4) * t61;
t124 = qJD(2) * t9;
t8 = t79 * t80;
t123 = t8 * qJD(2);
t120 = t122 * t104;
t118 = t122 * t48;
t3 = -t184 / 0.2e1 - t118;
t87 = t178 / 0.2e1;
t37 = t87 + t120;
t84 = (t142 - t112) * t180;
t119 = -qJD(1) * t37 - qJD(2) * t3 - qJD(3) * t84;
t55 = 0.2e1 * t183 + t161;
t90 = t102 * qJD(5);
t74 = t102 * t147;
t68 = t76 * qJD(3);
t67 = t76 * qJD(4);
t59 = qJD(3) * t180 - t151;
t58 = pkin(3) * t130;
t47 = -t143 * t180 + t151;
t46 = (-t130 - t129) * pkin(3);
t36 = t86 + t120;
t34 = t179 / 0.2e1 - t122 * t102;
t32 = t55 + t164;
t21 = pkin(3) * t131 + t132 + t177 / 0.2e1 + t87;
t2 = t184 / 0.2e1 - t118;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t145, -t144, 0, 0, 0, 0, 0, -t64, t190, -t64, t190, 0, (-pkin(3) * t162 - t163) * qJD(3) + t36 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t64, t190, -t64, t190, 0, qJD(3) * t36 - t139; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, t116 * t144, t110 * qJD(3), 0, 0, 0, -pkin(2) * t145, -pkin(2) * t144, -t102 * t64, t143 * t53, 0, 0, 0, qJD(3) * t56 + t104 * t149, qJD(3) * t57 - t102 * t149, qJD(3) * t17 - qJD(4) * t19, qJD(3) * t18 - qJD(4) * t20, qJD(5) * t69, qJD(3) * t8 + qJD(4) * t9 + qJD(5) * t11; 0, 0, 0, 0, t133, t146, t144, -t145, 0, -pkin(6) * t144 - t141, pkin(6) * t145 - t140, -t74, t155, -t190, -t64, 0, qJD(3) * t188 + t55 * qJD(4) + t154, t153 + t191, qJD(3) * t48 + qJD(4) * t32 + t169, t168 + t193, (t112 * t102 - t104 * t180) * qJD(3) + t34 * qJD(4), (t112 * t48 - t180 * t192) * qJD(3) + t2 * qJD(4) + t21 * qJD(5) + t123; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t74, t155, -t190, -t64, 0, t55 * qJD(3) + qJD(4) * t188 + t134, -t135 + t191, qJD(3) * t32 + t165 - t167, -t166 + t193, qJD(3) * t34 + t102 * t175, pkin(4) * t165 + qJD(3) * t2 + t124; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t152, qJD(3) * t21 + t170; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4) * t37; 0, 0, 0, 0, -t133, -t146, 0, 0, 0, t141, t140, t74, -t155, 0, 0, 0, t67 - t154, -t153, t67 - t94 - t169, t90 - t168, qJD(4) * t35, qJD(4) * t3 + qJD(5) * t22 - t123; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t138, -t126, -t138, -t126, 0, t84 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, t46, t47, t46, t156, -pkin(4) * t138 - t119; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t147, t148, 0, t157; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3) * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74, -t155, 0, 0, 0, -t68 - t134, t135, -t68 - t94 + t167, t90 + t166, -qJD(3) * t35, -pkin(4) * t94 - qJD(3) * t3 - t124; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, t58, t59, t58, -t156, t119; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t147, t148, 0, -t137; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64, -t190, -t152, -qJD(3) * t22 + t139 - t170; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t147, -t148, 0, -t157; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t147, -t148, 0, t137; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t1;
