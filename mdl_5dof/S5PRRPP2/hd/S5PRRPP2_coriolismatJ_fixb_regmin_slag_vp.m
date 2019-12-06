% Calculate minimal parameter regressor of coriolis matrix for
% S5PRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x17]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5PRRPP2_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP2_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP2_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP2_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:10:11
% EndTime: 2019-12-05 16:10:16
% DurationCPUTime: 1.18s
% Computational Cost: add. (1346->133), mult. (3224->215), div. (0->0), fcn. (3432->6), ass. (0->111)
t115 = cos(qJ(2));
t111 = sin(pkin(8));
t114 = cos(qJ(3));
t112 = sin(qJ(3));
t173 = cos(pkin(8));
t141 = t173 * t112;
t88 = t111 * t114 + t141;
t73 = t115 * t88;
t211 = t73 / 0.2e1;
t210 = t88 / 0.2e1;
t113 = sin(qJ(2));
t185 = t113 / 0.2e1;
t140 = t173 * t114;
t172 = t111 * t112;
t123 = t140 - t172;
t71 = t123 * t113;
t72 = t88 * t113;
t191 = t72 * t210 + t71 * t123 / 0.2e1;
t127 = t185 + t191;
t74 = t123 * t115;
t125 = -t113 * t115 + t71 * t74 + t72 * t73;
t195 = t125 * qJD(1);
t208 = t127 * qJD(4) + t195;
t122 = t185 - t191;
t207 = -qJD(4) * t122 - t195;
t85 = t88 ^ 2;
t202 = t123 ^ 2 + t85;
t206 = qJD(4) * t202;
t205 = t202 * qJD(2);
t204 = t111 / 0.2e1;
t203 = qJ(4) + pkin(6);
t139 = t203 * t172;
t98 = t203 * t114;
t91 = t173 * t98;
t201 = t91 - t139;
t61 = t111 * t98 + t141 * t203;
t135 = t123 * t201 + t61 * t88;
t198 = qJD(4) * t135;
t196 = t122 * qJD(2);
t194 = t125 * qJD(2);
t193 = t127 * qJD(2);
t190 = -qJD(1) * t122 + qJD(2) * t135;
t137 = t91 / 0.2e1;
t188 = t88 * pkin(4);
t102 = pkin(3) * t111 + qJ(5);
t187 = -t102 / 0.2e1;
t104 = -pkin(3) * t173 - pkin(4);
t186 = -t104 / 0.2e1;
t184 = -t115 / 0.2e1;
t183 = t112 * pkin(3);
t177 = qJD(3) * pkin(3);
t174 = t123 * qJ(5);
t170 = t115 * t112;
t107 = t183 / 0.2e1;
t19 = t107 + (pkin(4) / 0.2e1 + t186) * t88 - (qJ(5) / 0.2e1 + t102 / 0.2e1) * t123;
t167 = t19 * qJD(2);
t144 = -t173 / 0.2e1;
t121 = t123 * t204 + t144 * t88;
t35 = (-t112 / 0.2e1 + t121) * pkin(3);
t166 = t35 * qJD(2);
t146 = -t73 / 0.2e1;
t38 = t146 + t211;
t165 = t38 * qJD(1);
t162 = t71 * qJD(3);
t161 = t85 * qJD(2);
t160 = t123 * qJD(2);
t159 = t123 * qJD(3);
t158 = t88 * qJD(2);
t157 = t88 * qJD(5);
t156 = qJD(2) * t114;
t101 = -t112 ^ 2 + t114 ^ 2;
t155 = t101 * qJD(2);
t154 = t112 * qJD(3);
t153 = t113 * qJD(2);
t152 = t114 * qJD(3);
t151 = t115 * qJD(2);
t150 = pkin(2) * t112 * qJD(2);
t149 = pkin(2) * t156;
t148 = t123 * t158;
t106 = -pkin(3) * t114 - pkin(2);
t143 = t112 * t156;
t142 = t201 * t74 + t61 * t73;
t134 = (t123 * t74 + t73 * t88) * qJD(2);
t45 = -t174 + t183 + t188;
t116 = t45 * t184;
t124 = t186 * t73 + t187 * t74;
t2 = t116 + t124;
t44 = -pkin(4) * t123 - qJ(5) * t88 + t106;
t5 = t44 * t45;
t133 = qJD(1) * t2 + qJD(2) * t5;
t12 = t106 * t183;
t120 = t144 * t73 + t204 * t74;
t3 = (t170 / 0.2e1 + t120) * pkin(3);
t132 = -t3 * qJD(1) + t12 * qJD(2);
t13 = -t123 * t45 + t44 * t88;
t131 = t13 * qJD(2);
t14 = -t123 * t44 - t45 * t88;
t117 = t140 / 0.2e1 - t172 / 0.2e1;
t40 = (-t123 / 0.2e1 + t117) * t115;
t130 = -qJD(1) * t40 + qJD(2) * t14;
t58 = t137 - t91 / 0.2e1;
t126 = qJD(2) * t58 + qJD(3) * t102;
t43 = 0.2e1 * t211;
t42 = t146 - t211;
t41 = t115 * t117 - t123 * t184;
t36 = 0.2e1 * t137 - t139;
t34 = pkin(3) * t121 + t107;
t20 = -t123 * t187 + t104 * t210 + t107 - t174 / 0.2e1 + t188 / 0.2e1;
t4 = (-t170 / 0.2e1 + t120) * pkin(3);
t1 = t116 - t124;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t194, 0, 0, 0, t194; 0, 0, -t153, -t151, 0, 0, 0, 0, 0, -t114 * t153 - t115 * t154, t112 * t153 - t115 * t152, t134, (t106 * t113 + t142) * qJD(2) + t4 * qJD(3) + t208, qJD(3) * t42 - t123 * t153, t134, qJD(3) * t41 - t153 * t88, (t113 * t44 + t142) * qJD(2) + t1 * qJD(3) + t43 * qJD(5) + t208; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t112 * t151 - t113 * t152, t113 * t154 - t114 * t151, 0, t4 * qJD(2) + (-t111 * t72 - t173 * t71) * t177, qJD(2) * t42 - t162, 0, qJD(2) * t41 - qJD(3) * t72, t1 * qJD(2) + (-t102 * t72 + t104 * t71) * qJD(3) + t71 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t193, 0, 0, 0, t193; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t43 + t162; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3) * t3 + t207, 0, 0, -t40 * qJD(3), qJD(3) * t2 - qJD(5) * t38 + t207; 0, 0, 0, 0, t112 * t152, t101 * qJD(3), 0, 0, 0, -pkin(2) * t154, -pkin(2) * t152, t206, qJD(3) * t12 + t198, qJD(3) * t13 + t123 * t157, t206, qJD(3) * t14 + qJD(5) * t85, qJD(3) * t5 - t157 * t44 + t198; 0, 0, 0, 0, t143, t155, t152, -t154, 0, -pkin(6) * t152 - t150, pkin(6) * t154 - t149, (-t111 * t88 - t123 * t173) * t177, (-t111 * t61 - t173 * t201) * t177 + t34 * qJD(4) + t132, -qJD(3) * t201 + t131, (-t102 * t88 + t104 * t123) * qJD(3) + qJD(5) * t123, -qJD(3) * t61 + t130, (-t102 * t61 + t104 * t201) * qJD(3) + t20 * qJD(4) + t36 * qJD(5) + t133; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t205, qJD(3) * t34 + t190, 0, t205, 0, qJD(3) * t20 + t190; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t148, t159, t161, qJD(3) * t36 - t158 * t44 - t165; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3 * qJD(2), 0, 0, t40 * qJD(2), -t2 * qJD(2); 0, 0, 0, 0, -t143, -t155, 0, 0, 0, t150, t149, 0, qJD(4) * t35 - t132, -qJD(4) * t88 - t131, 0, qJD(4) * t123 - t130, -qJD(4) * t19 + qJD(5) * t58 - t133; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(5), t102 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t166, -t158, 0, t160, -t167; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3), t126; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t196, 0, 0, 0, t196; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t205, -qJD(3) * t35 - t190, t88 * qJD(3), -t205, -t159, qJD(3) * t19 - t157 - t190; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t166, t158, 0, -t160, t167; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t158; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38 * qJD(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t148, 0, -t161, t165 - t58 * qJD(3) + (qJD(2) * t44 + qJD(4)) * t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3), -t126; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t158; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t6;
