% Calculate minimal parameter regressor of coriolis matrix for
% S5RPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x22]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:30
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RPRRP3_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP3_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP3_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP3_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:30:21
% EndTime: 2022-01-23 09:30:24
% DurationCPUTime: 1.20s
% Computational Cost: add. (2034->133), mult. (3684->175), div. (0->0), fcn. (3812->6), ass. (0->117)
t144 = qJD(3) + qJD(4);
t117 = sin(qJ(4));
t118 = sin(qJ(3));
t119 = cos(qJ(3));
t184 = cos(qJ(4));
t107 = t117 * t119 + t184 * t118;
t111 = sin(pkin(8)) * pkin(1) + pkin(6);
t179 = pkin(7) + t111;
t102 = t179 * t118;
t103 = t179 * t119;
t193 = t184 * t102 + t117 * t103;
t50 = -t107 * qJ(5) - t193;
t195 = t144 * t50;
t143 = t184 * pkin(3);
t105 = t117 * t118 - t184 * t119;
t165 = t105 * qJ(5);
t162 = t117 * t102;
t89 = t184 * t103;
t192 = -t89 + t162;
t48 = t192 + t165;
t194 = t144 * t193;
t191 = t144 * t105;
t70 = t144 * t107;
t104 = t105 ^ 2;
t190 = t107 ^ 2;
t189 = pkin(4) / 0.2e1;
t187 = -t89 / 0.2e1;
t186 = t48 * pkin(4);
t114 = t143 + pkin(4);
t185 = -t114 / 0.2e1;
t183 = pkin(3) * t117;
t182 = t105 * pkin(4);
t181 = t107 * pkin(4);
t180 = t118 * pkin(3);
t177 = pkin(4) * qJD(4);
t112 = -cos(pkin(8)) * pkin(1) - pkin(2);
t109 = -pkin(3) * t119 + t112;
t77 = t109 + t182;
t65 = t77 * t107;
t11 = t48 * t105 - t50 * t107;
t172 = qJD(1) * t11;
t80 = t180 + t181;
t17 = t105 * t80 + t65;
t171 = qJD(1) * t17;
t64 = t77 * t105;
t18 = t107 * t80 - t64;
t170 = qJD(1) * t18;
t21 = -t105 * t181 - t65;
t169 = qJD(1) * t21;
t22 = -pkin(4) * t190 + t64;
t168 = qJD(1) * t22;
t67 = t187 + t89 / 0.2e1;
t167 = qJD(1) * t67;
t166 = qJD(4) * t48;
t164 = t114 * t105;
t163 = t114 * t107;
t160 = t117 * t105;
t30 = (t189 + t114 / 0.2e1) * t107 + (t118 / 0.2e1 + t160 / 0.2e1) * pkin(3);
t157 = t30 * qJD(1);
t127 = t143 / 0.2e1;
t123 = t127 + t189;
t43 = (t185 + t123) * t105;
t156 = t43 * qJD(1);
t57 = t105 * t180 + t107 * t109;
t155 = t57 * qJD(1);
t58 = -t105 * t109 + t107 * t180;
t154 = t58 * qJD(1);
t59 = t104 - t190;
t153 = t59 * qJD(1);
t71 = t104 + t190;
t152 = t71 * qJD(1);
t151 = qJD(1) * t119;
t150 = qJD(4) * t109;
t149 = t105 * qJD(1);
t148 = t107 * qJD(1);
t99 = t107 * qJD(5);
t110 = -t118 ^ 2 + t119 ^ 2;
t147 = t110 * qJD(1);
t146 = t118 * qJD(3);
t145 = t119 * qJD(3);
t142 = t107 * t177;
t141 = qJD(4) * t183;
t140 = pkin(4) * t148;
t139 = -pkin(3) * t105 / 0.2e1;
t137 = t109 * t149;
t136 = t109 * t148;
t135 = t112 * t118 * qJD(1);
t134 = t112 * t151;
t133 = t118 * t151;
t132 = t184 * qJD(3);
t131 = t184 * qJD(4);
t128 = pkin(3) * t131;
t9 = pkin(4) * t65;
t125 = qJD(1) * t9;
t8 = t77 * t80;
t124 = t8 * qJD(1);
t120 = (-t143 / 0.2e1 - t185) * t48;
t3 = -t186 / 0.2e1 + t120;
t79 = -t163 / 0.2e1;
t42 = t123 * t107 + t79;
t90 = (t143 - t114) * t183;
t121 = -qJD(1) * t3 - qJD(2) * t42 - qJD(3) * t90;
t36 = 0.2e1 * t187 + t162;
t95 = t105 * qJD(5);
t76 = t105 * t148;
t63 = t67 * qJD(3);
t62 = t67 * qJD(4);
t54 = qJD(3) * t183 - t167;
t53 = pkin(3) * t132;
t47 = -t144 * t183 + t167;
t46 = (-t132 - t131) * pkin(3);
t44 = t164 / 0.2e1 + t184 * t139 + t182 / 0.2e1;
t41 = t107 * t127 + t79 - t181 / 0.2e1;
t31 = t79 + t117 * t139 + t180 / 0.2e1 + t181 / 0.2e1;
t20 = t36 + t165;
t2 = t186 / 0.2e1 + t120;
t1 = [0, 0, 0, 0, t118 * t145, t110 * qJD(3), 0, 0, 0, t112 * t146, t112 * t145, -t105 * t70, t144 * t59, 0, 0, 0, qJD(3) * t57 + t107 * t150, qJD(3) * t58 - t105 * t150, qJD(3) * t17 - qJD(4) * t21, qJD(3) * t18 - qJD(4) * t22, qJD(5) * t71, qJD(3) * t8 + qJD(4) * t9 + qJD(5) * t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, t133, t147, t145, -t146, 0, -t111 * t145 + t135, t111 * t146 + t134, -t76, t153, -t191, -t70, 0, qJD(3) * t192 + t36 * qJD(4) + t155, t154 + t194, qJD(3) * t48 + qJD(4) * t20 + t171, t170 - t195, (-t107 * t183 + t164) * qJD(3) + t44 * qJD(4), (t114 * t48 + t183 * t50) * qJD(3) + t2 * qJD(4) + t31 * qJD(5) + t124; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t76, t153, -t191, -t70, 0, t36 * qJD(3) + qJD(4) * t192 + t136, -t137 + t194, qJD(3) * t20 + t166 - t169, -t168 - t195, qJD(3) * t44 + t105 * t177, pkin(4) * t166 + qJD(3) * t2 + t125; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t152, qJD(3) * t31 + t172; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t146, -t145, 0, 0, 0, 0, 0, -t70, t191, -t70, t191, 0, (-pkin(3) * t160 - t163) * qJD(3) + t41 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t70, t191, -t70, t191, 0, qJD(3) * t41 - t142; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, -t133, -t147, 0, 0, 0, -t135, -t134, t76, -t153, 0, 0, 0, t62 - t155, -t154, t62 - t99 - t171, t95 - t170, -qJD(4) * t43, qJD(4) * t3 - qJD(5) * t30 - t124; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4) * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t141, -t128, -t141, -t128, 0, t90 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, t46, t47, t46, -t156, -pkin(4) * t141 - t121; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t148, t149, 0, -t157; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t76, -t153, 0, 0, 0, -t63 - t136, t137, -t63 - t99 + t169, t95 + t168, qJD(3) * t43, -pkin(4) * t99 - qJD(3) * t3 - t125; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3) * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, t53, t54, t53, t156, t121; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t148, t149, 0, -t140; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70, -t191, -t152, qJD(3) * t30 + t142 - t172; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t148, -t149, 0, t157; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t148, -t149, 0, t140; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t1;
