% Calculate minimal parameter regressor of coriolis matrix for
% S4RRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% cmat_reg [(4*%NQJ)%x21]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 14:30
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4RRRP4_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP4_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP4_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP4_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 14:30:25
% EndTime: 2021-01-15 14:30:28
% DurationCPUTime: 0.83s
% Computational Cost: add. (1066->116), mult. (2119->158), div. (0->0), fcn. (2116->4), ass. (0->108)
t107 = qJD(2) + qJD(3);
t137 = cos(qJ(3));
t143 = pkin(5) + pkin(6);
t86 = sin(qJ(2));
t78 = t143 * t86;
t87 = cos(qJ(2));
t79 = t143 * t87;
t85 = sin(qJ(3));
t149 = t137 * t78 + t85 * t79;
t74 = t137 * t86 + t85 * t87;
t14 = t74 * qJ(4) + t149;
t152 = t107 * t14;
t102 = t137 * pkin(2);
t151 = -t102 / 0.2e1;
t72 = -t137 * t87 + t85 * t86;
t125 = t72 * qJ(4);
t136 = t85 * t78;
t77 = t137 * t79;
t148 = -t77 + t136;
t27 = t148 + t125;
t150 = t107 * t149;
t42 = t107 * t74;
t147 = t74 ^ 2;
t146 = -pkin(3) / 0.2e1;
t145 = -t77 / 0.2e1;
t82 = t102 + pkin(3);
t144 = -t82 / 0.2e1;
t142 = pkin(2) * t85;
t141 = t27 * pkin(3);
t140 = t72 * pkin(3);
t139 = t74 * pkin(3);
t138 = t86 * pkin(2);
t83 = -pkin(2) * t87 - pkin(1);
t53 = t83 + t140;
t40 = t53 * t74;
t133 = pkin(3) * qJD(3);
t6 = pkin(3) * t40;
t131 = qJD(1) * t6;
t7 = t14 * t74 + t27 * t72;
t130 = qJD(1) * t7;
t54 = t138 + t139;
t8 = t54 * t72 + t40;
t129 = qJD(1) * t8;
t39 = t53 * t72;
t9 = t54 * t74 - t39;
t128 = qJD(1) * t9;
t5 = t53 * t54;
t126 = t5 * qJD(1);
t10 = -t72 * t139 - t40;
t124 = qJD(1) * t10;
t11 = -pkin(3) * t147 + t39;
t123 = qJD(1) * t11;
t122 = qJD(1) * t87;
t121 = qJD(3) * t27;
t120 = qJD(3) * t83;
t98 = -t85 * t72 / 0.2e1;
t13 = (t144 + t146) * t74 + (t98 - t86 / 0.2e1) * pkin(2);
t119 = t13 * qJD(1);
t90 = t151 + t82 / 0.2e1;
t17 = (t146 + t90) * t72;
t118 = t17 * qJD(1);
t71 = t72 ^ 2;
t32 = t71 - t147;
t117 = t32 * qJD(1);
t35 = t72 * t138 + t74 * t83;
t116 = t35 * qJD(1);
t36 = t74 * t138 - t72 * t83;
t115 = t36 * qJD(1);
t47 = t71 + t147;
t114 = t47 * qJD(1);
t50 = t145 + t77 / 0.2e1;
t113 = t50 * qJD(1);
t112 = t72 * qJD(1);
t111 = t74 * qJD(1);
t64 = t74 * qJD(4);
t80 = -t86 ^ 2 + t87 ^ 2;
t110 = t80 * qJD(1);
t109 = t86 * qJD(2);
t108 = t87 * qJD(2);
t106 = pkin(1) * t86 * qJD(1);
t105 = pkin(1) * t122;
t104 = qJD(3) * t142;
t103 = pkin(3) * t111;
t101 = t83 * t112;
t100 = t83 * t111;
t99 = t86 * t122;
t95 = t137 * qJD(2);
t94 = t137 * qJD(3);
t93 = pkin(2) * t94;
t88 = (-t144 + t151) * t27;
t2 = -t141 / 0.2e1 + t88;
t57 = (t102 - t82) * t142;
t91 = -qJD(1) * t2 - qJD(2) * t57;
t34 = 0.2e1 * t145 + t136;
t60 = t72 * qJD(4);
t48 = t72 * t111;
t46 = t50 * qJD(2);
t45 = t50 * qJD(3);
t41 = t107 * t72;
t38 = qJD(2) * t142 - t113;
t37 = pkin(2) * t95;
t26 = -t107 * t142 + t113;
t25 = (-t95 - t94) * pkin(2);
t16 = t140 / 0.2e1 + t90 * t72;
t15 = t34 + t125;
t12 = pkin(2) * t98 + t74 * t144 + t138 / 0.2e1 + t139 / 0.2e1;
t1 = t141 / 0.2e1 + t88;
t3 = [0, 0, 0, t86 * t108, t80 * qJD(2), 0, 0, 0, -pkin(1) * t109, -pkin(1) * t108, -t72 * t42, t107 * t32, 0, 0, 0, qJD(2) * t35 + t74 * t120, qJD(2) * t36 - t72 * t120, qJD(2) * t8 - qJD(3) * t10, qJD(2) * t9 - qJD(3) * t11, qJD(4) * t47, qJD(2) * t5 + qJD(3) * t6 + qJD(4) * t7; 0, 0, 0, t99, t110, t108, -t109, 0, -pkin(5) * t108 - t106, pkin(5) * t109 - t105, -t48, t117, -t41, -t42, 0, qJD(2) * t148 + t34 * qJD(3) + t116, t115 + t150, qJD(2) * t27 + qJD(3) * t15 + t129, t128 + t152, (-t74 * t142 + t82 * t72) * qJD(2) + t16 * qJD(3), t126 + (-t14 * t142 + t27 * t82) * qJD(2) + t1 * qJD(3) + t12 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48, t117, -t41, -t42, 0, t34 * qJD(2) + qJD(3) * t148 + t100, -t101 + t150, qJD(2) * t15 + t121 - t124, -t123 + t152, qJD(2) * t16 + t72 * t133, pkin(3) * t121 + qJD(2) * t1 + t131; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t114, qJD(2) * t12 + t130; 0, 0, 0, -t99, -t110, 0, 0, 0, t106, t105, t48, -t117, 0, 0, 0, t45 - t116, -t115, t45 - t64 - t129, t60 - t128, qJD(3) * t17, qJD(3) * t2 + qJD(4) * t13 - t126; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t104, -t93, -t104, -t93, 0, t57 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, t25, t26, t25, t118, -pkin(3) * t104 - t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t111, t112, 0, t119; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, -t117, 0, 0, 0, -t46 - t100, t101, -t46 - t64 + t124, t60 + t123, -qJD(2) * t17, -pkin(3) * t64 - qJD(2) * t2 - t131; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, t37, t38, t37, -t118, t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t111, t112, 0, -t103; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, -t41, -t114, -qJD(2) * t13 + t74 * t133 - t130; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t111, -t112, 0, -t119; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t111, -t112, 0, t103; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t3;
