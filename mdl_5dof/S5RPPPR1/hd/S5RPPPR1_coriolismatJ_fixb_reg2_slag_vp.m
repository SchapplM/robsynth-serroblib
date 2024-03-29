% Calculate inertial parameters regressor of coriolis matrix for
% S5RPPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% 
% Output:
% cmat_reg [(5*5)x(5*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 09:13
% Revision: 008671b0a00594318b890887636eaaff83cd5e2f (2021-12-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RPPPR1_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR1_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 09:12:50
% EndTime: 2022-01-20 09:12:52
% DurationCPUTime: 1.32s
% Computational Cost: add. (2113->107), mult. (4421->176), div. (0->0), fcn. (4621->8), ass. (0->109)
t112 = sin(pkin(9));
t171 = sin(qJ(5));
t130 = t171 * t112;
t114 = cos(pkin(9));
t172 = cos(qJ(5));
t131 = t172 * t114;
t93 = t130 - t131;
t180 = -t93 / 0.2e1;
t129 = t171 * t114;
t132 = t172 * t112;
t94 = t132 + t129;
t179 = -t94 / 0.2e1;
t175 = t114 ^ 2;
t176 = t112 ^ 2;
t178 = t175 + t176;
t113 = sin(pkin(8));
t78 = t94 * t113;
t81 = t93 * t113;
t168 = t78 * t179 + t81 * t180;
t77 = t78 ^ 2;
t177 = t81 ^ 2;
t173 = t113 / 0.2e1;
t115 = cos(pkin(8));
t138 = qJD(1) * t115;
t71 = t81 * t138;
t75 = t81 * qJD(5);
t167 = -t71 + t75;
t107 = sin(pkin(7)) * pkin(1) + qJ(3);
t120 = -cos(pkin(7)) * pkin(1) - t115 * pkin(3) - pkin(2);
t117 = -t113 * qJ(4) + t120;
t65 = t114 * t115 * t107 + t112 * t117;
t158 = t112 * t107;
t116 = (-pkin(4) - t158) * t115 + ((-pkin(6) - qJ(4)) * t113 + t120) * t114;
t43 = -t112 * t113 * pkin(6) + t65;
t25 = -t172 * t116 + t171 * t43;
t26 = t171 * t116 + t172 * t43;
t8 = -t25 * t81 - t26 * t78;
t166 = qJD(1) * t8;
t157 = t113 * t115;
t80 = t94 * t115;
t83 = t93 * t115;
t28 = t83 * t81 / 0.2e1 + t80 * t78 / 0.2e1 - t157 / 0.2e1;
t165 = qJD(1) * t28;
t32 = t83 * t78 - t80 * t81;
t164 = qJD(1) * t32;
t40 = -t113 * t78 - t80 * t115;
t163 = qJD(1) * t40;
t41 = -t113 * t81 - t83 * t115;
t162 = qJD(1) * t41;
t87 = (pkin(4) * t112 + t107) * t113;
t13 = -t25 * t115 - t87 * t78;
t156 = t13 * qJD(1);
t14 = -t26 * t115 + t87 * t81;
t155 = t14 * qJD(1);
t64 = t114 * t117 - t115 * t158;
t30 = (t112 * t65 + t114 * t64) * t113;
t153 = t30 * qJD(1);
t33 = t173 - t168;
t152 = t33 * qJD(1);
t35 = t77 - t177;
t151 = t35 * qJD(1);
t44 = t77 + t177;
t150 = t44 * qJD(1);
t118 = -t129 / 0.2e1 - t132 / 0.2e1;
t51 = (t179 + t118) * t115;
t149 = t51 * qJD(1);
t119 = -t131 / 0.2e1 + t130 / 0.2e1;
t52 = (t93 / 0.2e1 + t119) * t115;
t148 = t52 * qJD(1);
t124 = t176 / 0.2e1 + t175 / 0.2e1;
t67 = (-0.1e1 / 0.2e1 + t124) * t157;
t147 = t67 * qJD(1);
t110 = t113 ^ 2;
t103 = t115 ^ 2 + t110;
t76 = t103 * t107;
t146 = t76 * qJD(1);
t145 = t78 * qJD(1);
t144 = t78 * qJD(5);
t143 = t81 * qJD(1);
t85 = (0.1e1 / 0.2e1 + t124) * t113;
t142 = t85 * qJD(1);
t90 = t178 * t110;
t141 = t90 * qJD(1);
t91 = t103 * t112;
t140 = t91 * qJD(1);
t92 = t103 * t114;
t139 = t92 * qJD(1);
t137 = qJD(4) * t115;
t136 = qJD(5) * t115;
t135 = t103 * qJD(1);
t134 = t78 * t143;
t133 = t78 * t75;
t128 = t113 * t138;
t127 = t113 * t137;
t126 = t112 * t128;
t125 = t114 * t128;
t7 = t87 * t113 + t25 * t80 - t26 * t83;
t123 = t7 * qJD(1) + t28 * qJD(2);
t27 = t110 * t107 + (-t112 * t64 + t114 * t65) * t115;
t122 = -t27 * qJD(1) - t67 * qJD(2);
t86 = t173 - t178 * t113 / 0.2e1;
t70 = t78 * t138;
t66 = t67 * qJD(3);
t54 = (t94 / 0.2e1 + t118) * t115;
t53 = (t180 + t119) * t115;
t42 = t70 - t144;
t34 = t173 + t168;
t1 = qJD(3) * t28;
t2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t103 * qJD(3), t76 * qJD(3), 0, 0, 0, 0, 0, 0, t91 * qJD(3) + t114 * t127, t92 * qJD(3) - t112 * t127, t90 * qJD(4), qJD(3) * t27 - qJD(4) * t30, t133, t35 * qJD(5), t78 * t136, -t133, -t81 * t136, 0, -t40 * qJD(3) - t14 * qJD(5) - t81 * t137, t41 * qJD(3) + t13 * qJD(5) - t78 * t137, qJD(3) * t32 + qJD(4) * t44, qJD(3) * t7 + qJD(4) * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t135, t146, 0, 0, 0, 0, 0, 0, t140, t139, 0, t86 * qJD(4) - t122, 0, 0, 0, 0, 0, 0, qJD(5) * t54 - t163, qJD(5) * t53 + t162, t164, (t80 * t93 - t83 * t94) * qJD(3) + t34 * qJD(4) + t123; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t125, -t126, t141, t86 * qJD(3) - t153, 0, 0, 0, 0, 0, 0, -t71, -t70, t150, qJD(3) * t34 + t166; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t134, t151, t42, -t134, t167, 0, qJD(3) * t54 - qJD(5) * t26 - t155, qJD(3) * t53 + qJD(5) * t25 + t156, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t147, 0, 0, 0, 0, 0, 0, 0, 0, 0, t165; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, t144, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t135, -t146, 0, 0, 0, 0, 0, 0, -t140, -t139, 0, -t85 * qJD(4) + t122, 0, 0, 0, 0, 0, 0, -qJD(5) * t51 + t163, -qJD(5) * t52 - t162, -t164, -qJD(4) * t33 - t123; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t147, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t165; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t142, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t152; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t94 * qJD(5) - t149, t93 * qJD(5) - t148, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t125, t126, -t141, t85 * qJD(3) + t153, 0, 0, 0, 0, 0, 0, -t167, t42, -t150, qJD(3) * t33 - t166; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t142, 0, 0, 0, 0, 0, 0, 0, 0, 0, t152; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t143, -t145, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t134, -t151, -t70, t134, t71, 0, t51 * qJD(3) + t81 * qJD(4) + t155, t52 * qJD(3) + qJD(4) * t78 - t156, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t149, t148, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t143, t145, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t2;
