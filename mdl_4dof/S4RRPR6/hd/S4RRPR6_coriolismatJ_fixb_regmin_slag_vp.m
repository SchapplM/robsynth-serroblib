% Calculate minimal parameter regressor of coriolis matrix for
% S4RRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% cmat_reg [(4*%NQJ)%x21]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 10:46
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4RRPR6_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR6_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR6_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR6_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 10:46:25
% EndTime: 2021-01-15 10:46:28
% DurationCPUTime: 0.75s
% Computational Cost: add. (948->97), mult. (1946->155), div. (0->0), fcn. (2166->6), ass. (0->93)
t133 = cos(qJ(4));
t120 = cos(pkin(7));
t127 = -qJ(3) - pkin(5);
t74 = sin(qJ(2));
t65 = t127 * t74;
t75 = cos(qJ(2));
t66 = t127 * t75;
t72 = sin(pkin(7));
t43 = -t120 * t66 + t72 * t65;
t57 = -t120 * t75 + t72 * t74;
t25 = t57 * pkin(6) - t43;
t144 = t133 * t25;
t83 = t144 / 0.2e1;
t59 = t120 * t74 + t72 * t75;
t85 = t120 * t65 + t72 * t66;
t77 = -t59 * pkin(6) + t85;
t137 = t133 * t77;
t73 = sin(qJ(4));
t145 = t73 * t25;
t147 = -t137 - t145;
t138 = t73 * t77;
t146 = t144 - t138;
t35 = t133 * t57 + t73 * t59;
t107 = t35 * qJD(4);
t11 = -t35 * qJD(2) - t107;
t128 = t73 * t57;
t56 = t133 * t59;
t136 = t56 - t128;
t139 = -t136 ^ 2 + t35 ^ 2;
t143 = t139 * qJD(1);
t82 = -t137 / 0.2e1;
t141 = qJD(1) * t35;
t140 = qJD(3) * t35;
t108 = t136 * qJD(1);
t81 = t56 / 0.2e1;
t135 = pkin(2) * t72;
t134 = t74 * pkin(2);
t126 = qJD(2) * pkin(2);
t70 = -t75 * pkin(2) - pkin(1);
t7 = t70 * t134;
t123 = t7 * qJD(1);
t44 = t57 * pkin(3) + t70;
t45 = t59 * pkin(3) + t134;
t8 = t136 * t44 + t45 * t35;
t122 = t8 * qJD(1);
t9 = t136 * t45 - t35 * t44;
t121 = t9 * qJD(1);
t119 = qJD(1) * t44;
t118 = qJD(1) * t75;
t12 = -t43 * t57 - t59 * t85;
t116 = t12 * qJD(1);
t14 = 0.2e1 * t81 - t128;
t114 = t14 * qJD(1);
t76 = -t72 * t57 / 0.2e1 - t120 * t59 / 0.2e1;
t24 = (-t74 / 0.2e1 + t76) * pkin(2);
t113 = t24 * qJD(1);
t29 = t57 * t134 + t70 * t59;
t112 = t29 * qJD(1);
t30 = t59 * t134 - t70 * t57;
t111 = t30 * qJD(1);
t33 = t81 - t56 / 0.2e1;
t110 = t33 * qJD(1);
t109 = t33 * qJD(4);
t104 = t136 * qJD(4);
t39 = t57 ^ 2 + t59 ^ 2;
t103 = t39 * qJD(1);
t102 = t57 * qJD(1);
t101 = t59 * qJD(1);
t67 = -t74 ^ 2 + t75 ^ 2;
t100 = t67 * qJD(1);
t99 = t74 * qJD(2);
t98 = t75 * qJD(2);
t97 = pkin(1) * t74 * qJD(1);
t96 = pkin(1) * t118;
t95 = t35 * t108;
t94 = t136 * t141;
t93 = t35 * t119;
t92 = t136 * t119;
t91 = t74 * t118;
t1 = t83 - t144 / 0.2e1;
t69 = t120 * pkin(2) + pkin(3);
t54 = t133 * t135 + t73 * t69;
t80 = t1 * qJD(1) + t54 * qJD(2);
t2 = t82 + t137 / 0.2e1;
t53 = -t133 * t69 + t73 * t135;
t79 = -t2 * qJD(1) - t53 * qJD(2);
t78 = qJD(2) * t136 + t14 * qJD(4);
t47 = t54 * qJD(4);
t46 = t53 * qJD(4);
t23 = t134 / 0.2e1 + t76 * pkin(2);
t4 = 0.2e1 * t83 - t138;
t3 = -t145 + 0.2e1 * t82;
t5 = [0, 0, 0, t74 * t98, t67 * qJD(2), 0, 0, 0, -pkin(1) * t99, -pkin(1) * t98, t29 * qJD(2), t30 * qJD(2), t39 * qJD(3), t7 * qJD(2) + t12 * qJD(3), t11 * t136, (qJD(2) + qJD(4)) * t139, 0, 0, 0, t8 * qJD(2) + t44 * t104, t9 * qJD(2) - t44 * t107; 0, 0, 0, t91, t100, t98, -t99, 0, -pkin(5) * t98 - t97, pkin(5) * t99 - t96, -qJD(2) * t43 + t112, -qJD(2) * t85 + t111, (t120 * t57 - t59 * t72) * t126, t123 + (-t120 * t43 + t72 * t85) * t126 + t23 * qJD(3), -t94, t143, t11, -t78, 0, t146 * qJD(2) + t4 * qJD(4) + t122, t147 * qJD(2) + t3 * qJD(4) + t121; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t103, t23 * qJD(2) + t116, 0, 0, 0, 0, 0, t109, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t95, t143, t11, -t14 * qJD(2) - t104, 0, t4 * qJD(2) + t33 * qJD(3) + t146 * qJD(4) + t92, t3 * qJD(2) + t147 * qJD(4) - t93; 0, 0, 0, -t91, -t100, 0, 0, 0, t97, t96, -t59 * qJD(3) - t112, t57 * qJD(3) - t111, 0, t24 * qJD(3) - t123, t94, -t143, 0, -t109, 0, -qJD(3) * t136 - t1 * qJD(4) - t122, t2 * qJD(4) - t121 + t140; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t47, t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t101, t102, 0, t113, 0, 0, 0, 0, 0, -t108, t141; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t110, 0, -t47 - t80, t46 - t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59 * qJD(2), -t57 * qJD(2), -t103, -t24 * qJD(2) - t116, 0, 0, 0, 0, 0, t78, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t101, -t102, 0, -t113, 0, 0, 0, 0, 0, t108, -t141; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t114, -t141; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t95, -t143, 0, t33 * qJD(2), 0, t1 * qJD(2) - t14 * qJD(3) - t92, -t2 * qJD(2) + t140 + t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t110, 0, t80, t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t114, t141; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t5;
