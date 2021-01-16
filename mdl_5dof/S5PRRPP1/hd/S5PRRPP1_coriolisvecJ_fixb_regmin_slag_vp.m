% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5PRRPP1
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
% tauc_reg [5x19]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 15:23
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRRPP1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 15:22:40
% EndTime: 2021-01-15 15:22:42
% DurationCPUTime: 0.44s
% Computational Cost: add. (822->142), mult. (2101->194), div. (0->0), fcn. (1331->4), ass. (0->85)
t85 = (qJD(2) * qJD(3));
t105 = -2 * t85;
t62 = sin(pkin(8));
t63 = cos(pkin(8));
t64 = sin(qJ(3));
t65 = cos(qJ(3));
t47 = t62 * t65 + t63 * t64;
t92 = qJD(2) * t47;
t39 = t92 ^ 2;
t98 = t63 * t65;
t81 = qJD(2) * t98;
t90 = qJD(2) * t64;
t40 = t62 * t90 - t81;
t104 = -t40 ^ 2 - t39;
t103 = pkin(3) * t64;
t94 = -qJ(4) - pkin(6);
t50 = t94 * t65;
t80 = t94 * t64;
t26 = -t62 * t50 - t63 * t80;
t77 = qJD(3) * t94;
t36 = t65 * qJD(4) + t64 * t77;
t84 = qJD(3) * qJD(1);
t25 = t36 * qJD(2) + t65 * t84;
t71 = -t64 * qJD(4) + t65 * t77;
t68 = t71 * qJD(2) - t64 * t84;
t3 = t62 * t25 - t63 * t68;
t102 = t3 * t26;
t46 = t62 * t64 - t98;
t101 = t3 * t46;
t58 = -t65 * pkin(3) - pkin(2);
t91 = qJD(2) * t58;
t49 = qJD(4) + t91;
t10 = t40 * pkin(4) - qJ(5) * t92 + t49;
t100 = t10 * t92;
t37 = t64 * qJD(1) - qJD(2) * t50;
t99 = t62 * t37;
t29 = t63 * t37;
t67 = qJD(2) ^ 2;
t97 = t65 * t67;
t66 = qJD(3) ^ 2;
t96 = t66 * t64;
t95 = t66 * t65;
t4 = t63 * t25 + t62 * t68;
t35 = t65 * qJD(1) + qJD(2) * t80;
t32 = qJD(3) * pkin(3) + t35;
t12 = t62 * t32 + t29;
t93 = t64 ^ 2 - t65 ^ 2;
t42 = t47 * qJD(3);
t89 = t42 * qJD(3);
t87 = t64 * qJD(3);
t45 = qJD(3) * t98 - t62 * t87;
t88 = t45 * qJD(3);
t16 = t63 * t35 - t99;
t86 = qJD(5) - t16;
t83 = pkin(3) * t87;
t82 = pkin(3) * t90;
t79 = t64 * t85;
t78 = t65 * t85;
t76 = 0.2e1 * t92;
t75 = pkin(2) * t105;
t11 = t63 * t32 - t99;
t33 = qJD(2) * t42;
t51 = t62 * t79;
t34 = t63 * t78 - t51;
t54 = pkin(3) * t79;
t74 = t33 * pkin(4) - t34 * qJ(5) + t54;
t14 = t62 * t35 + t29;
t73 = t14 * qJD(3) - t3;
t72 = -t47 * t33 + t46 * t34 - t45 * t40 + t42 * t92;
t15 = t62 * t36 - t63 * t71;
t17 = t63 * t36 + t62 * t71;
t27 = -t63 * t50 + t62 * t80;
t70 = t15 * t92 - t17 * t40 + t26 * t34 - t27 * t33 + t3 * t47;
t69 = t76 * qJD(3);
t57 = -t63 * pkin(3) - pkin(4);
t55 = t62 * pkin(3) + qJ(5);
t19 = t46 * pkin(4) - t47 * qJ(5) + t58;
t18 = -t51 + (-t40 + t81) * qJD(3);
t13 = pkin(4) * t92 + t40 * qJ(5) + t82;
t8 = qJD(3) * qJ(5) + t12;
t7 = -qJD(3) * pkin(4) + qJD(5) - t11;
t6 = t42 * pkin(4) - t45 * qJ(5) - t47 * qJD(5) + t83;
t2 = qJD(3) * qJD(5) + t4;
t1 = -qJD(5) * t92 + t74;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -t96, -t95, -t89, -t88, t72, -t11 * t42 + t12 * t45 + t4 * t47 + t101, -t89, t72, t88, t2 * t47 + t7 * t42 + t8 * t45 + t101; 0, 0, 0, 0, 0.2e1 * t64 * t78, t93 * t105, t95, -t96, 0, -pkin(6) * t95 + t64 * t75, pkin(6) * t96 + t65 * t75, t58 * t33 + t49 * t42 + (-t15 + (qJD(2) * t46 + t40) * t103) * qJD(3), t58 * t34 + t49 * t45 + (t76 * t103 - t17) * qJD(3), -t11 * t45 - t12 * t42 - t4 * t46 + t70, -t11 * t15 + t12 * t17 + t102 + t4 * t27 + (t49 + t91) * t83, -t15 * qJD(3) + t1 * t46 + t10 * t42 + t19 * t33 + t6 * t40, -t2 * t46 - t8 * t42 + t7 * t45 + t70, t17 * qJD(3) - t1 * t47 - t10 * t45 - t19 * t34 - t6 * t92, t1 * t19 + t10 * t6 + t7 * t15 + t8 * t17 + t2 * t27 + t102; 0, 0, 0, 0, -t64 * t97, t93 * t67, 0, 0, 0, t67 * pkin(2) * t64, pkin(2) * t97, -t40 * t82 - t49 * t92 + t73, t16 * qJD(3) + t49 * t40 - t82 * t92 - t4, (t12 - t14) * t92 + (-t11 + t16) * t40 + (-t33 * t62 - t34 * t63) * pkin(3), t11 * t14 - t12 * t16 + (-t3 * t63 + t4 * t62 - t49 * t90) * pkin(3), -t13 * t40 - t100 + t73, -t55 * t33 + t57 * t34 + (-t14 + t8) * t92 + (t7 - t86) * t40, -t10 * t40 + t13 * t92 + (0.2e1 * qJD(5) - t16) * qJD(3) + t4, -t10 * t13 - t7 * t14 + t2 * t55 + t3 * t57 + t86 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, t18, t104, t11 * t92 + t12 * t40 + t54, t69, t104, -t18, t8 * t40 + (-qJD(5) - t7) * t92 + t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t92 * t40, -t51 + (t40 + t81) * qJD(3), -t39 - t66, -t8 * qJD(3) + t100 + t3;];
tauc_reg = t5;
