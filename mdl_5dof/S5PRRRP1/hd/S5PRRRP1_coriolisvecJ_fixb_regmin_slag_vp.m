% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5PRRRP1
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
% tauc_reg [5x18]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 16:15
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRRRP1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 16:14:40
% EndTime: 2021-01-15 16:14:43
% DurationCPUTime: 0.44s
% Computational Cost: add. (552->127), mult. (991->169), div. (0->0), fcn. (474->4), ass. (0->87)
t46 = sin(qJ(4));
t48 = cos(qJ(4));
t43 = qJD(2) + qJD(3);
t47 = sin(qJ(3));
t80 = pkin(2) * qJD(2);
t68 = t47 * t80;
t25 = t43 * pkin(7) + t68;
t78 = qJ(5) * t43;
t62 = t25 + t78;
t53 = t62 * t48;
t8 = t46 * qJD(1) + t53;
t39 = t48 * qJD(1);
t7 = -t62 * t46 + t39;
t77 = qJD(4) * pkin(4);
t4 = t7 + t77;
t95 = t4 - t7;
t44 = t46 ^ 2;
t94 = pkin(4) * t44;
t93 = t48 * pkin(4);
t49 = cos(qJ(3));
t92 = t49 * pkin(2);
t67 = t49 * t80;
t26 = -t43 * pkin(3) - t67;
t91 = t26 * t43;
t42 = t43 ^ 2;
t90 = t42 * t48;
t89 = t43 * t46;
t88 = t43 * t48;
t50 = qJD(4) ^ 2;
t87 = t50 * t46;
t40 = t50 * t48;
t86 = -qJ(5) - pkin(7);
t37 = -pkin(3) - t93;
t11 = t37 * t43 + qJD(5) - t67;
t79 = pkin(2) * qJD(3);
t63 = qJD(2) * t79;
t34 = t47 * t63;
t73 = t46 * qJD(4);
t66 = t43 * t73;
t16 = pkin(4) * t66 + t34;
t72 = t48 * qJD(4);
t85 = t11 * t72 + t16 * t46;
t84 = t26 * t72 + t46 * t34;
t57 = qJD(4) * t67;
t58 = t43 * t68;
t83 = t46 * t57 + t48 * t58;
t45 = t48 ^ 2;
t82 = -t44 - t45;
t81 = t44 - t45;
t35 = t47 * pkin(2) + pkin(7);
t76 = -qJ(5) - t35;
t75 = qJD(5) * t43;
t71 = -qJD(2) - t43;
t70 = t49 * t79;
t69 = t47 * t79;
t9 = t11 * t73;
t65 = t43 * t72;
t19 = t25 * t73;
t56 = t49 * t63;
t52 = qJD(4) * qJD(1) + t56;
t2 = -qJ(5) * t66 - t19 + (t52 + t75) * t48;
t3 = (-t56 - t75) * t46 - t8 * qJD(4);
t64 = t2 * t48 - t3 * t46;
t61 = qJD(4) * t86;
t60 = 0.2e1 * t65;
t59 = qJD(4) * t76;
t55 = -t4 * t48 - t46 * t8;
t54 = t4 * t46 - t48 * t8;
t51 = (-qJD(5) - t11) * t43 - t52;
t41 = t48 * qJ(5);
t38 = t48 * qJD(5);
t36 = -pkin(3) - t92;
t32 = t48 * pkin(7) + t41;
t31 = t86 * t46;
t29 = t48 * t57;
t27 = t37 - t92;
t24 = t46 * t60;
t22 = pkin(4) * t73 + t69;
t21 = t48 * t35 + t41;
t20 = t76 * t46;
t17 = t26 * t73;
t15 = -t46 * qJD(5) + t48 * t61;
t14 = t46 * t61 + t38;
t13 = -0.2e1 * t81 * t43 * qJD(4);
t6 = (-qJD(5) - t70) * t46 + t48 * t59;
t5 = t46 * t59 + t48 * t70 + t38;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t87, -t40, -t87, -t40, 0, -t54 * qJD(4) + t2 * t46 + t3 * t48; 0, 0, 0, 0, 0, -t43 * t69 - t34, t71 * t70, t24, t13, t40, -t87, 0, t36 * t66 - t35 * t40 + t17 + (t71 * t48 * t47 - t49 * t73) * t79, t36 * t65 + t35 * t87 + (t47 * t89 - t49 * t72) * t79 + t84, t9 + (-t22 * t43 - t16) * t48 + (t27 * t89 + t6) * qJD(4), t22 * t89 + (t27 * t88 - t5) * qJD(4) + t85, (-t46 * t6 + t48 * t5) * t43 + ((-t20 * t48 - t21 * t46) * t43 + t55) * qJD(4) + t64, t11 * t22 + t16 * t27 + t2 * t21 + t3 * t20 + t4 * t6 + t8 * t5; 0, 0, 0, 0, 0, -t34 + t58, (-qJD(3) + t43) * t67, t24, t13, t40, -t87, 0, -pkin(3) * t66 + t17 + (-pkin(7) * t50 - t34) * t48 + t83, pkin(7) * t87 + t29 + (-pkin(3) * t72 - t46 * t68) * t43 + t84, -t16 * t48 + t9 + (t15 + (t37 - t93) * t89) * qJD(4) + t83, -t46 * t58 + t29 + (-t14 + (t37 * t48 + t94) * t43) * qJD(4) + t85, t55 * qJD(4) + (t14 * t48 - t15 * t46 + (-t31 * t48 - t32 * t46) * qJD(4) + t82 * t67) * t43 + t64, pkin(4) * t9 + t8 * t14 + t4 * t15 + t16 * t37 + t2 * t32 + t3 * t31 + (-t11 * t47 + t54 * t49) * t80; 0, 0, 0, 0, 0, 0, 0, -t46 * t90, t81 * t42, 0, 0, 0, (-t56 - t91) * t46, t19 + (-t46 * t25 + t39) * qJD(4) + (-t52 - t91) * t48, (t8 - t53) * qJD(4) + (pkin(4) * t90 + t51) * t46, -t42 * t94 + t19 + (t46 * t78 + t7) * qJD(4) + t51 * t48, (-t77 + t95) * t88, t95 * t8 + (-t11 * t89 + t3) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t66, t60, t82 * t42, t54 * t43 + t16;];
tauc_reg = t1;
