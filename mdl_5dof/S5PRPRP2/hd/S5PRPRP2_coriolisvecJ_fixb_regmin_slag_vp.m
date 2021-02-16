% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5PRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% 
% Output:
% tauc_reg [5x18]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 15:05
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRPRP2_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP2_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP2_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP2_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 15:04:51
% EndTime: 2021-01-15 15:04:54
% DurationCPUTime: 0.52s
% Computational Cost: add. (500->110), mult. (1288->184), div. (0->0), fcn. (786->4), ass. (0->84)
t41 = cos(pkin(8));
t42 = sin(qJ(4));
t43 = cos(qJ(4));
t40 = sin(pkin(8));
t27 = -t41 * pkin(3) - t40 * pkin(6) - pkin(2);
t84 = qJ(5) * t40;
t57 = -t27 + t84;
t95 = -t43 * t41 * qJ(3) + t57 * t42;
t78 = qJD(4) * t43;
t61 = qJ(5) * t78;
t18 = t27 * qJD(2) + qJD(3);
t74 = qJ(3) * qJD(2);
t25 = t40 * qJD(1) + t41 * t74;
t73 = qJD(2) * qJD(3);
t60 = t41 * t73;
t79 = qJD(4) * t42;
t70 = -t18 * t78 + t25 * t79 - t43 * t60;
t77 = qJD(5) * t42;
t82 = qJD(2) * t40;
t1 = (-t61 - t77) * t82 - t70;
t66 = t40 * t79;
t54 = qJD(2) * t66;
t28 = qJ(5) * t54;
t81 = qJD(3) * t42;
t64 = t41 * t81;
t89 = t40 * t43;
t46 = -qJD(5) * t89 - t64;
t49 = -t42 * t18 - t43 * t25;
t47 = t49 * qJD(4);
t2 = t46 * qJD(2) + t28 + t47;
t76 = t41 * qJD(2);
t33 = -qJD(4) + t76;
t58 = t43 * t18 - t42 * t25;
t62 = qJ(5) * t82;
t6 = -t43 * t62 + t58;
t3 = -t33 * pkin(4) + t6;
t7 = -t42 * t62 - t49;
t50 = t3 * t42 - t43 * t7;
t94 = -t50 * qJD(4) + t1 * t42 + t2 * t43;
t37 = t41 ^ 2;
t36 = t40 ^ 2;
t93 = 0.2e1 * t36;
t92 = t3 - t6;
t91 = t33 * t41;
t44 = qJD(2) ^ 2;
t90 = t36 * t44;
t80 = qJD(3) * t43;
t88 = t27 * t78 + t41 * t80;
t72 = qJD(2) * qJD(4);
t59 = t43 * t72;
t53 = t40 * t59;
t17 = pkin(4) * t53 + t40 * t73;
t87 = t36 + t37;
t38 = t42 ^ 2;
t39 = t43 ^ 2;
t86 = t38 - t39;
t85 = qJ(3) * t42;
t83 = qJD(2) * t36;
t23 = (pkin(4) * t42 + qJ(3)) * t40;
t35 = t41 * qJD(1);
t14 = qJD(2) * t23 + qJD(5) - t35;
t75 = qJD(5) + t14;
t71 = t42 * t90;
t68 = t42 * t82;
t67 = t43 * t82;
t65 = t40 * t78;
t63 = qJ(3) * t79;
t56 = t87 * qJD(2);
t55 = t33 * t66;
t51 = t3 * t43 + t42 * t7;
t24 = t40 * t74 - t35;
t48 = t24 * t40 + t25 * t41;
t45 = -t33 ^ 2 - t90;
t26 = t41 * t54;
t22 = (pkin(4) * t78 + qJD(3)) * t40;
t21 = t33 * t67;
t13 = (t33 - t76) * t65;
t12 = t26 - t55;
t10 = -t57 * t43 + (-pkin(4) - t85) * t41;
t9 = t45 * t43;
t8 = t45 * t42;
t5 = t95 * qJD(4) + t46;
t4 = -t40 * t77 + (-t41 * t85 - t43 * t84) * qJD(4) + t88;
t11 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, t12, t13, t12, 0, -t17 * t41 + (-t51 * qJD(4) + t1 * t43 - t2 * t42) * t40; 0, 0, 0, 0, 0, 0.2e1 * qJD(3) * t56, (qJ(3) * t56 + t48) * qJD(3), -0.2e1 * t36 * t42 * t59, t86 * t72 * t93, t26 + t55, (t33 + t76) * t65, 0, (t91 + (t93 + t37) * qJD(2)) * t81 + ((t18 * t41 + t27 * t33) * t42 + ((t83 + t91) * qJ(3) + t48) * t43) * qJD(4), (-t41 * t63 + t88) * t33 - t70 * t41 - t24 * t66 + (-t63 + 0.2e1 * t80) * t83, -t2 * t41 - t5 * t33 + (t14 * t78 + t17 * t42 + (t22 * t42 + t23 * t78) * qJD(2)) * t40, t1 * t41 + t4 * t33 + (-t14 * t79 + t17 * t43 + (t22 * t43 - t23 * t79) * qJD(2)) * t40, ((-t4 * t42 - t43 * t5 + (t10 * t42 + t43 * t95) * qJD(4)) * qJD(2) - t94) * t40, -t1 * t95 + t2 * t10 + t14 * t22 + t17 * t23 + t3 * t5 + t7 * t4; 0, 0, 0, 0, 0, -t87 * t44, -t48 * qJD(2), 0, 0, 0, 0, 0, t8, t9, t8, t9, 0, (-t14 * t40 + t50 * t41) * qJD(2) + t94; 0, 0, 0, 0, 0, 0, 0, t43 * t71, -t86 * t90, (-qJD(4) - t33) * t68, -t21 - t53, 0, t49 * t33 + t47 + (-t24 * t89 - t64) * qJD(2), t24 * t68 - t58 * t33 + t70, -t7 * t33 + t28 + (-qJD(4) * t18 - t60) * t42 + (-pkin(4) * t71 - qJD(4) * t25 - t75 * t82) * t43, -t39 * pkin(4) * t90 - t6 * t33 + (t75 * t42 + t61) * t82 + t70, (pkin(4) * qJD(4) - t92) * t68, t92 * t7 + (-t14 * t67 + t2) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21 + t53, (-qJD(4) + t33) * t68, (-t38 - t39) * t90, t51 * t82 + t17;];
tauc_reg = t11;
