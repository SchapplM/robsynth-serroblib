% Calculate minimal parameter regressor of coriolis joint torque vector for
% S4RRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,theta3]';
% 
% Output:
% tauc_reg [4x18]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 10:36
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RRPP3_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP3_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP3_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP3_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 10:36:04
% EndTime: 2021-01-15 10:36:07
% DurationCPUTime: 0.42s
% Computational Cost: add. (612->124), mult. (1676->176), div. (0->0), fcn. (1035->4), ass. (0->80)
t80 = (qJD(1) * qJD(2));
t97 = -2 * t80;
t60 = sin(pkin(6));
t61 = cos(pkin(6));
t62 = sin(qJ(2));
t63 = cos(qJ(2));
t44 = t60 * t63 + t61 * t62;
t85 = qJD(1) * t44;
t33 = t85 ^ 2;
t91 = t61 * t63;
t77 = qJD(1) * t91;
t83 = qJD(1) * t62;
t34 = t60 * t83 - t77;
t96 = -t34 ^ 2 - t33;
t95 = pkin(2) * t62;
t87 = -qJ(3) - pkin(5);
t49 = t87 * t63;
t76 = t87 * t62;
t20 = -t60 * t49 - t61 * t76;
t73 = qJD(2) * t87;
t31 = t63 * qJD(3) + t62 * t73;
t25 = t31 * qJD(1);
t68 = -t62 * qJD(3) + t63 * t73;
t26 = t68 * qJD(1);
t5 = t60 * t25 - t61 * t26;
t94 = t5 * t20;
t57 = -t63 * pkin(2) - pkin(1);
t84 = qJD(1) * t57;
t48 = qJD(3) + t84;
t8 = t34 * pkin(3) - qJ(4) * t85 + t48;
t93 = t8 * t85;
t47 = qJD(1) * t49;
t92 = t60 * t47;
t40 = t61 * t47;
t65 = qJD(1) ^ 2;
t90 = t63 * t65;
t64 = qJD(2) ^ 2;
t89 = t64 * t62;
t88 = t64 * t63;
t6 = t61 * t25 + t60 * t26;
t46 = qJD(1) * t76;
t42 = qJD(2) * pkin(2) + t46;
t17 = t60 * t42 - t40;
t86 = t62 ^ 2 - t63 ^ 2;
t82 = qJD(2) * t62;
t19 = t61 * t46 + t92;
t81 = qJD(4) - t19;
t79 = pkin(2) * t82;
t78 = pkin(2) * t83;
t75 = t62 * t80;
t74 = t63 * t80;
t72 = 0.2e1 * t85;
t71 = pkin(1) * t97;
t16 = t61 * t42 + t92;
t36 = t44 * qJD(2);
t28 = qJD(1) * t36;
t50 = t60 * t75;
t29 = t61 * t74 - t50;
t53 = pkin(2) * t75;
t70 = t28 * pkin(3) - t29 * qJ(4) + t53;
t18 = t60 * t46 - t40;
t69 = t18 * qJD(2) - t5;
t10 = t60 * t31 - t61 * t68;
t11 = t61 * t31 + t60 * t68;
t21 = -t61 * t49 + t60 * t76;
t67 = t10 * t85 - t11 * t34 + t20 * t29 - t21 * t28 + t5 * t44;
t66 = t72 * qJD(2);
t56 = -t61 * pkin(2) - pkin(3);
t54 = t60 * pkin(2) + qJ(4);
t43 = t60 * t62 - t91;
t39 = qJD(2) * t91 - t60 * t82;
t15 = t43 * pkin(3) - t44 * qJ(4) + t57;
t14 = -t50 + (-t34 + t77) * qJD(2);
t13 = qJD(2) * qJ(4) + t17;
t12 = -qJD(2) * pkin(3) + qJD(4) - t16;
t9 = pkin(3) * t85 + t34 * qJ(4) + t78;
t4 = qJD(2) * qJD(4) + t6;
t3 = t36 * pkin(3) - t39 * qJ(4) - t44 * qJD(4) + t79;
t1 = -qJD(4) * t85 + t70;
t2 = [0, 0, 0, 0.2e1 * t62 * t74, t86 * t97, t88, -t89, 0, -pkin(5) * t88 + t62 * t71, pkin(5) * t89 + t63 * t71, t57 * t28 + t48 * t36 + (-t10 + (qJD(1) * t43 + t34) * t95) * qJD(2), t57 * t29 + t48 * t39 + (t72 * t95 - t11) * qJD(2), -t16 * t39 - t17 * t36 - t6 * t43 + t67, -t16 * t10 + t17 * t11 + t94 + t6 * t21 + (t48 + t84) * t79, -t10 * qJD(2) + t1 * t43 + t15 * t28 + t3 * t34 + t8 * t36, t12 * t39 - t13 * t36 - t4 * t43 + t67, t11 * qJD(2) - t1 * t44 - t15 * t29 - t3 * t85 - t8 * t39, t1 * t15 + t12 * t10 + t13 * t11 + t4 * t21 + t8 * t3 + t94; 0, 0, 0, -t62 * t90, t86 * t65, 0, 0, 0, t65 * pkin(1) * t62, pkin(1) * t90, -t34 * t78 - t48 * t85 + t69, t19 * qJD(2) + t48 * t34 - t78 * t85 - t6, (t17 - t18) * t85 + (-t16 + t19) * t34 + (-t28 * t60 - t29 * t61) * pkin(2), t16 * t18 - t17 * t19 + (-t48 * t83 - t5 * t61 + t6 * t60) * pkin(2), -t9 * t34 + t69 - t93, -t54 * t28 + t56 * t29 + (t13 - t18) * t85 + (t12 - t81) * t34, -t8 * t34 + t9 * t85 + (0.2e1 * qJD(4) - t19) * qJD(2) + t6, -t12 * t18 + t81 * t13 + t4 * t54 + t5 * t56 - t8 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66, t14, t96, t16 * t85 + t17 * t34 + t53, t66, t96, -t14, t13 * t34 + (-qJD(4) - t12) * t85 + t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t85 * t34, -t50 + (t34 + t77) * qJD(2), -t33 - t64, -t13 * qJD(2) + t5 + t93;];
tauc_reg = t2;
