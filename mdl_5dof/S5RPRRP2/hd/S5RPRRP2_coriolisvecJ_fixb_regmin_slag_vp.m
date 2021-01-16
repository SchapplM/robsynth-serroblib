% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RPRRP2
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
% tauc_reg [5x18]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 12:36
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRRP2_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP2_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP2_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP2_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 12:36:27
% EndTime: 2021-01-15 12:36:30
% DurationCPUTime: 0.47s
% Computational Cost: add. (820->135), mult. (1642->171), div. (0->0), fcn. (876->6), ass. (0->92)
t60 = sin(qJ(4));
t62 = cos(qJ(4));
t48 = cos(pkin(8)) * pkin(1) + pkin(2);
t42 = t48 * qJD(1);
t61 = sin(qJ(3));
t63 = cos(qJ(3));
t106 = pkin(1) * sin(pkin(8));
t81 = qJD(1) * t106;
t30 = t61 * t42 + t63 * t81;
t55 = qJD(1) + qJD(3);
t20 = t55 * pkin(7) + t30;
t90 = qJ(5) * t55;
t77 = t20 + t90;
t68 = t77 * t62;
t8 = t60 * qJD(2) + t68;
t108 = -t61 * t106 + t63 * t48;
t51 = t62 * qJD(2);
t7 = -t77 * t60 + t51;
t89 = qJD(4) * pkin(4);
t6 = t7 + t89;
t107 = t6 - t7;
t56 = t60 ^ 2;
t105 = pkin(4) * t56;
t104 = t62 * pkin(4);
t29 = t63 * t42 - t61 * t81;
t19 = -t55 * pkin(3) - t29;
t103 = t19 * t55;
t66 = t63 * t106 + t61 * t48;
t32 = t66 * qJD(3);
t102 = t32 * t55;
t54 = t55 ^ 2;
t101 = t54 * t62;
t100 = t55 * t60;
t99 = t55 * t62;
t64 = qJD(4) ^ 2;
t97 = t64 * t60;
t96 = -qJ(5) - pkin(7);
t49 = -pkin(3) - t104;
t12 = t49 * t55 + qJD(5) - t29;
t72 = qJD(3) * t81;
t87 = qJD(3) * t42;
t27 = t61 * t87 + t63 * t72;
t84 = t60 * qJD(4);
t79 = t55 * t84;
t13 = pkin(4) * t79 + t27;
t83 = t62 * qJD(4);
t95 = t12 * t83 + t13 * t60;
t94 = t19 * t83 + t27 * t60;
t93 = t29 * t84 + t30 * t99;
t57 = t62 ^ 2;
t92 = -t56 - t57;
t91 = t56 - t57;
t35 = pkin(7) + t66;
t88 = -qJ(5) - t35;
t86 = qJD(5) * t55;
t82 = 0.2e1 * qJD(4) * t55;
t80 = pkin(4) * t84;
t16 = t20 * t84;
t26 = -t61 * t72 + t63 * t87;
t73 = qJD(4) * qJD(2) + t26;
t2 = -qJ(5) * t79 - t16 + (t73 + t86) * t62;
t3 = (-t26 - t86) * t60 - t8 * qJD(4);
t78 = t2 * t62 - t3 * t60;
t76 = qJD(4) * t96;
t75 = t62 * t82;
t74 = qJD(4) * t88;
t34 = -pkin(3) - t108;
t71 = -t6 * t62 - t60 * t8;
t70 = t6 * t60 - t62 * t8;
t69 = t35 * t64 + t102;
t31 = t108 * qJD(3);
t67 = qJD(4) * (t34 * t55 - t31);
t65 = (-qJD(5) - t12) * t55 - t73;
t53 = t62 * qJ(5);
t52 = t64 * t62;
t50 = t62 * qJD(5);
t44 = t62 * pkin(7) + t53;
t43 = t96 * t60;
t39 = t60 * t75;
t37 = -t60 * qJD(5) + t62 * t76;
t36 = t60 * t76 + t50;
t33 = t91 * t82;
t28 = t34 - t104;
t25 = t32 + t80;
t24 = t62 * t35 + t53;
t23 = t88 * t60;
t22 = t29 * t83;
t14 = t19 * t84;
t9 = t12 * t84;
t5 = (-qJD(5) - t31) * t60 + t62 * t74;
t4 = t62 * t31 + t60 * t74 + t50;
t1 = [0, 0, 0, 0, 0, -t27 - t102, -t31 * t55 - t26, t39, -t33, t52, -t97, 0, t14 + t60 * t67 + (-t27 - t69) * t62, t69 * t60 + t62 * t67 + t94, t9 + (-t25 * t55 - t13) * t62 + (t28 * t100 + t5) * qJD(4), t25 * t100 + (t28 * t99 - t4) * qJD(4) + t95, (t4 * t62 - t5 * t60) * t55 + ((-t23 * t62 - t24 * t60) * t55 + t71) * qJD(4) + t78, t12 * t25 + t13 * t28 + t2 * t24 + t3 * t23 + t8 * t4 + t6 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t97, -t52, -t97, -t52, 0, -t70 * qJD(4) + t2 * t60 + t3 * t62; 0, 0, 0, 0, 0, t30 * t55 - t27, t29 * t55 - t26, t39, -t33, t52, -t97, 0, -pkin(3) * t79 + t14 + (-pkin(7) * t64 - t27) * t62 + t93, pkin(7) * t97 + t22 + (-pkin(3) * t83 - t30 * t60) * t55 + t94, -t13 * t62 + t9 + (t37 + (t49 - t104) * t100) * qJD(4) + t93, -t30 * t100 + t22 + (-t36 + (t49 * t62 + t105) * t55) * qJD(4) + t95, t71 * qJD(4) + (t36 * t62 - t37 * t60 + t92 * t29 + (-t43 * t62 - t44 * t60) * qJD(4)) * t55 + t78, t13 * t49 + t2 * t44 + t3 * t43 + t8 * t36 + t6 * t37 + t70 * t29 + (-t30 + t80) * t12; 0, 0, 0, 0, 0, 0, 0, -t60 * t101, t91 * t54, 0, 0, 0, (-t26 - t103) * t60, t16 + (-t60 * t20 + t51) * qJD(4) + (-t73 - t103) * t62, (t8 - t68) * qJD(4) + (pkin(4) * t101 + t65) * t60, -t54 * t105 + t16 + (t60 * t90 + t7) * qJD(4) + t65 * t62, (-t89 + t107) * t99, t107 * t8 + (-t12 * t100 + t3) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t79, t75, t92 * t54, t70 * t55 + t13;];
tauc_reg = t1;
