% Calculate inertial parameters regressor of coriolis joint torque vector for
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
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 09:13
% Revision: 008671b0a00594318b890887636eaaff83cd5e2f (2021-12-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPPPR1_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR1_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 09:12:45
% EndTime: 2022-01-20 09:12:48
% DurationCPUTime: 0.74s
% Computational Cost: add. (1157->127), mult. (2884->204), div. (0->0), fcn. (2045->8), ass. (0->85)
t70 = sin(qJ(5));
t64 = sin(pkin(9));
t65 = sin(pkin(8));
t93 = qJD(1) * t65;
t83 = t64 * t93;
t67 = cos(pkin(9));
t71 = cos(qJ(5));
t99 = t71 * t67;
t85 = t65 * t99;
t30 = qJD(1) * t85 - t70 * t83;
t68 = cos(pkin(8));
t59 = t68 * qJD(1) - qJD(5);
t90 = qJD(4) * t65;
t91 = qJD(3) * t68;
t45 = -t64 * t91 - t67 * t90;
t36 = qJD(1) * t45;
t46 = -t64 * t90 + t67 * t91;
t37 = qJD(1) * t46;
t49 = -cos(pkin(7)) * pkin(1) - t68 * pkin(3) - t65 * qJ(4) - pkin(2);
t23 = qJD(1) * t49 + qJD(3);
t60 = sin(pkin(7)) * pkin(1) + qJ(3);
t58 = qJD(1) * t60;
t39 = t65 * qJD(2) + t68 * t58;
t12 = t64 * t23 + t67 * t39;
t10 = -pkin(6) * t83 + t12;
t11 = t67 * t23 - t64 * t39;
t108 = pkin(6) * t65;
t87 = t67 * t108;
t9 = (-t68 * pkin(4) - t87) * qJD(1) + t11;
t79 = t70 * t10 - t71 * t9;
t1 = -qJD(5) * t79 + t70 * t36 + t71 * t37;
t109 = t30 ^ 2;
t52 = t71 * t64 + t70 * t67;
t40 = t52 * t65;
t33 = qJD(5) * t40;
t21 = qJD(1) * t33;
t107 = t21 * t68;
t74 = qJD(1) * t52;
t27 = t65 * t74;
t106 = t27 * t59;
t105 = t30 * t27;
t104 = t30 * t59;
t103 = t33 * t59;
t102 = t60 * t68;
t22 = t30 * qJD(5);
t101 = t68 * t22;
t100 = t70 * t64;
t51 = t99 - t100;
t41 = t51 * t65;
t98 = -t41 * t22 + t33 * t27;
t97 = t52 * qJD(5) - t68 * t74;
t96 = t59 * t51;
t95 = t67 * t102 + t64 * t49;
t62 = t65 ^ 2;
t63 = t68 ^ 2;
t94 = t62 + t63;
t92 = qJD(3) * t65;
t88 = qJD(1) * qJD(3);
t72 = qJD(1) ^ 2;
t86 = t65 * t68 * t72;
t84 = 0.2e1 * qJD(3) * t62;
t82 = t94 * t72;
t61 = t65 * t88;
t38 = t68 * qJD(2) - t65 * t58;
t6 = t71 * t10 + t70 * t9;
t43 = t67 * t49;
t13 = -t87 + t43 + (-t60 * t64 - pkin(4)) * t68;
t16 = -t64 * t108 + t95;
t7 = t71 * t13 - t70 * t16;
t8 = t70 * t13 + t71 * t16;
t34 = (-t100 * t65 + t85) * qJD(5);
t78 = -t40 * t21 + t30 * t34;
t77 = -t36 * t67 - t37 * t64;
t75 = t38 * t65 - t39 * t68;
t35 = qJD(4) - t38;
t73 = qJD(5) * t27;
t2 = -qJD(5) * t6 + t71 * t36 - t70 * t37;
t54 = t62 * t60 * t88;
t44 = (pkin(4) * t64 + t60) * t65;
t26 = t27 ^ 2;
t18 = pkin(4) * t83 + t35;
t17 = t34 * t59;
t4 = -qJD(5) * t8 + t71 * t45 - t70 * t46;
t3 = qJD(5) * t7 + t70 * t45 + t71 * t46;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t94 * t88, t54 + (t63 * t58 - t75) * qJD(3), 0, 0, 0, 0, 0, 0, -t36 * t68 + (-t45 * t68 + t64 * t84) * qJD(1), t37 * t68 + (t46 * t68 + t67 * t84) * qJD(1), ((-t45 * t67 - t46 * t64) * qJD(1) + t77) * t65, t37 * t95 + t12 * t46 + t36 * (-t64 * t102 + t43) + t11 * t45 + t54 + t35 * t92, -t21 * t41 - t30 * t33, -t78 + t98, t103 + t107, t22 * t40 + t27 * t34, t17 + t101, 0, t18 * t34 - t2 * t68 + t44 * t22 - t4 * t59 + (qJD(1) * t40 + t27) * t92, t1 * t68 - t18 * t33 - t44 * t21 + t3 * t59 + (qJD(1) * t41 + t30) * t92, -t1 * t40 - t2 * t41 + t7 * t21 - t8 * t22 - t3 * t27 - t4 * t30 - t33 * t79 - t6 * t34, t1 * t8 + t2 * t7 + t6 * t3 - t79 * t4 + (qJD(1) * t44 + t18) * t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (-t36 * t64 + t37 * t67 - t68 * t88) * t65, 0, 0, 0, 0, 0, 0, t17 - t101, -t103 + t107, t78 + t98, t1 * t41 - t2 * t40 - t6 * t33 + t34 * t79 - t68 * t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t82, t75 * qJD(1), 0, 0, 0, 0, 0, 0, -t64 * t82, -t67 * t82, 0, (-t35 * t65 + (t11 * t64 - t12 * t67) * t68) * qJD(1) - t77, 0, 0, 0, 0, 0, 0, -t27 * t93 + t97 * t59, -t30 * t93 - t96 * t59, t51 * t21 - t52 * t22 + t96 * t27 + t97 * t30, t1 * t52 - t18 * t93 + t2 * t51 - t96 * t6 + t79 * t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t67 * t86, t64 * t86, (-t64 ^ 2 - t67 ^ 2) * t72 * t62, t61 + (t11 * t67 + t12 * t64) * t93, 0, 0, 0, 0, 0, 0, t22 - t104, -t73 + t106, -t26 - t109, t6 * t27 - t30 * t79 + t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t105, -t26 + t109, -t73 - t106, -t105, -t22 - t104, 0, -t18 * t30 - t6 * t59 + t2, t18 * t27 + t59 * t79 - t1, 0, 0;];
tauc_reg = t5;
