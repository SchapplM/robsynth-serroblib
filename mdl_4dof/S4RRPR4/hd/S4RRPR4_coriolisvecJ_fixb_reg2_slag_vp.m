% Calculate inertial parameters regressor of coriolis joint torque vector for
% S4RRPR4
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
% tauc_reg [4x(4*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RRPR4_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR4_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR4_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR4_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:02:33
% EndTime: 2019-12-31 17:02:35
% DurationCPUTime: 0.46s
% Computational Cost: add. (900->107), mult. (1660->157), div. (0->0), fcn. (1052->6), ass. (0->82)
t65 = qJD(1) + qJD(2);
t71 = cos(qJ(2));
t91 = pkin(1) * qJD(2);
t83 = qJD(1) * t91;
t44 = t65 * qJD(3) + t71 * t83;
t66 = sin(pkin(7));
t67 = cos(pkin(7));
t93 = t66 ^ 2 + t67 ^ 2;
t106 = t93 * t44;
t92 = pkin(1) * qJD(1);
t73 = -t71 * t92 + qJD(3);
t68 = sin(qJ(4));
t100 = t68 * t66;
t70 = cos(qJ(4));
t99 = t70 * t67;
t45 = -t99 + t100;
t105 = t45 * t44;
t46 = t70 * t66 + t68 * t67;
t32 = t46 * t65;
t69 = sin(qJ(2));
t104 = t67 * t69;
t103 = t32 ^ 2;
t102 = t71 * pkin(1);
t90 = t65 * t100;
t30 = -t65 * t99 + t90;
t101 = t32 * t30;
t51 = (-pkin(6) - qJ(3)) * t66;
t61 = t67 * pkin(6);
t52 = t67 * qJ(3) + t61;
t20 = t70 * t51 - t68 * t52;
t98 = t20 * qJD(4) - t73 * t45;
t21 = t68 * t51 + t70 * t52;
t97 = -t21 * qJD(4) - t73 * t46;
t60 = -t67 * pkin(3) - pkin(2);
t29 = t60 * t65 + t73;
t41 = t46 * qJD(4);
t57 = t69 * t83;
t96 = t29 * t41 + t45 * t57;
t85 = qJD(4) * t99;
t86 = qJD(4) * t100;
t40 = -t85 + t86;
t95 = -t29 * t40 + t46 * t57;
t89 = t69 * t91;
t88 = t69 * t92;
t48 = t65 * qJ(3) + t88;
t82 = pkin(6) * t65 + t48;
t22 = t82 * t66;
t23 = t82 * t67;
t12 = -t70 * t22 - t68 * t23;
t13 = -t68 * t22 + t70 * t23;
t4 = t12 * qJD(4) - t105;
t72 = t46 * t44;
t5 = -t13 * qJD(4) - t72;
t84 = t12 * t40 - t13 * t41 - t4 * t45 - t5 * t46;
t81 = t93 * t71;
t55 = t71 * t91 + qJD(3);
t79 = t93 * t55;
t78 = t93 * qJD(3);
t77 = t65 * t89;
t76 = t65 * t88;
t75 = (-qJD(2) + t65) * t92;
t74 = (-qJD(1) - t65) * t91;
t59 = t69 * pkin(1) + qJ(3);
t42 = (-pkin(6) - t59) * t66;
t43 = t67 * t59 + t61;
t18 = t70 * t42 - t68 * t43;
t19 = t68 * t42 + t70 * t43;
t53 = t66 * t57;
t50 = t60 - t102;
t49 = t65 * t85;
t47 = -t65 * pkin(2) + t73;
t37 = t41 * qJD(4);
t36 = t40 * qJD(4);
t28 = t30 ^ 2;
t27 = t65 * t41;
t26 = t65 * t86 - t49;
t11 = -t19 * qJD(4) - t46 * t55;
t10 = t18 * qJD(4) - t45 * t55;
t7 = t27 * t45 + t30 * t41;
t6 = -t26 * t46 - t32 * t40;
t1 = t26 * t45 - t46 * t27 + t40 * t30 - t32 * t41;
t2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t57 - t77, t71 * t74, 0, 0, 0, 0, 0, 0, 0, 0, t74 * t104, t66 * t77 + t53, t65 * t79 + t106, t48 * t79 + t59 * t106 + (t47 + (-pkin(2) - t102) * qJD(1)) * t89, t6, t1, -t36, t7, -t37, 0, t11 * qJD(4) + t50 * t27 + t30 * t89 + t96, -t10 * qJD(4) - t50 * t26 + t32 * t89 + t95, -t10 * t30 - t11 * t32 + t18 * t26 - t19 * t27 + t84, t13 * t10 + t12 * t11 + t5 * t18 + t4 * t19 + (qJD(1) * t50 + t29) * t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t57 + t76, t71 * t75, 0, 0, 0, 0, 0, 0, 0, 0, t75 * t104, -t66 * t76 + t53, (-t81 * t92 + t78) * t65 + t106, t48 * t78 + qJ(3) * t106 + ((-pkin(2) * qJD(2) - t47) * t69 - t48 * t81) * t92, t6, t1, -t36, t7, -t37, 0, t97 * qJD(4) + t60 * t27 - t30 * t88 + t96, -t98 * qJD(4) - t60 * t26 - t32 * t88 + t95, t20 * t26 - t21 * t27 - t98 * t30 - t97 * t32 + t84, t5 * t20 + t4 * t21 + t98 * t13 + t97 * t12 + (qJD(2) * t60 - t29) * t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t93 * t65 ^ 2, -t93 * t65 * t48 + t57, 0, 0, 0, 0, 0, 0, 0.2e1 * t32 * qJD(4), t49 + (-t30 - t90) * qJD(4), -t28 - t103, t12 * t32 + t13 * t30 + t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t101, -t28 + t103, t49 + (t30 - t90) * qJD(4), -t101, 0, 0, -t29 * t32 - t72, t29 * t30 + t105, 0, 0;];
tauc_reg = t2;
