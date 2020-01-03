% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RPRRP5
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
% Datum: 2019-12-31 18:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRRP5_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP5_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP5_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP5_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:40:57
% EndTime: 2019-12-31 18:40:58
% DurationCPUTime: 0.38s
% Computational Cost: add. (782->112), mult. (1609->135), div. (0->0), fcn. (842->6), ass. (0->84)
t56 = sin(qJ(3));
t58 = cos(qJ(3));
t103 = pkin(1) * sin(pkin(8));
t75 = qJD(1) * t103;
t69 = qJD(3) * t75;
t44 = cos(pkin(8)) * pkin(1) + pkin(2);
t41 = t44 * qJD(1);
t83 = qJD(3) * t41;
t25 = -t56 * t69 + t58 * t83;
t107 = qJD(2) * qJD(4) + t25;
t57 = cos(qJ(4));
t28 = t56 * t41 + t58 * t75;
t50 = qJD(1) + qJD(3);
t22 = t50 * pkin(7) + t28;
t55 = sin(qJ(4));
t92 = t55 * t22;
t12 = t57 * qJD(2) - t92;
t106 = qJD(5) - t12;
t87 = t107 * t57;
t2 = (qJD(5) - t92) * qJD(4) + t87;
t81 = qJD(4) * t57;
t4 = t107 * t55 + t22 * t81;
t105 = t2 * t57 + t4 * t55;
t9 = -qJD(4) * pkin(4) + t106;
t13 = t55 * qJD(2) + t57 * t22;
t78 = qJD(4) * qJ(5);
t10 = t13 + t78;
t104 = t56 * t103 - t58 * t44;
t59 = qJD(4) ^ 2;
t102 = pkin(7) * t59;
t101 = t50 * pkin(3);
t27 = t58 * t41 - t56 * t75;
t100 = t27 * t50;
t99 = t28 * t50;
t29 = t104 * qJD(3);
t98 = t29 * t50;
t63 = t58 * t103 + t56 * t44;
t30 = t63 * qJD(3);
t97 = t30 * t50;
t33 = pkin(7) + t63;
t96 = t33 * t59;
t38 = -t57 * pkin(4) - t55 * qJ(5) - pkin(3);
t95 = t38 * t50;
t94 = t50 * t55;
t93 = t50 * t57;
t90 = t59 * t55;
t21 = -t27 - t101;
t26 = t56 * t83 + t58 * t69;
t89 = t21 * t81 + t26 * t55;
t82 = qJD(4) * t55;
t88 = t27 * t82 + t28 * t93;
t79 = t55 * qJD(5);
t34 = pkin(4) * t82 - t57 * t78 - t79;
t86 = t28 - t34;
t51 = t55 ^ 2;
t52 = t57 ^ 2;
t85 = t51 - t52;
t84 = t51 + t52;
t80 = t10 * qJD(4);
t49 = t50 ^ 2;
t76 = t55 * t49 * t57;
t74 = t50 * t82;
t67 = pkin(4) * t55 - qJ(5) * t57;
t5 = (t67 * qJD(4) - t79) * t50 + t26;
t73 = -t5 - t102;
t72 = t12 + t92;
t23 = t38 + t104;
t71 = t23 * t50 + t29;
t68 = t10 * t57 + t55 * t9;
t66 = t96 + t97;
t65 = qJD(4) * ((-pkin(3) + t104) * t50 + t29);
t64 = t13 * qJD(4) - t4;
t11 = t30 + t34;
t62 = -t11 * t50 - t5 - t96;
t61 = -t55 * t80 + t9 * t81 + t105;
t60 = (-t10 * t55 + t57 * t9) * qJD(4) + t105;
t48 = t59 * t57;
t37 = 0.2e1 * t57 * t74;
t35 = t67 * t50;
t31 = -0.2e1 * t85 * t50 * qJD(4);
t14 = t21 * t82;
t8 = -t27 + t95;
t6 = t8 * t82;
t1 = [0, 0, 0, 0, 0, -t26 - t97, -t25 + t98, t37, t31, t48, -t90, 0, t14 + t55 * t65 + (-t26 - t66) * t57, t66 * t55 + t57 * t65 + t89, t62 * t57 + t71 * t82 + t6, -t84 * t98 + t61, t62 * t55 + (-t71 - t8) * t81, t8 * t11 + t5 * t23 - t68 * t29 + t60 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t90, -t48, -t90, 0, t48, t68 * qJD(4) + t2 * t55 - t4 * t57; 0, 0, 0, 0, 0, -t26 + t99, -t25 + t100, t37, t31, t48, -t90, 0, -pkin(3) * t74 + t14 + (-t26 - t102) * t57 + t88, (-t99 + t102) * t55 + (t27 - t101) * t81 + t89, t38 * t74 + t6 + (-t34 * t50 + t73) * t57 + t88, -t84 * t100 + t61, (-t27 - t8 - t95) * t81 + (t86 * t50 + t73) * t55, t60 * pkin(7) - t68 * t27 + t5 * t38 - t86 * t8; 0, 0, 0, 0, 0, 0, 0, -t76, t85 * t49, 0, 0, 0, -t21 * t94 + t64, t72 * qJD(4) - t21 * t93 - t87, (t35 * t57 - t55 * t8) * t50 + t64, 0, (t35 * t55 + t57 * t8) * t50 + (0.2e1 * qJD(5) - t72) * qJD(4) + t87, -t4 * pkin(4) + t2 * qJ(5) + t106 * t10 - t9 * t13 - t8 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t76, 0, -t51 * t49 - t59, t8 * t94 + t4 - t80;];
tauc_reg = t1;
