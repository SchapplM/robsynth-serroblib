% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RPRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2,theta4]';
% 
% Output:
% tauc_reg [5x17]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:09
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRPP1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPP1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:09:13
% EndTime: 2019-12-31 18:09:15
% DurationCPUTime: 0.47s
% Computational Cost: add. (916->130), mult. (2187->182), div. (0->0), fcn. (1368->6), ass. (0->80)
t104 = 2 * qJD(3);
t58 = sin(pkin(7)) * pkin(1) + pkin(6);
t95 = qJ(4) + t58;
t64 = sin(pkin(8));
t66 = cos(pkin(8));
t68 = sin(qJ(3));
t69 = cos(qJ(3));
t49 = t64 * t69 + t66 * t68;
t44 = t49 * qJD(1);
t39 = t44 ^ 2;
t94 = qJD(1) * t69;
t86 = t66 * t94;
t93 = t68 * qJD(1);
t41 = t64 * t93 - t86;
t103 = -t41 ^ 2 - t39;
t82 = t95 * qJD(1);
t32 = t69 * qJD(2) - t82 * t68;
t33 = t68 * qJD(2) + t82 * t69;
t88 = qJD(1) * qJD(4);
t24 = qJD(3) * t32 + t69 * t88;
t72 = -qJD(3) * t33 - t68 * t88;
t2 = t64 * t24 - t66 * t72;
t47 = t95 * t69;
t84 = t95 * t68;
t22 = t64 * t47 + t66 * t84;
t102 = t2 * t22;
t99 = t66 * t69;
t48 = t64 * t68 - t99;
t101 = t2 * t48;
t100 = t64 * t33;
t28 = t66 * t33;
t70 = qJD(3) ^ 2;
t98 = t70 * t68;
t97 = t70 * t69;
t3 = t66 * t24 + t64 * t72;
t30 = (qJD(3) * pkin(3)) + t32;
t9 = t64 * t30 + t28;
t96 = t68 ^ 2 - t69 ^ 2;
t60 = -cos(pkin(7)) * pkin(1) - pkin(2);
t51 = qJD(1) * t60;
t91 = t68 * qJD(3);
t13 = t66 * t32 - t100;
t90 = qJD(5) - t13;
t89 = qJD(1) * qJD(3);
t87 = pkin(3) * t91;
t85 = t68 * t89;
t83 = qJD(3) * t95;
t8 = t66 * t30 - t100;
t80 = -t69 * pkin(3) + t60;
t79 = t51 * t104;
t43 = t49 * qJD(3);
t35 = qJD(1) * t43;
t52 = t64 * t85;
t36 = qJD(3) * t86 - t52;
t55 = pkin(3) * t85;
t78 = t35 * pkin(4) - t36 * qJ(5) + t55;
t75 = t80 * qJD(1);
t40 = qJD(4) + t75;
t16 = t41 * pkin(4) - t44 * qJ(5) + t40;
t77 = t16 * t44 + t2;
t46 = qJD(3) * t99 - t64 * t91;
t76 = -t49 * t35 + t48 * t36 - t46 * t41 + t43 * t44;
t74 = -t68 * qJD(4) - t69 * t83;
t34 = t69 * qJD(4) - t68 * t83;
t14 = t64 * t34 - t66 * t74;
t15 = t66 * t34 + t64 * t74;
t23 = t66 * t47 - t64 * t84;
t73 = t14 * t44 - t15 * t41 + t2 * t49 + t22 * t36 - t23 * t35;
t71 = qJD(1) ^ 2;
t59 = -t66 * pkin(3) - pkin(4);
t56 = t64 * pkin(3) + qJ(5);
t18 = t48 * pkin(4) - t49 * qJ(5) + t80;
t17 = pkin(3) * t93 + t44 * pkin(4) + t41 * qJ(5);
t12 = t64 * t32 + t28;
t10 = t43 * pkin(4) - t46 * qJ(5) - t49 * qJD(5) + t87;
t7 = qJD(3) * qJ(5) + t9;
t6 = -qJD(3) * pkin(4) + qJD(5) - t8;
t5 = -t44 * qJD(5) + t78;
t1 = qJD(3) * qJD(5) + t3;
t4 = [0, 0, 0, 0, 0.2e1 * t69 * t85, -0.2e1 * t96 * t89, t97, -t98, 0, -t58 * t97 + t68 * t79, t58 * t98 + t69 * t79, -t3 * t48 - t9 * t43 - t8 * t46 + t73, -t8 * t14 + t9 * t15 + t102 + t3 * t23 + (t40 + t75) * t87, -t14 * qJD(3) + t10 * t41 + t16 * t43 + t18 * t35 + t5 * t48, -t1 * t48 - t7 * t43 + t6 * t46 + t73, t15 * qJD(3) - t10 * t44 - t16 * t46 - t18 * t36 - t5 * t49, t1 * t23 + t16 * t10 + t6 * t14 + t7 * t15 + t5 * t18 + t102; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t98, -t97, t76, t3 * t49 - t8 * t43 + t9 * t46 + t101, -t43 * qJD(3), t76, t46 * qJD(3), t1 * t49 + t6 * t43 + t7 * t46 + t101; 0, 0, 0, 0, -t68 * t71 * t69, t96 * t71, 0, 0, 0, -t51 * t93, -t51 * t94, (-t12 + t9) * t44 + (t13 - t8) * t41 + (-t35 * t64 - t36 * t66) * pkin(3), t8 * t12 - t9 * t13 + (-t2 * t66 + t3 * t64 - t40 * t93) * pkin(3), t12 * qJD(3) - t17 * t41 - t77, -t56 * t35 + t59 * t36 + (-t12 + t7) * t44 + (t6 - t90) * t41, -t16 * t41 + t17 * t44 + (0.2e1 * qJD(5) - t13) * qJD(3) + t3, t1 * t56 - t6 * t12 - t16 * t17 + t2 * t59 + t90 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t103, t9 * t41 + t8 * t44 + t55, t44 * t104, t103, t52 + (t41 - t86) * qJD(3), t7 * t41 + (-qJD(5) - t6) * t44 + t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44 * t41, -t52 + (t41 + t86) * qJD(3), -t39 - t70, -t7 * qJD(3) + t77;];
tauc_reg = t4;
