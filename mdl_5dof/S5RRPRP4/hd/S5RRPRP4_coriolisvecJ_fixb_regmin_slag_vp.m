% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% 
% Output:
% tauc_reg [5x20]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:53
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRPRP4_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP4_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP4_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP4_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:52:56
% EndTime: 2019-12-31 19:52:57
% DurationCPUTime: 0.45s
% Computational Cost: add. (644->127), mult. (973->149), div. (0->0), fcn. (401->4), ass. (0->80)
t72 = qJD(4) * qJ(5);
t42 = qJD(1) + qJD(2);
t49 = -pkin(2) - pkin(7);
t48 = cos(qJ(2));
t77 = pkin(1) * qJD(1);
t68 = t48 * t77;
t58 = qJD(3) - t68;
t15 = t49 * t42 + t58;
t45 = sin(qJ(4));
t86 = t45 * t15;
t10 = t72 + t86;
t46 = sin(qJ(2));
t76 = pkin(1) * qJD(2);
t66 = qJD(1) * t76;
t61 = t46 * t66;
t30 = t45 * t61;
t47 = cos(qJ(4));
t85 = t47 * t15;
t3 = t30 + (qJD(5) + t85) * qJD(4);
t31 = t47 * t61;
t74 = qJD(4) * t45;
t5 = t15 * t74 - t31;
t63 = qJD(4) * pkin(4) - qJD(5);
t9 = -t63 - t85;
t51 = t3 * t45 - t5 * t47 + (t10 * t47 + t45 * t9) * qJD(4);
t29 = t45 * pkin(4) - t47 * qJ(5) + qJ(3);
t95 = t29 * t42;
t43 = t45 ^ 2;
t44 = t47 ^ 2;
t94 = t42 * (t43 + t44);
t41 = t42 ^ 2;
t93 = 0.2e1 * qJD(4);
t35 = t48 * t66;
t55 = pkin(4) * t47 + qJ(5) * t45;
t62 = -t47 * qJD(5) + qJD(3);
t2 = t35 + (t55 * qJD(4) + t62) * t42;
t73 = qJD(4) * t47;
t69 = t46 * t77;
t8 = t69 + t95;
t92 = t2 * t45 + t8 * t73;
t90 = t46 * pkin(1);
t67 = -t48 * pkin(1) - pkin(2);
t36 = -pkin(7) + t67;
t50 = qJD(4) ^ 2;
t89 = t36 * t50;
t88 = t42 * t45;
t87 = t42 * t47;
t84 = t49 * t50;
t83 = t50 * t45;
t82 = t50 * t47;
t23 = t42 * qJD(3) + t35;
t75 = t42 * qJ(3);
t27 = t69 + t75;
t81 = t23 * t45 + t27 * t73;
t80 = t41 + t50;
t79 = t43 - t44;
t71 = t46 * t76;
t70 = t48 * t76;
t64 = t73 * t90;
t60 = qJD(2) * t64 - t36 * t83;
t14 = pkin(4) * t73 + t45 * t72 + t62;
t59 = -t14 + t68;
t56 = t10 * t45 - t47 * t9;
t54 = -qJD(1) * t64 - t49 * t83;
t37 = qJ(3) + t90;
t34 = qJD(3) + t70;
t33 = t47 * t41 * t45;
t28 = -0.2e1 * t73 * t88;
t26 = t80 * t47;
t25 = t80 * t45;
t24 = -t42 * pkin(2) + t58;
t22 = t29 + t90;
t20 = t55 * t42;
t19 = t23 * t47;
t17 = (qJD(1) + t42) * t71;
t16 = (qJD(2) - t42) * t69;
t13 = t79 * t42 * t93;
t11 = t14 + t70;
t6 = t8 * t74;
t1 = [0, 0, 0, 0, -t17, -t42 * t70 - t35, t17, t35 + (qJD(3) + t34) * t42, t23 * t37 + t27 * t34 + (t67 * qJD(1) + t24) * t71, t28, t13, -t83, -t82, 0, (t34 * t45 + t37 * t73) * t42 + t60 + t81, t19 + (t34 * t42 - t89) * t47 + (-t37 * t42 - t27 - t71) * t74, (t11 * t45 + t22 * t73) * t42 + t60 + t92, -t71 * t94 - t51, t6 + (t22 * t42 + t71) * t74 + (-t11 * t42 - t2 + t89) * t47, t8 * t11 + t2 * t22 + t51 * t36 + t56 * t71; 0, 0, 0, 0, -t16, t42 * t68 - t35, t16, t35 + (0.2e1 * qJD(3) - t68) * t42, t23 * qJ(3) + t27 * qJD(3) + (-t27 * t48 + (-pkin(2) * qJD(2) - t24) * t46) * t77, t28, t13, -t83, -t82, 0, (qJ(3) * t73 + t58 * t45) * t42 + t54 + t81, t19 + (t58 * t42 - t84) * t47 + (-t27 + t69 - t75) * t74, (t29 * t73 - t59 * t45) * t42 + t54 + t92, t69 * t94 - t51, t6 + (-t69 + t95) * t74 + (t59 * t42 - t2 + t84) * t47, t8 * t14 + t2 * t29 + (-t56 * t46 - t48 * t8) * t77 + t51 * t49; 0, 0, 0, 0, 0, 0, 0, -t41, -t27 * t42 + t61, 0, 0, 0, 0, 0, -t25, -t26, -t25, 0, t26, -t8 * t42 + t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, -t79 * t41, 0, 0, 0, -t27 * t87 + t31, t27 * t88 - t30, t31 + (-t20 * t45 - t47 * t8) * t42, ((t10 - t72) * t47 + (t63 + t9) * t45) * t42, qJD(5) * t93 + t30 + (t20 * t47 - t45 * t8) * t42, -t9 * t86 - t5 * pkin(4) + t3 * qJ(5) - t8 * t20 + (qJD(5) - t85) * t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, 0, -t44 * t41 - t50, -t10 * qJD(4) + t8 * t87 + t5;];
tauc_reg = t1;
