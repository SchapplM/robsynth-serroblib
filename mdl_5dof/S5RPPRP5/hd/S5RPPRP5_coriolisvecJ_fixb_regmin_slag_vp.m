% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RPPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
% 
% Output:
% tauc_reg [5x22]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPPRP5_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP5_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP5_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP5_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:53:43
% EndTime: 2019-12-31 17:53:44
% DurationCPUTime: 0.44s
% Computational Cost: add. (594->132), mult. (1545->180), div. (0->0), fcn. (1029->4), ass. (0->78)
t59 = sin(pkin(7));
t60 = cos(pkin(7));
t62 = sin(qJ(4));
t63 = cos(qJ(4));
t35 = t59 * t62 + t60 * t63;
t90 = qJD(1) * t35;
t99 = t90 ^ 2;
t89 = qJD(1) * t59;
t78 = t63 * t89;
t88 = qJD(1) * t60;
t79 = t62 * t88;
t31 = t78 - t79;
t98 = t31 ^ 2;
t25 = -qJD(1) * pkin(1) - pkin(2) * t88 - qJ(3) * t89 + qJD(2);
t19 = pkin(3) * t88 - t25;
t3 = pkin(4) * t90 - t31 * qJ(5) + t19;
t97 = t3 * t31;
t96 = t31 * t90;
t94 = -pkin(6) + qJ(2);
t40 = t94 * t60;
t37 = qJD(1) * t40;
t95 = t62 * t37;
t57 = t59 ^ 2;
t58 = t60 ^ 2;
t93 = t57 + t58;
t39 = t94 * t59;
t70 = t63 * t39 - t62 * t40;
t6 = t35 * qJD(2) + t70 * qJD(4);
t92 = t6 * qJD(4);
t17 = t62 * t39 + t63 * t40;
t36 = t59 * t63 - t60 * t62;
t7 = -t36 * qJD(2) + t17 * qJD(4);
t91 = t7 * qJD(4);
t87 = qJD(4) * t62;
t86 = qJD(4) * t63;
t85 = t59 * qJD(3);
t46 = qJ(2) * t89 + qJD(3);
t34 = -pkin(6) * t89 + t46;
t14 = t63 * t34 - t95;
t84 = qJD(5) - t14;
t83 = qJD(1) * qJD(2);
t82 = qJD(1) * qJD(3);
t81 = -t60 * pkin(2) - t59 * qJ(3) - pkin(1);
t76 = t59 * t83;
t77 = t60 * t83;
t80 = t34 * t86 + t62 * t76 + t63 * t77;
t65 = qJD(1) ^ 2;
t38 = t93 * t65;
t75 = t59 * t82;
t74 = t14 + t95;
t73 = qJ(2) * t83;
t33 = t60 * pkin(3) - t81;
t5 = t34 * t87 + t37 * t86 + t62 * t77 - t63 * t76;
t72 = 0.2e1 * t90;
t71 = t58 * t73;
t15 = t62 * t34 + t63 * t37;
t64 = qJD(4) ^ 2;
t69 = t31 * t89 + t64 * t63;
t68 = t15 * qJD(4) - t5;
t27 = t35 * qJD(4);
t45 = qJD(4) * t79;
t67 = t45 + (-t31 - t78) * qJD(4);
t22 = qJD(1) * t27;
t23 = qJD(4) * t78 - t45;
t66 = -t23 * pkin(4) - t22 * qJ(5) - t75;
t47 = t57 * t73;
t28 = t59 * t86 - t60 * t87;
t26 = 0.2e1 * t93 * t83;
t18 = -t64 * t62 - t89 * t90;
t13 = t31 * pkin(4) + qJ(5) * t90;
t11 = t72 * qJD(4);
t10 = qJD(4) * qJ(5) + t15;
t9 = -qJD(4) * pkin(4) + t84;
t8 = t35 * pkin(4) - t36 * qJ(5) + t33;
t4 = t28 * pkin(4) + t27 * qJ(5) - t36 * qJD(5) + t85;
t2 = (qJD(5) - t95) * qJD(4) + t80;
t1 = -t31 * qJD(5) - t66;
t12 = [0, 0, 0, 0, 0, t26, 0.2e1 * t47 + 0.2e1 * t71, 0.2e1 * t60 * t75, t26, 0.2e1 * t57 * t82, 0.2e1 * t71 + t47 + (t46 * qJD(2) + (-qJD(1) * t81 - t25) * qJD(3)) * t59, -t22 * t36 - t31 * t27, t22 * t35 - t36 * t23 + t27 * t90 - t31 * t28, -t27 * qJD(4), -t28 * qJD(4), 0, t19 * t28 + t33 * t23 + t72 * t85 - t91, -t92 - t19 * t27 - t33 * t22 + (qJD(1) * t36 + t31) * t85, t1 * t35 + t8 * t23 + t3 * t28 + t4 * t90 - t91, -t10 * t28 - t17 * t23 - t2 * t35 + t22 * t70 - t9 * t27 + t7 * t31 + t5 * t36 - t6 * t90, -t1 * t36 + t8 * t22 + t3 * t27 - t4 * t31 + t92, t1 * t8 + t10 * t6 + t2 * t17 + t3 * t4 - t5 * t70 + t9 * t7; 0, 0, 0, 0, 0, -t38, -qJ(2) * t38, 0, -t38, 0, -t58 * t65 * qJ(2) + (-qJD(3) - t46) * t89, 0, 0, 0, 0, 0, t67, t11, t67, t98 + t99, -t11, -t10 * t90 + (qJD(5) + t9) * t31 + t66; 0, 0, 0, 0, 0, 0, 0, -t59 * t65 * t60, 0, -t57 * t65, (qJD(2) + t25) * t89, 0, 0, 0, 0, 0, t18, -t69, t18, t63 * t22 - t62 * t23 + (t31 * t62 - t63 * t90) * qJD(4), t69, -t3 * t89 + t2 * t62 - t5 * t63 + (t10 * t63 + t62 * t9) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t96, t98 - t99, 0, t45 + (t31 - t78) * qJD(4), 0, -t19 * t31 + t68, t74 * qJD(4) + t19 * t90 - t80, -t13 * t90 + t68 - t97, pkin(4) * t22 - t23 * qJ(5) + (t10 - t15) * t31 + (t9 - t84) * t90, t13 * t31 - t3 * t90 + (0.2e1 * qJD(5) - t74) * qJD(4) + t80, -t5 * pkin(4) + t2 * qJ(5) + t84 * t10 - t3 * t13 - t9 * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t96, 0, -t64 - t98, -t10 * qJD(4) + t5 + t97;];
tauc_reg = t12;
