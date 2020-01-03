% Calculate minimal parameter regressor of coriolis joint torque vector for
% S4RPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% 
% Output:
% tauc_reg [4x18]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RPRR4_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR4_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR4_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR4_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:50:35
% EndTime: 2019-12-31 16:50:37
% DurationCPUTime: 0.52s
% Computational Cost: add. (367->111), mult. (936->180), div. (0->0), fcn. (568->6), ass. (0->73)
t36 = cos(qJ(3));
t59 = t36 * qJD(1);
t24 = -qJD(4) + t59;
t85 = qJD(4) + t24;
t34 = sin(qJ(3));
t29 = t34 ^ 2;
t42 = qJD(1) * t29 - t24 * t36;
t33 = sin(qJ(4));
t63 = qJD(4) * t33;
t56 = t34 * t63;
t35 = cos(qJ(4));
t60 = t35 * qJD(3);
t84 = -t24 * t56 - t42 * t60;
t62 = qJD(4) * t35;
t55 = t34 * t62;
t57 = qJD(3) * qJD(4);
t61 = t33 * qJD(3);
t4 = (t36 * t61 + t55) * qJD(1) + t33 * t57;
t53 = t36 * t60;
t65 = qJD(1) * t34;
t54 = t33 * t65;
t3 = qJD(1) * t53 - qJD(4) * t54 + t35 * t57;
t83 = t3 * t33;
t25 = sin(pkin(7)) * pkin(1) + pkin(5);
t21 = t25 * qJD(1);
t10 = t36 * qJD(2) - t34 * t21;
t5 = -qJD(3) * pkin(3) - t10;
t82 = t33 * t5;
t26 = -cos(pkin(7)) * pkin(1) - pkin(2);
t15 = -t36 * pkin(3) - t34 * pkin(6) + t26;
t9 = t15 * qJD(1);
t81 = t33 * t9;
t80 = t35 * t5;
t79 = t35 * t9;
t78 = t36 * t4;
t11 = t34 * qJD(2) + t36 * t21;
t8 = t11 * qJD(3);
t77 = t8 * t33;
t76 = t8 * t35;
t16 = t54 - t60;
t75 = t16 * t24;
t74 = t16 * t34;
t18 = t35 * t65 + t61;
t73 = t18 * t24;
t72 = t24 * t33;
t71 = t33 * t36;
t70 = t35 * t24;
t69 = t35 * t36;
t37 = qJD(3) ^ 2;
t68 = t37 * t34;
t67 = t37 * t36;
t66 = -t36 ^ 2 + t29;
t22 = qJD(1) * t26;
t64 = qJD(3) * t34;
t58 = qJD(1) * qJD(3);
t45 = pkin(3) * t34 - pkin(6) * t36;
t20 = t45 * qJD(3);
t14 = qJD(1) * t20;
t7 = t10 * qJD(3);
t52 = -t35 * t14 + t33 * t7;
t51 = t18 * t64 - t3 * t36;
t6 = qJD(3) * pkin(6) + t11;
t50 = t24 * t25 + t6;
t48 = t34 * t58;
t46 = t24 * t55;
t2 = t35 * t6 + t81;
t44 = t33 * t6 - t79;
t43 = t33 * t14 + t35 * t7;
t41 = 0.2e1 * qJD(3) * t22;
t40 = t42 * t33;
t38 = qJD(1) ^ 2;
t19 = t45 * qJD(1);
t1 = [0, 0, 0, 0, 0.2e1 * t36 * t48, -0.2e1 * t66 * t58, t67, -t68, 0, -t25 * t67 + t34 * t41, t25 * t68 + t36 * t41, t3 * t35 * t34 + (t53 - t56) * t18, (-t16 * t35 - t18 * t33) * t36 * qJD(3) + (-t83 - t35 * t4 + (t16 * t33 - t18 * t35) * qJD(4)) * t34, t51 - t84, t46 + t78 + (-t40 - t74) * qJD(3), (-t24 - t59) * t64, -(-t15 * t63 + t35 * t20) * t24 + ((t16 * t25 + t82) * qJD(3) + (t50 * t35 + t81) * qJD(4) + t52) * t36 + (t5 * t62 + t25 * t4 + t77 + (-t25 * t72 + (t35 * t15 - t25 * t71) * qJD(1) - t44) * qJD(3)) * t34, (t15 * t62 + t33 * t20) * t24 + ((t18 * t25 + t80) * qJD(3) + (-t50 * t33 + t79) * qJD(4) + t43) * t36 + (-t5 * t63 + t25 * t3 + t76 + (-t25 * t70 - (t33 * t15 + t25 * t69) * qJD(1) - t2) * qJD(3)) * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t68, -t67, 0, 0, 0, 0, 0, t46 - t78 + (-t40 + t74) * qJD(3), t51 + t84; 0, 0, 0, 0, -t34 * t38 * t36, t66 * t38, 0, 0, 0, -t22 * t65, -t22 * t59, -t18 * t70 + t83, (t3 + t75) * t35 + (-t4 + t73) * t33, -t24 * t62 + (t24 * t69 + (-t18 + t61) * t34) * qJD(1), t24 * t63 + (-t24 * t71 + (t16 + t60) * t34) * qJD(1), t24 * t65, -pkin(3) * t4 - t76 + (-t33 * t10 + t35 * t19) * t24 - t11 * t16 + (pkin(6) * t70 + t82) * qJD(4) + (t44 * t34 + (-pkin(6) * t64 - t36 * t5) * t33) * qJD(1), -pkin(3) * t3 + t77 - (t35 * t10 + t33 * t19) * t24 - t11 * t18 + (-pkin(6) * t72 + t80) * qJD(4) + (-t5 * t69 + (-pkin(6) * t60 + t2) * t34) * qJD(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18 * t16, -t16 ^ 2 + t18 ^ 2, t3 - t75, -t4 - t73, t48, -t5 * t18 - t85 * t2 - t52, t5 * t16 + t85 * t44 - t43;];
tauc_reg = t1;
