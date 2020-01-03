% Calculate minimal parameter regressor of coriolis joint torque vector for
% S4RPRR7
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
% tauc_reg [4x21]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RPRR7_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR7_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR7_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR7_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:54:04
% EndTime: 2019-12-31 16:54:07
% DurationCPUTime: 0.57s
% Computational Cost: add. (694->126), mult. (1951->187), div. (0->0), fcn. (1456->6), ass. (0->73)
t50 = cos(pkin(7));
t54 = cos(qJ(3));
t78 = t54 * t50;
t42 = qJD(1) * t78;
t49 = sin(pkin(7));
t52 = sin(qJ(3));
t79 = t52 * t49;
t71 = qJD(1) * t79;
t30 = t42 - t71;
t29 = qJD(4) - t30;
t91 = qJD(4) - t29;
t35 = t54 * t49 + t52 * t50;
t90 = t35 * qJD(1);
t77 = pkin(5) + qJ(2);
t39 = t77 * t49;
t36 = qJD(1) * t39;
t40 = t77 * t50;
t37 = qJD(1) * t40;
t19 = -t52 * t36 + t54 * t37;
t59 = t35 * qJD(2);
t4 = qJD(1) * t59 + t19 * qJD(3);
t89 = (t90 * pkin(3) + t29 * pkin(6)) * t29 + t4;
t51 = sin(qJ(4));
t53 = cos(qJ(4));
t24 = t51 * qJD(3) + t53 * t90;
t41 = qJD(3) * t42;
t27 = -qJD(3) * t71 + t41;
t6 = qJD(4) * t24 + t51 * t27;
t18 = -t54 * t36 - t52 * t37;
t13 = -qJD(3) * pkin(3) - t18;
t34 = -t78 + t79;
t44 = -t50 * pkin(2) - pkin(1);
t17 = t34 * pkin(3) - t35 * pkin(6) + t44;
t21 = -t52 * t39 + t54 * t40;
t33 = t35 * qJD(3);
t28 = qJD(1) * t33;
t58 = t34 * qJD(2);
t3 = -qJD(1) * t58 + t18 * qJD(3);
t32 = t34 * qJD(3);
t38 = t44 * qJD(1) + qJD(2);
t7 = -t30 * pkin(3) - pkin(6) * t90 + t38;
t64 = -t54 * t39 - t52 * t40;
t8 = t64 * qJD(3) - t58;
t88 = -t13 * t32 - t21 * t28 - (qJD(4) * t17 + t8) * t29 - (qJD(4) * t7 + t3) * t34 + t4 * t35;
t73 = t53 * qJD(3);
t74 = qJD(4) * t51;
t5 = qJD(4) * t73 + t53 * t27 - t74 * t90;
t87 = t5 * t51;
t86 = t17 * t28;
t22 = t51 * t90 - t73;
t85 = t22 * t29;
t84 = t24 * t29;
t83 = t24 * t90;
t82 = t90 * t22;
t80 = t51 * t28;
t26 = t53 * t28;
t76 = t49 ^ 2 + t50 ^ 2;
t75 = qJD(4) * t35;
t72 = qJD(1) * qJD(2);
t70 = t76 * qJD(1) ^ 2;
t66 = t29 * t53;
t14 = qJD(3) * pkin(6) + t19;
t2 = t53 * t14 + t51 * t7;
t65 = t51 * t14 - t53 * t7;
t63 = 0.2e1 * t76 * t72;
t62 = t26 + (t30 * t51 - t74) * t29;
t61 = -t53 * t32 - t35 * t74;
t56 = -pkin(6) * t28 + (t13 + t18) * t29;
t16 = t33 * pkin(3) + t32 * pkin(6);
t11 = t28 * pkin(3) - t27 * pkin(6);
t10 = t53 * t11;
t9 = t21 * qJD(3) + t59;
t1 = [0, 0, 0, 0, 0, t63, qJ(2) * t63, t27 * t35 - t32 * t90, -t27 * t34 - t35 * t28 - t32 * t30 - t33 * t90, -t32 * qJD(3), -t33 * qJD(3), 0, -t9 * qJD(3) + t44 * t28 + t38 * t33, -t8 * qJD(3) + t44 * t27 - t38 * t32, t5 * t53 * t35 + t61 * t24, -(-t22 * t53 - t24 * t51) * t32 + (-t87 - t53 * t6 + (t22 * t51 - t24 * t53) * qJD(4)) * t35, t24 * t33 + t35 * t26 + t61 * t29 + t5 * t34, -t35 * t80 - t22 * t33 - t6 * t34 + (t51 * t32 - t53 * t75) * t29, t28 * t34 + t29 * t33, -t65 * t33 + t10 * t34 - t64 * t6 + t9 * t22 + (t16 * t29 + t86 + (t13 * t35 - t14 * t34 - t21 * t29) * qJD(4)) * t53 + t88 * t51, -t2 * t33 - t64 * t5 + t9 * t24 + (-(-qJD(4) * t21 + t16) * t29 - t86 - (-qJD(4) * t14 + t11) * t34 - t13 * t75) * t51 + t88 * t53; 0, 0, 0, 0, 0, -t70, -qJ(2) * t70, 0, 0, 0, 0, 0, 0.2e1 * t90 * qJD(3), t41 + (t30 - t71) * qJD(3), 0, 0, 0, 0, 0, t62 - t82, -t29 ^ 2 * t53 - t80 - t83; 0, 0, 0, 0, 0, 0, 0, -t90 * t30, -t30 ^ 2 + t90 ^ 2, t41 + (-t30 - t71) * qJD(3), 0, 0, (-qJD(2) - t38) * t90, -t38 * t30 + t34 * t72, t24 * t66 + t87, (t5 - t85) * t53 + (-t6 - t84) * t51, t29 * t66 + t80 - t83, t62 + t82, -t29 * t90, -pkin(3) * t6 - t19 * t22 + t56 * t51 - t89 * t53 + t65 * t90, -pkin(3) * t5 - t19 * t24 + t2 * t90 + t89 * t51 + t56 * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24 * t22, -t22 ^ 2 + t24 ^ 2, t5 + t85, -t6 + t84, t28, -t13 * t24 - t91 * t2 - t51 * t3 + t10, -t51 * t11 + t13 * t22 - t53 * t3 + t91 * t65;];
tauc_reg = t1;
