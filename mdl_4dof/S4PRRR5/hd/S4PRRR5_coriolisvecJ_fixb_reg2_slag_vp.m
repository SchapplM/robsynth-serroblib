% Calculate inertial parameters regressor of coriolis joint torque vector for
% S4PRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% 
% Output:
% tauc_reg [4x(4*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4PRRR5_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR5_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR5_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR5_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:33:46
% EndTime: 2019-12-31 16:33:47
% DurationCPUTime: 0.34s
% Computational Cost: add. (417->72), mult. (1016->113), div. (0->0), fcn. (676->6), ass. (0->67)
t36 = cos(qJ(2));
t61 = t36 * qJD(1);
t23 = qJD(2) * pkin(2) + t61;
t33 = sin(qJ(2));
t32 = sin(qJ(3));
t35 = cos(qJ(3));
t18 = t32 * t36 + t35 * t33;
t43 = t18 * qJD(2);
t63 = qJD(3) * t35;
t39 = (t33 * t63 + t43) * qJD(1);
t64 = qJD(3) * t32;
t5 = t23 * t64 + t39;
t51 = qJD(2) * t61;
t28 = qJD(2) + qJD(3);
t65 = qJD(1) * t33;
t56 = t32 * t65;
t68 = t28 * t56;
t4 = (qJD(3) * t23 + t51) * t35 - t68;
t31 = sin(qJ(4));
t29 = t31 ^ 2;
t34 = cos(qJ(4));
t30 = t34 ^ 2;
t66 = t29 + t30;
t81 = t66 * t4;
t17 = t32 * t33 - t35 * t36;
t80 = t17 * t28;
t16 = t17 * qJD(1);
t79 = pkin(2) * t63 + t16;
t24 = t32 * pkin(2) + pkin(6);
t37 = qJD(4) ^ 2;
t15 = t18 * qJD(1);
t49 = pkin(2) * t64 - t15;
t77 = t24 * t37 + t49 * t28;
t12 = t35 * t23 - t56;
t74 = t28 * pkin(3);
t10 = -t12 - t74;
t62 = qJD(4) * t34;
t75 = t10 * t62 + t5 * t31;
t73 = t5 * t17;
t7 = t18 * qJD(3) + t43;
t72 = t7 * t28;
t13 = t32 * t23 + t35 * t65;
t71 = t13 * t28;
t69 = t37 * t31;
t67 = t29 - t30;
t27 = t28 ^ 2;
t60 = t31 * t27 * t34;
t54 = t66 * t80;
t53 = -t10 * t28 - t4;
t52 = t66 * t12;
t50 = t31 * t28 * t62;
t48 = pkin(6) * t37 - t71;
t47 = t18 * t37 + t72;
t46 = qJD(4) * (t12 - t74);
t45 = 0.2e1 * qJD(4) * t80;
t44 = (-pkin(2) * t28 - t23) * qJD(3);
t25 = -t35 * pkin(2) - pkin(3);
t41 = qJD(4) * (t25 * t28 - t79);
t40 = t79 * t66;
t38 = qJD(2) ^ 2;
t26 = t37 * t34;
t20 = -0.2e1 * t50;
t19 = 0.2e1 * t50;
t14 = -0.2e1 * t67 * t28 * qJD(4);
t11 = t28 * pkin(6) + t13;
t8 = t10 * qJD(4) * t31;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38 * t33, -t38 * t36, 0, 0, 0, 0, 0, 0, 0, 0, -t72, t80 * t28, 0, -t12 * t7 - t13 * t80 + t4 * t18 + t73, 0, 0, 0, 0, 0, 0, t31 * t45 - t47 * t34, t47 * t31 + t34 * t45, -t28 * t54, t10 * t7 - t11 * t54 + t18 * t81 + t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15 * t28 + t32 * t44 - t39, -t16 * t28 + (t44 - t51) * t35 + t68, 0, t12 * t15 + t13 * t16 + (t32 * t4 - t35 * t5 + (-t12 * t32 + t13 * t35) * qJD(3)) * pkin(2), t19, t14, t26, t20, -t69, 0, t8 + t31 * t41 + (-t5 - t77) * t34, t77 * t31 + t34 * t41 + t75, t40 * t28 + t81, t49 * t10 + t40 * t11 + t24 * t81 + t5 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71 - t5, t12 * t28 - t4, 0, 0, t19, t14, t26, t20, -t69, 0, t8 + t31 * t46 + (-t48 - t5) * t34, t48 * t31 + t34 * t46 + t75, -t28 * t52 + t81, -t5 * pkin(3) + pkin(6) * t81 - t10 * t13 - t11 * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t60, t67 * t27, 0, t60, 0, 0, t53 * t31, t53 * t34, 0, 0;];
tauc_reg = t1;
