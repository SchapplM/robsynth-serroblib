% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RPPRP3
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
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPPRP3_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP3_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP3_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP3_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:51:02
% EndTime: 2019-12-31 17:51:03
% DurationCPUTime: 0.37s
% Computational Cost: add. (459->91), mult. (892->126), div. (0->0), fcn. (403->4), ass. (0->67)
t30 = -cos(pkin(7)) * pkin(1) - pkin(2) - pkin(6);
t18 = t30 * qJD(1) + qJD(3);
t41 = sin(qJ(4));
t42 = cos(qJ(4));
t10 = -t41 * qJD(2) + t42 * t18;
t75 = t41 * t18;
t11 = t42 * qJD(2) + t75;
t7 = t11 * qJD(4);
t58 = qJD(2) * qJD(4);
t66 = qJD(4) * t42;
t72 = -t18 * t66 + t41 * t58;
t79 = (t10 * t41 - t11 * t42) * qJD(4) + t72 * t41 + t7 * t42;
t57 = qJ(5) * t66;
t63 = t41 * qJD(5);
t1 = (-t57 - t63) * qJD(1) - t72;
t59 = qJD(1) * qJD(4);
t56 = t41 * t59;
t28 = qJ(5) * t56;
t62 = t42 * qJD(5);
t2 = -qJD(1) * t62 + t28 - t7;
t60 = qJ(5) * qJD(1);
t4 = -t42 * t60 + t10;
t70 = qJD(4) * pkin(4);
t3 = t4 + t70;
t5 = -t41 * t60 + t11;
t78 = -t1 * t41 - t2 * t42 + (t3 * t41 - t42 * t5) * qJD(4);
t77 = 0.2e1 * qJD(3);
t76 = t3 - t4;
t43 = qJD(4) ^ 2;
t35 = t43 * t41;
t74 = t43 * t42;
t44 = qJD(1) ^ 2;
t73 = t44 * t41;
t36 = qJD(3) * qJD(1);
t55 = t42 * t59;
t19 = pkin(4) * t55 + t36;
t37 = t41 ^ 2;
t38 = t42 ^ 2;
t71 = t37 - t38;
t69 = qJ(5) - t30;
t32 = sin(pkin(7)) * pkin(1) + qJ(3);
t20 = t41 * pkin(4) + t32;
t68 = qJD(1) * t20;
t24 = qJD(1) * t32;
t67 = qJD(4) * t41;
t14 = qJD(5) + t68;
t65 = t14 * qJD(1);
t64 = t24 * qJD(1);
t61 = qJD(5) + t14;
t16 = t69 * t42;
t54 = t14 + t68;
t27 = pkin(4) * t66 + qJD(3);
t53 = qJD(1) * t27 + t19;
t52 = t41 * t55;
t50 = t3 * t42 + t41 * t5;
t46 = t24 * t77;
t29 = t42 * t73;
t26 = -0.2e1 * t52;
t25 = 0.2e1 * t52;
t23 = (-t43 - t44) * t42;
t22 = -t35 - t73;
t21 = t71 * t44;
t17 = 0.2e1 * t71 * t59;
t15 = t69 * t41;
t9 = -qJD(4) * t16 - t63;
t8 = t69 * t67 - t62;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t36, t46, t26, t17, -t35, t25, -t74, 0, t24 * t66 - t30 * t35 + (t32 * t66 + t41 * t77) * qJD(1), -t24 * t67 - t30 * t74 + (-t32 * t67 + t42 * t77) * qJD(1), t79, -t30 * t79 + t46, t26, t17, -t35, t25, -t74, 0, t53 * t41 + (t54 * t42 + t8) * qJD(4), t53 * t42 + (-t54 * t41 - t9) * qJD(4), (-t41 * t9 - t42 * t8 + (t15 * t42 - t16 * t41) * qJD(4)) * qJD(1) + t78, -t1 * t15 + t14 * t27 - t2 * t16 + t19 * t20 + t3 * t8 + t5 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t74, t35, 0, t7 * t41 - t72 * t42 + (-t10 * t42 - t11 * t41) * qJD(4), 0, 0, 0, 0, 0, 0, -t74, t35, 0, -t50 * qJD(4) + t1 * t42 - t2 * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44, -t64, 0, 0, 0, 0, 0, 0, t22, t23, 0, -t79 - t64, 0, 0, 0, 0, 0, 0, t22, t23, 0, -t65 - t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, -t21, 0, -t29, 0, 0, -t42 * t64, t10 * qJD(4) + t41 * t64 + t72, 0, 0, t29, -t21, 0, -t29, 0, 0, t28 + (t5 - t75) * qJD(4) + (-pkin(4) * t73 - t61 * qJD(1) - t58) * t42, -t38 * t44 * pkin(4) + t4 * qJD(4) + (t61 * t41 + t57) * qJD(1) + t72, (t70 - t76) * t41 * qJD(1), t76 * t5 + (-t42 * t65 + t2) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t55, -0.2e1 * t56, (-t37 - t38) * t44, t50 * qJD(1) + t19;];
tauc_reg = t6;
