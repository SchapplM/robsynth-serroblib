% Calculate inertial parameters regressor of coriolis joint torque vector for
% S4RRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% 
% Output:
% tauc_reg [4x(4*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RRRR1_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR1_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR1_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR1_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:22:14
% EndTime: 2019-12-31 17:22:15
% DurationCPUTime: 0.38s
% Computational Cost: add. (639->88), mult. (1466->130), div. (0->0), fcn. (752->6), ass. (0->79)
t32 = qJD(1) + qJD(2);
t40 = cos(qJ(2));
t76 = pkin(1) * qJD(1);
t69 = t40 * t76;
t21 = t32 * pkin(2) + t69;
t39 = cos(qJ(3));
t60 = qJD(2) * t69;
t36 = sin(qJ(3));
t37 = sin(qJ(2));
t71 = t37 * t76;
t62 = t36 * t71;
t80 = (qJD(2) + qJD(3)) * t62;
t4 = (qJD(3) * t21 + t60) * t39 - t80;
t35 = sin(qJ(4));
t33 = t35 ^ 2;
t38 = cos(qJ(4));
t34 = t38 ^ 2;
t77 = t33 + t34;
t93 = t77 * t4;
t83 = t36 * t37;
t53 = t39 * t40 - t83;
t18 = t53 * t76;
t74 = qJD(3) * t39;
t92 = pkin(2) * t74 - t18;
t26 = t36 * pkin(2) + pkin(7);
t30 = qJD(3) + t32;
t41 = qJD(4) ^ 2;
t82 = t37 * t39;
t54 = t36 * t40 + t82;
t17 = t54 * t76;
t75 = qJD(3) * t36;
t59 = pkin(2) * t75 - t17;
t91 = t26 * t41 + t59 * t30;
t90 = t54 * qJD(2) + t37 * t74;
t12 = t39 * t21 - t62;
t87 = t30 * pkin(3);
t10 = -t12 - t87;
t44 = t90 * pkin(1);
t68 = t21 * t75;
t5 = qJD(1) * t44 + t68;
t73 = qJD(4) * t38;
t88 = t10 * t73 + t5 * t35;
t28 = t40 * pkin(1) + pkin(2);
t7 = t28 * t75 + t44;
t86 = t7 * t30;
t13 = t36 * t21 + t39 * t71;
t85 = t13 * t30;
t81 = t41 * t35;
t79 = pkin(1) * t82 + t36 * t28;
t78 = t33 - t34;
t29 = t30 ^ 2;
t72 = t35 * t29 * t38;
t6 = t28 * t74 + (t53 * qJD(2) - t37 * t75) * pkin(1);
t65 = t77 * t6;
t64 = -t10 * t30 - t4;
t63 = t77 * t12;
t61 = t35 * t30 * t73;
t58 = (-qJD(2) + t32) * t76;
t57 = pkin(1) * qJD(2) * (-qJD(1) - t32);
t56 = pkin(7) * t41 - t85;
t16 = pkin(7) + t79;
t55 = t16 * t41 + t86;
t52 = qJD(4) * (t12 - t87);
t49 = -pkin(1) * t83 + t39 * t28;
t15 = -pkin(3) - t49;
t51 = qJD(4) * (t15 * t30 - t6);
t50 = (-pkin(2) * t30 - t21) * qJD(3);
t27 = -t39 * pkin(2) - pkin(3);
t47 = qJD(4) * (t27 * t30 - t92);
t46 = t92 * t77;
t43 = t90 * t76;
t42 = -t43 - t68;
t31 = t41 * t38;
t20 = -0.2e1 * t61;
t19 = 0.2e1 * t61;
t14 = -0.2e1 * t78 * t30 * qJD(4);
t11 = t30 * pkin(7) + t13;
t8 = t10 * qJD(4) * t35;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37 * t57, t40 * t57, 0, 0, 0, 0, 0, 0, 0, 0, t42 - t86, -t6 * t30 - t4, 0, -t12 * t7 + t13 * t6 + t4 * t79 - t5 * t49, t19, t14, t31, t20, -t81, 0, t8 + t35 * t51 + (-t5 - t55) * t38, t55 * t35 + t38 * t51 + t88, t30 * t65 + t93, t10 * t7 + t11 * t65 + t5 * t15 + t16 * t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37 * t58, t40 * t58, 0, 0, 0, 0, 0, 0, 0, 0, t17 * t30 + t36 * t50 - t43, t18 * t30 + (t50 - t60) * t39 + t80, 0, t12 * t17 - t13 * t18 + (t36 * t4 - t39 * t5 + (-t12 * t36 + t13 * t39) * qJD(3)) * pkin(2), t19, t14, t31, t20, -t81, 0, t8 + t35 * t47 + (-t5 - t91) * t38, t91 * t35 + t38 * t47 + t88, t46 * t30 + t93, t59 * t10 + t46 * t11 + t26 * t93 + t5 * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42 + t85, t12 * t30 - t4, 0, 0, t19, t14, t31, t20, -t81, 0, t8 + t35 * t52 + (-t5 - t56) * t38, t56 * t35 + t38 * t52 + t88, -t30 * t63 + t93, -t5 * pkin(3) + pkin(7) * t93 - t10 * t13 - t11 * t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t72, t78 * t29, 0, t72, 0, 0, t64 * t35, t64 * t38, 0, 0;];
tauc_reg = t1;
