% Calculate minimal parameter regressor of coriolis joint torque vector for
% S4RPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
% 
% Output:
% tauc_reg [4x18]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RPRP5_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP5_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP5_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP5_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:45:01
% EndTime: 2019-12-31 16:45:02
% DurationCPUTime: 0.32s
% Computational Cost: add. (442->99), mult. (1236->129), div. (0->0), fcn. (842->4), ass. (0->63)
t47 = sin(pkin(6));
t48 = cos(pkin(6));
t49 = sin(qJ(3));
t50 = cos(qJ(3));
t31 = t50 * t47 + t49 * t48;
t26 = t31 * qJD(1);
t77 = t26 ^ 2;
t72 = t50 * t48;
t62 = qJD(1) * t72;
t73 = t49 * t47;
t63 = qJD(1) * t73;
t24 = -t62 + t63;
t44 = -t48 * pkin(2) - pkin(1);
t34 = t44 * qJD(1) + qJD(2);
t5 = t24 * pkin(3) - t26 * qJ(4) + t34;
t76 = t5 * t26;
t75 = t26 * t24;
t71 = pkin(5) + qJ(2);
t36 = t71 * t48;
t33 = qJD(1) * t36;
t74 = t49 * t33;
t70 = t47 ^ 2 + t48 ^ 2;
t30 = -t72 + t73;
t35 = t71 * t47;
t56 = -t50 * t35 - t49 * t36;
t6 = -t30 * qJD(2) + t56 * qJD(3);
t69 = t6 * qJD(3);
t17 = -t49 * t35 + t50 * t36;
t7 = t31 * qJD(2) + t17 * qJD(3);
t68 = t7 * qJD(3);
t67 = qJD(3) * t49;
t66 = qJD(3) * t50;
t32 = qJD(1) * t35;
t14 = -t50 * t32 - t74;
t65 = qJD(4) - t14;
t64 = qJD(1) * qJD(2);
t61 = t70 * qJD(1) ^ 2;
t60 = t49 * t64;
t59 = t50 * t64;
t58 = t14 + t74;
t3 = -t32 * t67 + t33 * t66 + t47 * t59 + t48 * t60;
t40 = qJD(3) * t62;
t18 = qJD(3) * t63 - t40;
t29 = t31 * qJD(3);
t19 = qJD(1) * t29;
t57 = t19 * pkin(3) + t18 * qJ(4);
t15 = -t49 * t32 + t50 * t33;
t55 = 0.2e1 * t70 * t64;
t54 = t15 * qJD(3) - t3;
t53 = t32 * t66 + t47 * t60 - t48 * t59;
t52 = 0.2e1 * t26 * qJD(3);
t28 = t47 * t67 - t48 * t66;
t23 = t24 ^ 2;
t13 = t30 * pkin(3) - t31 * qJ(4) + t44;
t12 = t26 * pkin(3) + t24 * qJ(4);
t11 = qJD(3) * qJ(4) + t15;
t10 = t40 + (t24 - t63) * qJD(3);
t9 = -t40 + (t24 + t63) * qJD(3);
t8 = -qJD(3) * pkin(3) + t65;
t4 = t29 * pkin(3) + t28 * qJ(4) - t31 * qJD(4);
t2 = (qJD(4) - t74) * qJD(3) - t53;
t1 = -t26 * qJD(4) + t57;
t16 = [0, 0, 0, 0, 0, t55, qJ(2) * t55, -t18 * t31 - t26 * t28, t18 * t30 - t31 * t19 + t28 * t24 - t26 * t29, -t28 * qJD(3), -t29 * qJD(3), 0, t44 * t19 + t34 * t29 - t68, -t44 * t18 - t34 * t28 - t69, t1 * t30 + t13 * t19 + t4 * t24 + t5 * t29 - t68, -t11 * t29 - t17 * t19 + t18 * t56 - t2 * t30 - t6 * t24 + t7 * t26 - t8 * t28 + t3 * t31, -t1 * t31 + t13 * t18 - t4 * t26 + t5 * t28 + t69, t1 * t13 + t11 * t6 + t2 * t17 - t3 * t56 + t5 * t4 + t8 * t7; 0, 0, 0, 0, 0, -t61, -qJ(2) * t61, 0, 0, 0, 0, 0, t52, -t9, t52, -t23 - t77, t9, t11 * t24 + (-qJD(4) - t8) * t26 + t57; 0, 0, 0, 0, 0, 0, 0, t75, -t23 + t77, t10, 0, 0, -t34 * t26 + t54, t58 * qJD(3) + t34 * t24 + t53, -t12 * t24 + t54 - t76, pkin(3) * t18 - t19 * qJ(4) + (t11 - t15) * t26 + (t8 - t65) * t24, t12 * t26 - t5 * t24 + (0.2e1 * qJD(4) - t58) * qJD(3) - t53, -t3 * pkin(3) + t2 * qJ(4) + t65 * t11 - t5 * t12 - t8 * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, t10, -qJD(3) ^ 2 - t77, -t11 * qJD(3) + t3 + t76;];
tauc_reg = t16;
