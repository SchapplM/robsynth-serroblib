% Calculate inertial parameters regressor of coriolis joint torque vector for
% S4RRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
% 
% Output:
% tauc_reg [4x(4*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 18:16
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RRPR2_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR2_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR2_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPR2_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 18:16:34
% EndTime: 2019-07-18 18:16:35
% DurationCPUTime: 0.22s
% Computational Cost: add. (433->65), mult. (633->90), div. (0->0), fcn. (275->4), ass. (0->47)
t30 = qJD(1) + qJD(2);
t35 = -pkin(2) - pkin(3);
t34 = cos(qJ(2));
t53 = pkin(1) * qJD(1);
t46 = t34 * t53;
t41 = qJD(3) - t46;
t11 = t35 * t30 + t41;
t52 = pkin(1) * qJD(2);
t44 = qJD(1) * t52;
t25 = t34 * t44;
t28 = t30 * qJD(3);
t19 = t25 + t28;
t32 = sin(qJ(2));
t47 = t32 * t53;
t21 = t30 * qJ(3) + t47;
t31 = sin(qJ(4));
t33 = cos(qJ(4));
t42 = t32 * t44;
t50 = qJD(4) * t33;
t51 = qJD(4) * t31;
t1 = t11 * t50 + t33 * t19 - t21 * t51 + t31 * t42;
t29 = qJD(4) - t30;
t5 = t33 * t11 - t31 * t21;
t59 = -t5 * t29 + t1;
t2 = -t11 * t51 - t31 * t19 - t21 * t50 + t33 * t42;
t6 = t31 * t11 + t33 * t21;
t58 = t6 * t29 + t2;
t37 = -t31 * qJ(3) + t33 * t35;
t55 = t33 * qJD(3) + t37 * qJD(4) - (t31 * t32 + t33 * t34) * t53;
t38 = t33 * qJ(3) + t31 * t35;
t54 = -t31 * qJD(3) - t38 * qJD(4) - (-t31 * t34 + t32 * t33) * t53;
t49 = t32 * t52;
t48 = t34 * t52;
t45 = -t34 * pkin(1) - pkin(2);
t43 = t29 ^ 2;
t26 = -pkin(3) + t45;
t27 = t32 * pkin(1) + qJ(3);
t40 = t33 * t26 - t31 * t27;
t39 = t31 * t26 + t33 * t27;
t36 = t30 * t46 - t25;
t24 = qJD(3) + t48;
t20 = -t30 * pkin(2) + t41;
t15 = (-qJD(1) - t30) * t49;
t14 = (-qJD(2) + t30) * t47;
t4 = -t39 * qJD(4) - t31 * t24 + t33 * t49;
t3 = t40 * qJD(4) + t33 * t24 + t31 * t49;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, -t30 * t48 - t25, 0, 0, 0, 0, 0, 0, 0, 0, t15, 0, t24 * t30 + t19, t19 * t27 + t21 * t24 + (t45 * qJD(1) + t20) * t49, 0, 0, 0, 0, 0, 0, t4 * t29 - t2, -t3 * t29 + t1, 0, t1 * t39 + t2 * t40 + t6 * t3 + t5 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, t36, 0, 0, 0, 0, 0, 0, 0, 0, t14, 0, 0.2e1 * t28 - t36, t19 * qJ(3) + t21 * qJD(3) + (-t21 * t34 + (-pkin(2) * qJD(2) - t20) * t32) * t53, 0, 0, 0, 0, 0, 0, t54 * t29 - t2, -t55 * t29 + t1, 0, t1 * t38 + t2 * t37 + t54 * t5 + t55 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30 ^ 2, -t21 * t30 + t42, 0, 0, 0, 0, 0, 0, -t31 * t43, -t33 * t43, 0, t59 * t31 + t58 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, -t59, 0, 0;];
tauc_reg  = t7;
