% Calculate minimal parameter regressor of coriolis joint torque vector for
% S4PRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% 
% Output:
% tauc_reg [4x15]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-14 22:36
% Revision: beb2ba9bd8c5bd556f42a244985f3dab86917626 (2021-01-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4PRRP5_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP5_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP5_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP5_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-14 22:36:08
% EndTime: 2021-01-14 22:36:10
% DurationCPUTime: 0.27s
% Computational Cost: add. (215->76), mult. (581->118), div. (0->0), fcn. (301->4), ass. (0->65)
t24 = sin(qJ(2));
t27 = qJD(3) ^ 2;
t28 = qJD(2) ^ 2;
t66 = (t27 + t28) * t24;
t49 = t24 * qJD(1);
t19 = qJD(2) * t49;
t23 = sin(qJ(3));
t46 = qJD(2) * qJD(3);
t40 = t23 * t46;
t11 = pkin(3) * t40 + t19;
t25 = cos(qJ(3));
t52 = qJD(3) * pkin(3);
t16 = qJD(2) * pkin(5) + t49;
t38 = qJ(4) * qJD(2) + t16;
t6 = t38 * t23;
t5 = -t6 + t52;
t7 = t38 * t25;
t34 = t23 * t5 - t25 * t7;
t65 = t34 * qJD(2) + t11;
t64 = t5 + t6;
t21 = t23 ^ 2;
t63 = pkin(3) * t21;
t62 = t25 * pkin(3);
t61 = t25 * t16;
t60 = t27 * t23;
t59 = t27 * t25;
t26 = cos(qJ(2));
t58 = t28 * t26;
t57 = qJ(4) + pkin(5);
t22 = t25 ^ 2;
t56 = t21 - t22;
t55 = t21 + t22;
t53 = qJD(2) * pkin(2);
t20 = -pkin(2) - t62;
t48 = t26 * qJD(1);
t10 = qJD(2) * t20 + qJD(4) - t48;
t51 = qJD(2) * t10;
t50 = qJD(3) * t23;
t47 = qJ(4) * qJD(3);
t45 = t23 * t28 * t25;
t44 = 0.2e1 * t46;
t43 = t23 * t47;
t42 = t25 * t47;
t41 = qJD(3) * t48;
t39 = qJD(3) * t57;
t37 = -0.2e1 * t26 * t46;
t36 = t25 * t44;
t35 = qJD(4) + t48;
t33 = -t11 + t19;
t32 = -t10 - t35;
t17 = -t48 - t53;
t31 = qJD(3) * (t17 - t53);
t30 = qJD(2) * (-t17 - t48);
t1 = -t16 * t50 + (t35 * t25 - t43) * qJD(2);
t2 = -qJD(3) * t61 + (-t35 * t23 - t42) * qJD(2);
t29 = t1 * t25 - t2 * t23 + (-t23 * t7 - t25 * t5) * qJD(3);
t15 = t25 * t41;
t14 = t23 * t41;
t13 = t57 * t25;
t12 = t57 * t23;
t9 = -t23 * qJD(4) - t25 * t39;
t8 = t25 * qJD(4) - t23 * t39;
t4 = t23 * t37 - t25 * t66;
t3 = t23 * t66 + t25 * t37;
t18 = [0, 0, -t28 * t24, -t58, 0, 0, 0, 0, 0, t4, t3, t4, t3, t55 * t58, -t65 * t26 + (t29 + t51) * t24; 0, 0, 0, 0, t23 * t36, -t56 * t44, t59, -t60, 0, -pkin(5) * t59 + t23 * t31 + t14, pkin(5) * t60 + t25 * t31 + t15, t14 + t33 * t25 + (t9 + (t10 + (t20 - t62) * qJD(2)) * t23) * qJD(3), t15 - t33 * t23 + (t10 * t25 - t8 + (t20 * t25 + t63) * qJD(2)) * qJD(3), (-t23 * t9 + t25 * t8 + (t12 * t25 - t13 * t23) * qJD(3) - t55 * t48) * qJD(2) + t29, t10 * pkin(3) * t50 + t1 * t13 + t11 * t20 - t2 * t12 + t5 * t9 + t7 * t8 + (-t10 * t24 + t34 * t26) * qJD(1); 0, 0, 0, 0, -t45, t56 * t28, 0, 0, 0, t23 * t30, t25 * t30, pkin(3) * t45 + (t7 - t61) * qJD(3) + (t32 * t23 - t42) * qJD(2), -t28 * t63 + (t23 * t16 - t6) * qJD(3) + (t32 * t25 + t43) * qJD(2), (-t52 + t64) * t25 * qJD(2), t64 * t7 + (-t23 * t51 + t2) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t40, t36, -t55 * t28, t65;];
tauc_reg = t18;
