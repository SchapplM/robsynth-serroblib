% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RPPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
% 
% Output:
% tauc_reg [5x17]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPPRR5_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR5_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR5_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR5_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:56:38
% EndTime: 2019-12-31 17:56:39
% DurationCPUTime: 0.22s
% Computational Cost: add. (301->44), mult. (521->75), div. (0->0), fcn. (245->6), ass. (0->39)
t30 = sin(qJ(4));
t33 = qJD(5) ^ 2;
t45 = qJD(1) - qJD(4);
t58 = t45 ^ 2;
t59 = t30 * (t33 + t58);
t47 = qJD(5) * t45;
t17 = -cos(pkin(8)) * pkin(1) - pkin(2) - pkin(3);
t15 = t17 * qJD(1) + qJD(3);
t18 = sin(pkin(8)) * pkin(1) + qJ(3);
t16 = qJD(1) * t18;
t46 = qJD(3) * qJD(1);
t32 = cos(qJ(4));
t48 = qJD(4) * t32;
t49 = qJD(4) * t30;
t2 = t15 * t49 + t16 * t48 + t30 * t46;
t56 = -t2 - (t30 * t15 + t32 * t16) * t45;
t40 = t30 * t17 + t32 * t18;
t55 = t2 + (t30 * qJD(3) + qJD(4) * t40) * t45;
t54 = t45 * pkin(4);
t29 = sin(qJ(5));
t31 = cos(qJ(5));
t51 = t29 * t31;
t50 = t29 ^ 2 - t31 ^ 2;
t5 = t32 * t15 - t30 * t16;
t3 = -t5 + t54;
t38 = -t15 * t48 + t16 * t49 - t32 * t46;
t44 = t3 * t45 + t38;
t42 = t47 * t51;
t41 = t32 * t17 - t30 * t18;
t39 = pkin(7) * t33 - t56;
t37 = (-pkin(7) + t40) * t33 - t55;
t36 = qJD(5) * (t3 + t5 + t54);
t7 = t32 * qJD(3) + qJD(4) * t41;
t35 = qJD(5) * (-t45 * (pkin(4) - t41) - t3 - t7);
t34 = 0.2e1 * t32 * t47;
t22 = t33 * t31;
t21 = t33 * t29;
t13 = t50 * t47;
t1 = [0, 0, 0, 0, 0, 0.2e1 * t46, 0.2e1 * t16 * qJD(3), 0, t55, t45 * t7 - t38, 0.2e1 * t42, -0.2e1 * t13, -t22, t21, 0, t29 * t35 - t31 * t37, t29 * t37 + t31 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, t22; 0, 0, 0, 0, 0, -qJD(1) ^ 2, -t16 * qJD(1), 0, -t30 * t58, -t32 * t58, 0, 0, 0, 0, 0, t29 * t34 - t31 * t59, t29 * t59 + t31 * t34; 0, 0, 0, 0, 0, 0, 0, 0, t56, -t45 * t5 + t38, -0.2e1 * t42, 0.2e1 * t13, t22, -t21, 0, t29 * t36 - t31 * t39, t29 * t39 + t31 * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t58 * t51, t50 * t58, 0, 0, 0, t44 * t29, t44 * t31;];
tauc_reg = t1;
