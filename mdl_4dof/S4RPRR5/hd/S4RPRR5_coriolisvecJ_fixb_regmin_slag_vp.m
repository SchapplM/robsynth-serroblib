% Calculate minimal parameter regressor of coriolis joint torque vector for
% S4RPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
% 
% Output:
% tauc_reg [4x16]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RPRR5_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR5_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR5_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR5_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:51:40
% EndTime: 2019-12-31 16:51:40
% DurationCPUTime: 0.19s
% Computational Cost: add. (197->41), mult. (349->70), div. (0->0), fcn. (152->4), ass. (0->40)
t22 = sin(qJ(3));
t26 = qJD(4) ^ 2;
t40 = qJD(1) - qJD(3);
t56 = t40 ^ 2;
t57 = t22 * (t26 + t56);
t25 = -pkin(1) - pkin(2);
t13 = t25 * qJD(1) + qJD(2);
t41 = (qJD(1) * qJD(2));
t54 = qJD(3) * t13 + t41;
t43 = qJD(4) * t40;
t24 = cos(qJ(3));
t42 = qJ(2) * qJD(1);
t37 = qJD(3) * t42;
t2 = t54 * t22 + t24 * t37;
t53 = -t2 - (t22 * t13 + t24 * t42) * t40;
t34 = t24 * qJ(2) + t22 * t25;
t52 = t2 + (t22 * qJD(2) + t34 * qJD(3)) * t40;
t51 = t40 * pkin(3);
t21 = sin(qJ(4));
t23 = cos(qJ(4));
t48 = t21 * t23;
t47 = t26 * t21;
t46 = t26 * t23;
t45 = t21 ^ 2 - t23 ^ 2;
t39 = 2 * t41;
t1 = -t22 * t37 + t54 * t24;
t6 = t24 * t13 - t22 * t42;
t3 = -t6 + t51;
t38 = t3 * t40 - t1;
t35 = t43 * t48;
t33 = -t22 * qJ(2) + t24 * t25;
t32 = pkin(6) * t26 - t53;
t31 = (-pkin(6) + t34) * t26 - t52;
t30 = qJD(4) * (t3 + t6 + t51);
t4 = t24 * qJD(2) + t33 * qJD(3);
t29 = qJD(4) * (-(pkin(3) - t33) * t40 - t3 - t4);
t28 = 0.2e1 * t24 * t43;
t27 = qJD(1) ^ 2;
t8 = t45 * t43;
t5 = [0, 0, 0, 0, t39, qJ(2) * t39, 0, t52, t4 * t40 + t1, 0.2e1 * t35, -0.2e1 * t8, -t46, t47, 0, t21 * t29 - t31 * t23, t31 * t21 + t23 * t29; 0, 0, 0, 0, -t27, -t27 * qJ(2), 0, -t22 * t56, -t24 * t56, 0, 0, 0, 0, 0, t21 * t28 - t23 * t57, t21 * t57 + t23 * t28; 0, 0, 0, 0, 0, 0, 0, t53, -t40 * t6 - t1, -0.2e1 * t35, 0.2e1 * t8, t46, -t47, 0, t21 * t30 - t32 * t23, t32 * t21 + t23 * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t56 * t48, t45 * t56, 0, 0, 0, t38 * t21, t38 * t23;];
tauc_reg = t5;
