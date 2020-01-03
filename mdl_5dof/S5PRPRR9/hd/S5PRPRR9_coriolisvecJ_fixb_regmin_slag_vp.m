% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5PRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
% 
% Output:
% tauc_reg [5x17]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRPRR9_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR9_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR9_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR9_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:39:45
% EndTime: 2019-12-31 17:39:46
% DurationCPUTime: 0.21s
% Computational Cost: add. (225->41), mult. (391->70), div. (0->0), fcn. (170->4), ass. (0->40)
t25 = sin(qJ(4));
t29 = qJD(5) ^ 2;
t43 = qJD(2) - qJD(4);
t57 = t43 ^ 2;
t58 = t25 * (t29 + t57);
t28 = -pkin(2) - pkin(3);
t14 = t28 * qJD(2) + qJD(3);
t44 = (qJD(2) * qJD(3));
t55 = qJD(4) * t14 + t44;
t46 = qJD(5) * t43;
t27 = cos(qJ(4));
t45 = qJ(3) * qJD(2);
t40 = qJD(4) * t45;
t2 = t55 * t25 + t27 * t40;
t54 = -t2 - (t25 * t14 + t27 * t45) * t43;
t37 = t27 * qJ(3) + t25 * t28;
t53 = t2 + (t25 * qJD(3) + t37 * qJD(4)) * t43;
t52 = t43 * pkin(4);
t24 = sin(qJ(5));
t26 = cos(qJ(5));
t49 = t24 * t26;
t48 = t24 ^ 2 - t26 ^ 2;
t42 = 2 * t44;
t1 = -t25 * t40 + t55 * t27;
t7 = t27 * t14 - t25 * t45;
t3 = -t7 + t52;
t41 = t3 * t43 - t1;
t38 = t46 * t49;
t36 = -t25 * qJ(3) + t27 * t28;
t35 = pkin(7) * t29 - t54;
t34 = (-pkin(7) + t37) * t29 - t53;
t33 = qJD(5) * (t3 + t7 + t52);
t5 = t27 * qJD(3) + t36 * qJD(4);
t32 = qJD(5) * (-(pkin(4) - t36) * t43 - t3 - t5);
t31 = 0.2e1 * t27 * t46;
t30 = qJD(2) ^ 2;
t19 = t29 * t26;
t18 = t29 * t24;
t9 = t48 * t46;
t4 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, t19; 0, 0, 0, 0, 0, t42, qJ(3) * t42, 0, t53, t43 * t5 + t1, 0.2e1 * t38, -0.2e1 * t9, -t19, t18, 0, t24 * t32 - t34 * t26, t34 * t24 + t26 * t32; 0, 0, 0, 0, 0, -t30, -t30 * qJ(3), 0, -t25 * t57, -t27 * t57, 0, 0, 0, 0, 0, t24 * t31 - t26 * t58, t24 * t58 + t26 * t31; 0, 0, 0, 0, 0, 0, 0, 0, t54, -t43 * t7 - t1, -0.2e1 * t38, 0.2e1 * t9, t19, -t18, 0, t24 * t33 - t35 * t26, t35 * t24 + t26 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t57 * t49, t48 * t57, 0, 0, 0, t41 * t24, t41 * t26;];
tauc_reg = t4;
