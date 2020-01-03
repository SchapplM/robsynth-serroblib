% Calculate minimal parameter regressor of coriolis joint torque vector for
% S4RPRR2
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
% tauc_reg [4x14]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RPRR2_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR2_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR2_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR2_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:48:13
% EndTime: 2019-12-31 16:48:13
% DurationCPUTime: 0.14s
% Computational Cost: add. (161->37), mult. (385->52), div. (0->0), fcn. (198->6), ass. (0->38)
t18 = cos(pkin(7)) * pkin(1) + pkin(2);
t17 = t18 * qJD(1);
t27 = sin(qJ(3));
t29 = cos(qJ(3));
t50 = pkin(1) * sin(pkin(7));
t40 = qJD(1) * t50;
t8 = t29 * t17 - t27 * t40;
t9 = t27 * t17 + t29 * t40;
t52 = t9 * qJD(3);
t26 = sin(qJ(4));
t28 = cos(qJ(4));
t21 = qJD(1) + qJD(3);
t49 = t21 * pkin(3);
t4 = -t8 - t49;
t42 = qJD(4) * t4;
t51 = t26 * t52 + t28 * t42;
t48 = t9 * t21;
t33 = t27 * t18 + t29 * t50;
t47 = t33 * qJD(3) * t21;
t46 = t26 * t28;
t30 = qJD(4) ^ 2;
t44 = t30 * t26;
t43 = t26 ^ 2 - t28 ^ 2;
t41 = 0.2e1 * qJD(4) * t21;
t6 = t8 * qJD(3);
t39 = -t4 * t21 - t6;
t37 = pkin(6) * t30 - t48;
t36 = qJD(4) * (t8 - t49);
t35 = (pkin(6) + t33) * t30 + t47;
t32 = t29 * t18 - t27 * t50;
t10 = t32 * qJD(3);
t34 = qJD(4) * ((-pkin(3) - t32) * t21 - t10);
t20 = t21 ^ 2;
t19 = t30 * t28;
t15 = t41 * t46;
t12 = t43 * t41;
t1 = t26 * t42;
t2 = [0, 0, 0, 0, 0, -t52 - t47, -t10 * t21 - t6, t15, -t12, t19, -t44, 0, t1 + t26 * t34 + (-t35 - t52) * t28, t35 * t26 + t28 * t34 + t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44, -t19; 0, 0, 0, 0, 0, -t52 + t48, t8 * t21 - t6, t15, -t12, t19, -t44, 0, t1 + t26 * t36 + (-t37 - t52) * t28, t37 * t26 + t28 * t36 + t51; 0, 0, 0, 0, 0, 0, 0, -t20 * t46, t43 * t20, 0, 0, 0, t39 * t26, t39 * t28;];
tauc_reg = t2;
