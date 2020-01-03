% Calculate minimal parameter regressor of coriolis joint torque vector for
% S4RPRP4
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
% tauc_reg [4x15]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RPRP4_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP4_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP4_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP4_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:43:53
% EndTime: 2019-12-31 16:43:54
% DurationCPUTime: 0.18s
% Computational Cost: add. (202->53), mult. (483->86), div. (0->0), fcn. (226->4), ass. (0->41)
t25 = cos(qJ(3));
t15 = sin(pkin(6)) * pkin(1) + pkin(5);
t13 = t15 * qJD(1);
t24 = sin(qJ(3));
t45 = t24 * t13;
t7 = t25 * qJD(2) - t45;
t47 = qJD(4) - t7;
t46 = 2 * qJD(3);
t2 = -(qJD(3) * pkin(3)) + t47;
t8 = t24 * qJD(2) + t25 * t13;
t3 = (qJD(3) * qJ(4)) + t8;
t44 = t24 * t25;
t26 = qJD(3) ^ 2;
t43 = t26 * t24;
t19 = t26 * t25;
t38 = qJD(2) * qJD(3);
t40 = qJD(3) * t25;
t6 = t13 * t40 + t24 * t38;
t20 = t24 ^ 2;
t42 = -t25 ^ 2 + t20;
t31 = pkin(3) * t24 - qJ(4) * t25;
t9 = t31 * qJD(3) - t24 * qJD(4);
t5 = qJD(1) * t9;
t16 = -cos(pkin(6)) * pkin(1) - pkin(2);
t10 = -t25 * pkin(3) - t24 * qJ(4) + t16;
t4 = qJD(1) * t10;
t14 = qJD(1) * t16;
t41 = qJD(1) * t24;
t27 = qJD(1) ^ 2;
t37 = t27 * t44;
t36 = qJD(1) * t46;
t35 = t7 + t45;
t34 = 0.2e1 * t4;
t32 = t8 * qJD(3) - t6;
t30 = t14 * t46;
t29 = -t15 * t26 - 0.2e1 * t5;
t18 = t25 * t38;
t1 = t18 + (qJD(4) - t45) * qJD(3);
t28 = t1 * t25 + t6 * t24 + (t2 * t25 - t24 * t3) * qJD(3);
t12 = t31 * qJD(1);
t11 = [0, 0, 0, 0, t36 * t44, -t42 * t36, t19, -t43, 0, -t15 * t19 + t24 * t30, t15 * t43 + t25 * t30, t34 * t24 * qJD(3) + t29 * t25, t28, t29 * t24 - t34 * t40, t5 * t10 + t28 * t15 + t4 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43, -t19, -t43, 0, t19, t1 * t24 - t6 * t25 + (t2 * t24 + t25 * t3) * qJD(3); 0, 0, 0, 0, -t37, t42 * t27, 0, 0, 0, -t14 * t41 + t32, -t14 * t25 * qJD(1) + t35 * qJD(3) - t18, (t12 * t25 - t24 * t4) * qJD(1) + t32, 0, t18 + (t12 * t24 + t25 * t4) * qJD(1) + (0.2e1 * qJD(4) - t35) * qJD(3), -t6 * pkin(3) + t1 * qJ(4) - t4 * t12 - t2 * t8 + t47 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, 0, -t20 * t27 - t26, -t3 * qJD(3) + t4 * t41 + t6;];
tauc_reg = t11;
