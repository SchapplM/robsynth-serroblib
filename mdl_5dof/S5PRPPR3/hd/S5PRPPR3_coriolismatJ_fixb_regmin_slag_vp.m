% Calculate minimal parameter regressor of coriolis matrix for
% S5PRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x15]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:27
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5PRPPR3_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR3_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR3_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR3_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:26:57
% EndTime: 2019-12-05 15:26:58
% DurationCPUTime: 0.19s
% Computational Cost: add. (90->33), mult. (187->46), div. (0->0), fcn. (217->6), ass. (0->23)
t25 = cos(qJ(2));
t11 = sin(qJ(5));
t13 = cos(qJ(5));
t6 = t11 ^ 2 - t13 ^ 2;
t24 = t6 * qJD(2);
t10 = sin(pkin(8));
t8 = t10 * pkin(2) + qJ(4);
t23 = t8 * qJD(2);
t22 = cos(pkin(8));
t21 = qJD(5) * t11;
t20 = qJD(5) * t13;
t19 = t11 * qJD(2);
t18 = t13 * qJD(2);
t17 = t8 * t19;
t16 = t8 * t18;
t15 = t11 * t18;
t14 = -t22 * pkin(2) - pkin(3);
t12 = sin(qJ(2));
t7 = -pkin(6) + t14;
t5 = t10 * t25 + t22 * t12;
t4 = t10 * t12 - t22 * t25;
t3 = t5 * qJD(2);
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -qJD(2) * t12, -qJD(2) * t25, (-t10 * t4 - t22 * t5) * qJD(2) * pkin(2), t3, -t4 * qJD(2), (t5 * t14 - t4 * t8) * qJD(2) + t5 * qJD(4), 0, 0, 0, 0, 0, -t4 * t19 + t5 * t20, -t4 * t18 - t5 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, t3, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5 * t18 - t4 * t21, -t5 * t19 - t4 * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, qJD(4), t8 * qJD(4), -t11 * t20, t6 * qJD(5), 0, 0, 0, qJD(4) * t11 + t8 * t20, qJD(4) * t13 - t8 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, qJD(2), t23, 0, 0, 0, 0, 0, t19, t18; 0, 0, 0, 0, 0, 0, 0, 0, -t15, t24, -t21, -t20, 0, -t7 * t21 + t16, -t7 * t20 - t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20, t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, -qJD(2), -t23, 0, 0, 0, 0, 0, -t19, -t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, -t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, t15, -t24, 0, 0, 0, -t16, t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t1;
