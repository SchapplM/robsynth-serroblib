% Calculate minimal parameter regressor of coriolis matrix for
% S5PPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d4,d5,theta1,theta2,theta3]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x13]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 14:59
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5PPPRR2_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR2_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPPRR2_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR2_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 14:59:39
% EndTime: 2019-12-05 14:59:40
% DurationCPUTime: 0.18s
% Computational Cost: add. (74->29), mult. (245->60), div. (0->0), fcn. (256->8), ass. (0->30)
t12 = sin(pkin(9));
t13 = sin(pkin(8));
t35 = t12 * t13;
t19 = cos(qJ(4));
t34 = t12 * t19;
t14 = cos(pkin(9));
t33 = t13 * t14;
t18 = cos(qJ(5));
t32 = qJD(4) * t18;
t15 = cos(pkin(8));
t17 = sin(qJ(4));
t10 = -t15 * t17 + t19 * t33;
t31 = t10 * qJD(4);
t16 = sin(qJ(5));
t11 = -t16 ^ 2 + t18 ^ 2;
t30 = t11 * qJD(4);
t29 = t16 * qJD(5);
t28 = t17 * qJD(4);
t27 = t18 * qJD(5);
t26 = t19 * qJD(4);
t25 = pkin(4) * t16 * qJD(4);
t24 = pkin(4) * t32;
t23 = t16 * t32;
t22 = t12 * t28;
t21 = t17 * t29 - t18 * t26;
t20 = t16 * t26 + t17 * t27;
t9 = t15 * t19 + t17 * t33;
t4 = t9 * t18;
t2 = t9 * t16;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, -t31, t9 * qJD(4), 0, 0, 0, 0, 0, t2 * qJD(5) - t18 * t31, t4 * qJD(5) + t16 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2 * qJD(4) + (-t10 * t18 - t16 * t35) * qJD(5), t4 * qJD(4) + (t10 * t16 - t18 * t35) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, -t12 * t26, t22, 0, 0, 0, 0, 0, t21 * t12, t20 * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16 * t22 + (t16 * t14 - t18 * t34) * qJD(5), t18 * t22 + (t18 * t14 + t16 * t34) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, -t28, -t26, 0, 0, 0, 0, 0, -t18 * t28 - t19 * t29, t16 * t28 - t19 * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20, t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t16 * t27, t11 * qJD(5), 0, 0, 0, -pkin(4) * t29, -pkin(4) * t27; 0, 0, 0, 0, 0, 0, t23, t30, t27, -t29, 0, -pkin(6) * t27 - t25, pkin(6) * t29 - t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, -t23, -t30, 0, 0, 0, t25, t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t1;
