% Calculate minimal parameter regressor of coriolis matrix for
% S4RPRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
% 
% Output:
% cmat_reg [(4*%NQJ)%x17]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4RPRP7_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP7_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP7_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP7_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:47:18
% EndTime: 2019-12-31 16:47:19
% DurationCPUTime: 0.17s
% Computational Cost: add. (96->46), mult. (171->48), div. (0->0), fcn. (124->2), ass. (0->35)
t13 = sin(qJ(3));
t26 = t13 * qJD(4);
t14 = cos(qJ(3));
t6 = -t13 * pkin(3) + t14 * qJ(4);
t16 = t6 * qJD(3) + t26;
t4 = qJ(2) - t6;
t5 = t14 * pkin(3) + t13 * qJ(4);
t1 = t5 * t13 + t4 * t14;
t36 = t1 * qJD(1);
t2 = -t4 * t13 + t5 * t14;
t35 = t2 * qJD(1);
t34 = t4 * qJD(1);
t12 = t14 ^ 2;
t7 = t13 ^ 2 - t12;
t32 = t7 * qJD(1);
t31 = qJD(2) * t14;
t30 = qJD(3) * t13;
t29 = qJD(3) * t14;
t15 = -pkin(1) - pkin(5);
t28 = qJD(3) * t15;
t27 = t12 * qJD(1);
t10 = t14 * qJD(1);
t25 = qJ(2) * qJD(3);
t24 = qJD(1) * qJ(2);
t23 = qJD(3) * qJ(4);
t22 = t5 * t34;
t21 = t4 * t10;
t20 = t13 * t28;
t19 = t14 * t28;
t18 = t13 * t24;
t17 = t14 * t24;
t11 = qJD(2) * t13;
t9 = t13 * qJD(1);
t8 = t13 * t10;
t3 = [0, 0, 0, 0, qJD(2), qJ(2) * qJD(2), -t13 * t29, t7 * qJD(3), 0, 0, 0, t14 * t25 + t11, -t13 * t25 + t31, t1 * qJD(3) - t14 * t26 + t11, 0, -t2 * qJD(3) + t12 * qJD(4) - t31, (qJD(3) * t5 - qJD(4) * t14 + qJD(2)) * t4; 0, 0, 0, 0, qJD(1), t24, 0, 0, 0, 0, 0, t9, t10, t9, 0, -t10, t34; 0, 0, 0, 0, 0, 0, -t8, t32, -t30, -t29, 0, t17 - t20, -t18 - t19, -t20 + t36, -t16, t19 - t35, t16 * t15 + t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, -t30, t27, t20 - t21; 0, 0, 0, 0, -qJD(1), -t24, 0, 0, 0, 0, 0, -t9, -t10, -t9, 0, t10, -t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30, -t29, -t30, 0, t29, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30; 0, 0, 0, 0, 0, 0, t8, -t32, 0, 0, 0, -t17, t18, -t36, 0, t35, -t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4), qJ(4) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3), t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, 0, -t27, t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3), -t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t3;
