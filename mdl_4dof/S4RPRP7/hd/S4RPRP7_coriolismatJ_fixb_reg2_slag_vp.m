% Calculate inertial parameters regressor of coriolis matrix for
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
% cmat_reg [(4*4)x(4*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4RPRP7_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP7_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP7_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP7_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:47:19
% EndTime: 2019-12-31 16:47:20
% DurationCPUTime: 0.39s
% Computational Cost: add. (105->54), mult. (203->47), div. (0->0), fcn. (150->2), ass. (0->37)
t19 = sin(qJ(3));
t31 = t19 * qJD(4);
t20 = cos(qJ(3));
t6 = -t19 * pkin(3) + t20 * qJ(4);
t22 = t6 * qJD(3) + t31;
t4 = qJ(2) - t6;
t5 = t20 * pkin(3) + t19 * qJ(4);
t1 = t5 * t19 + t4 * t20;
t40 = t1 * qJD(1);
t2 = -t4 * t19 + t5 * t20;
t39 = t2 * qJD(1);
t38 = t4 * qJD(1);
t18 = t20 ^ 2;
t9 = t19 ^ 2 - t18;
t36 = t9 * qJD(1);
t35 = t9 * qJD(3);
t34 = qJD(2) * t20;
t33 = t18 * qJD(1);
t32 = t19 * qJD(3);
t14 = t20 * qJD(1);
t13 = t20 * qJD(3);
t30 = qJ(2) * qJD(3);
t16 = qJD(1) * qJ(2);
t29 = qJD(3) * qJ(4);
t28 = t5 * t38;
t27 = t4 * t14;
t21 = -pkin(1) - pkin(5);
t26 = t21 * t32;
t25 = t21 * t13;
t24 = t19 * t16;
t23 = t20 * t16;
t17 = qJ(2) * qJD(2);
t15 = qJD(2) * t19;
t12 = t19 * qJD(1);
t11 = t19 * t13;
t10 = t19 * t14;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2), t17, -t11, t35, 0, t11, 0, 0, t20 * t30 + t15, -t19 * t30 + t34, 0, t17, -t11, 0, -t35, 0, 0, t11, t1 * qJD(3) - t20 * t31 + t15, 0, -t2 * qJD(3) + t18 * qJD(4) - t34, (qJD(3) * t5 - qJD(4) * t20 + qJD(2)) * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(1), t16, 0, 0, 0, 0, 0, 0, t12, t14, 0, t16, 0, 0, 0, 0, 0, 0, t12, 0, -t14, t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, t36, -t32, t10, -t13, 0, t23 - t26, -t24 - t25, 0, 0, -t10, -t32, -t36, 0, t13, t10, -t26 + t40, -t22, t25 - t39, t22 * t21 + t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, -t32, t33, t26 - t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(1), -t16, 0, 0, 0, 0, 0, 0, -t12, -t14, 0, -t16, 0, 0, 0, 0, 0, 0, -t12, 0, t14, -t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32, -t13, 0, 0, 0, 0, 0, 0, 0, 0, -t32, 0, t13, t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t36, 0, -t10, 0, 0, -t23, t24, 0, 0, t10, 0, t36, 0, 0, -t10, -t40, 0, t39, -t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4), qJ(4) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3), t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, 0, -t33, t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3), -t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t3;
