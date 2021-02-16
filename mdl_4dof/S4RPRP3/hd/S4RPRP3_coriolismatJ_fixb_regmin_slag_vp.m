% Calculate minimal parameter regressor of coriolis matrix for
% S4RPRP3
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
% cmat_reg [(4*%NQJ)%x15]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 10:20
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4RPRP3_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP3_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP3_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP3_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 10:20:38
% EndTime: 2021-01-15 10:20:39
% DurationCPUTime: 0.15s
% Computational Cost: add. (160->35), mult. (300->47), div. (0->0), fcn. (231->4), ass. (0->36)
t23 = cos(qJ(3));
t42 = t23 * pkin(3);
t17 = -cos(pkin(6)) * pkin(1) - pkin(2);
t13 = t17 - t42;
t22 = sin(qJ(3));
t41 = t13 * t22;
t1 = pkin(3) * t41;
t40 = t1 * qJD(1);
t16 = sin(pkin(6)) * pkin(1) + pkin(5);
t36 = qJ(4) + t16;
t10 = t36 * t22;
t11 = t36 * t23;
t5 = t10 * t22 + t11 * t23;
t39 = t5 * qJD(1);
t6 = t22 * t42 - t41;
t38 = t6 * qJD(1);
t19 = t22 ^ 2;
t9 = t19 * pkin(3) + t13 * t23;
t37 = t9 * qJD(1);
t35 = qJD(3) * t22;
t18 = qJD(3) * t23;
t34 = t11 * qJD(3);
t20 = t23 ^ 2;
t14 = t19 + t20;
t33 = t14 * qJD(1);
t15 = t20 - t19;
t32 = t15 * qJD(1);
t31 = t22 * qJD(1);
t30 = t22 * qJD(4);
t29 = t23 * qJD(1);
t28 = pkin(3) * t35;
t27 = pkin(3) * t31;
t26 = t17 * t31;
t25 = t17 * t29;
t24 = t22 * t29;
t2 = [0, 0, 0, 0, t22 * t18, t15 * qJD(3), 0, 0, 0, t17 * t35, t17 * t18, -t6 * qJD(3), t9 * qJD(3), t14 * qJD(4), t1 * qJD(3) + t5 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, t24, t32, t18, -t35, 0, -t16 * t18 + t26, t16 * t35 + t25, -t34 - t38, t10 * qJD(3) + t37, -pkin(3) * t18, -pkin(3) * t34 + t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35, -t18, -t35, -t18, 0, -t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, -t24, -t32, 0, 0, 0, -t26, -t25, -t30 + t38, -t23 * qJD(4) - t37, 0, -pkin(3) * t30 - t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, -t29, 0, -t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, t18, -t33, t28 - t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, t29, 0, t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t2;
