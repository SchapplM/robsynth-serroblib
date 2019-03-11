% Calculate inertial parameters regressor of coriolis matrix for
% S3RRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2,d3]';
% 
% Output:
% cmat_reg [(3*3)x(3*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S3RRR1_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RRR1_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RRR1_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S3RRR1_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:08:09
% EndTime: 2019-03-08 18:08:09
% DurationCPUTime: 0.16s
% Computational Cost: add. (177->41), mult. (445->68), div. (0->0), fcn. (314->4), ass. (0->39)
t26 = sin(qJ(3));
t7 = t26 * pkin(2);
t29 = cos(qJ(2));
t49 = t29 * pkin(1);
t27 = sin(qJ(2));
t48 = t26 * t27;
t47 = t26 * t29;
t28 = cos(qJ(3));
t46 = t28 * t27;
t45 = t28 * t29;
t44 = pkin(1) * qJD(1);
t43 = pkin(1) * qJD(2);
t42 = pkin(2) * qJD(2);
t41 = pkin(2) * qJD(3);
t33 = pkin(2) + t49;
t21 = t28 * t33;
t13 = pkin(1) * t48 - t21;
t30 = t26 * t33;
t14 = pkin(1) * t46 + t30;
t17 = (t46 + t47) * pkin(1);
t18 = (t45 - t48) * pkin(1);
t2 = t13 * t17 + t14 * t18;
t40 = t2 * qJD(1);
t39 = t7 * qJD(1);
t8 = t21 / 0.2e1 + (-t49 / 0.2e1 + pkin(2) / 0.2e1) * t28;
t38 = t8 * qJD(1);
t37 = t13 * qJD(1);
t36 = t14 * qJD(1);
t35 = t17 * qJD(1);
t34 = t18 * qJD(1);
t32 = pkin(1) * (-qJD(1) - qJD(2));
t31 = pkin(2) * (-qJD(2) - qJD(3));
t16 = t18 * qJD(2);
t15 = t17 * qJD(2);
t12 = t14 * qJD(3);
t11 = t13 * qJD(3);
t6 = -t28 * pkin(2) / 0.2e1 - t21 / 0.2e1 + (t48 - t45 / 0.2e1) * pkin(1);
t5 = -t7 / 0.2e1 - t30 / 0.2e1 + (-t46 - t47 / 0.2e1) * pkin(1);
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t27 * t43, -t29 * t43, 0, 0, 0, 0, 0, 0, 0, 0, -t15 - t12, -t16 + t11, 0, t2 * qJD(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27 * t32, t29 * t32, 0, 0, 0, 0, 0, 0, 0, 0, t5 * qJD(3) - t15 - t35, t6 * qJD(3) - t16 - t34, 0, t40 + (-t17 * t28 + t18 * t26) * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5 * qJD(2) - t12 - t36, t6 * qJD(2) + t11 + t37, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27 * t44, t29 * t44, 0, 0, 0, 0, 0, 0, 0, 0, -t7 * qJD(3) + t35, -t8 * qJD(3) + t34, 0, -t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26 * t41, -t28 * t41, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26 * t31 - t39, t28 * t31 - t38, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7 * qJD(2) + t36, t8 * qJD(2) - t37, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26 * t42 + t39, t28 * t42 + t38, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t1;
