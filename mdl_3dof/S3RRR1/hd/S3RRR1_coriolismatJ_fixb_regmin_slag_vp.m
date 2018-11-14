% Calculate minimal parameter regressor of coriolis matrix for
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
% cmat_reg [(3*%NQJ)%x9]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 10:16
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function cmat_reg = S3RRR1_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RRR1_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RRR1_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S3RRR1_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:15:59
% EndTime: 2018-11-14 10:16:00
% DurationCPUTime: 0.11s
% Computational Cost: add. (90->37), mult. (228->57), div. (0->0), fcn. (148->4), ass. (0->39)
t44 = -pkin(2) / 0.2e1;
t24 = cos(qJ(2));
t43 = t24 * pkin(1);
t18 = pkin(2) + t43;
t21 = sin(qJ(3));
t42 = t21 * t18;
t41 = t21 * t24;
t23 = cos(qJ(3));
t40 = t23 * t18;
t22 = sin(qJ(2));
t39 = t23 * t22;
t38 = pkin(1) * qJD(1);
t37 = pkin(1) * qJD(2);
t36 = pkin(2) * qJD(2);
t35 = pkin(2) * qJD(3);
t28 = -t43 / 0.2e1;
t25 = t28 + pkin(2) / 0.2e1 + t18 / 0.2e1;
t3 = t25 * t21;
t34 = t3 * qJD(1);
t4 = t25 * t23;
t33 = t4 * qJD(1);
t17 = t21 * t22 * pkin(1);
t7 = t17 - t40;
t32 = t7 * qJD(1);
t8 = pkin(1) * t39 + t42;
t31 = t8 * qJD(1);
t11 = (t39 + t41) * pkin(1);
t30 = t11 * qJD(1);
t12 = t23 * t43 - t17;
t29 = t12 * qJD(1);
t27 = pkin(1) * (-qJD(1) - qJD(2));
t26 = pkin(2) * (-qJD(2) - qJD(3));
t10 = t12 * qJD(2);
t9 = t11 * qJD(2);
t6 = t8 * qJD(3);
t5 = t7 * qJD(3);
t2 = t17 - t40 / 0.2e1 + (t44 + t28) * t23;
t1 = t21 * t44 - t42 / 0.2e1 + (-t39 - t41 / 0.2e1) * pkin(1);
t13 = [0, 0, 0, 0, -t22 * t37, -t24 * t37, 0, -t9 - t6, -t10 + t5; 0, 0, 0, 0, t22 * t27, t24 * t27, 0, t1 * qJD(3) - t30 - t9, t2 * qJD(3) - t10 - t29; 0, 0, 0, 0, 0, 0, 0, t1 * qJD(2) - t31 - t6, t2 * qJD(2) + t32 + t5; 0, 0, 0, 0, t22 * t38, t24 * t38, 0, -t3 * qJD(3) + t30, -t4 * qJD(3) + t29; 0, 0, 0, 0, 0, 0, 0, -t21 * t35, -t23 * t35; 0, 0, 0, 0, 0, 0, 0, t21 * t26 - t34, t23 * t26 - t33; 0, 0, 0, 0, 0, 0, 0, t3 * qJD(2) + t31, t4 * qJD(2) - t32; 0, 0, 0, 0, 0, 0, 0, t21 * t36 + t34, t23 * t36 + t33; 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t13;
