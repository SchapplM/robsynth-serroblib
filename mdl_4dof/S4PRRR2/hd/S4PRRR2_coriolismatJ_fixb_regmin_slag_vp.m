% Calculate minimal parameter regressor of coriolis matrix for
% S4PRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4]';
% 
% Output:
% cmat_reg [(4*%NQJ)%x10]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 13:27
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4PRRR2_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR2_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR2_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S4PRRR2_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:27:25
% EndTime: 2019-07-18 13:27:25
% DurationCPUTime: 0.13s
% Computational Cost: add. (90->37), mult. (228->57), div. (0->0), fcn. (148->4), ass. (0->39)
t44 = -pkin(2) / 0.2e1;
t24 = cos(qJ(3));
t43 = t24 * pkin(1);
t18 = pkin(2) + t43;
t21 = sin(qJ(4));
t42 = t21 * t18;
t41 = t21 * t24;
t23 = cos(qJ(4));
t40 = t23 * t18;
t22 = sin(qJ(3));
t39 = t23 * t22;
t38 = pkin(1) * qJD(2);
t37 = pkin(1) * qJD(3);
t36 = pkin(2) * qJD(3);
t35 = pkin(2) * qJD(4);
t28 = -t43 / 0.2e1;
t25 = t28 + pkin(2) / 0.2e1 + t18 / 0.2e1;
t3 = t25 * t21;
t34 = t3 * qJD(2);
t4 = t25 * t23;
t33 = t4 * qJD(2);
t17 = t21 * t22 * pkin(1);
t7 = t17 - t40;
t32 = t7 * qJD(2);
t8 = pkin(1) * t39 + t42;
t31 = t8 * qJD(2);
t11 = (t39 + t41) * pkin(1);
t30 = t11 * qJD(2);
t12 = t23 * t43 - t17;
t29 = t12 * qJD(2);
t27 = pkin(1) * (-qJD(2) - qJD(3));
t26 = pkin(2) * (-qJD(3) - qJD(4));
t10 = t12 * qJD(3);
t9 = t11 * qJD(3);
t6 = t8 * qJD(4);
t5 = t7 * qJD(4);
t2 = t17 - t40 / 0.2e1 + (t44 + t28) * t23;
t1 = t21 * t44 - t42 / 0.2e1 + (-t39 - t41 / 0.2e1) * pkin(1);
t13 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -t22 * t37, -t24 * t37, 0, -t9 - t6, -t10 + t5; 0, 0, 0, 0, 0, t22 * t27, t24 * t27, 0, t1 * qJD(4) - t30 - t9, t2 * qJD(4) - t10 - t29; 0, 0, 0, 0, 0, 0, 0, 0, t1 * qJD(3) - t31 - t6, t2 * qJD(3) + t32 + t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, t22 * t38, t24 * t38, 0, -t3 * qJD(4) + t30, -t4 * qJD(4) + t29; 0, 0, 0, 0, 0, 0, 0, 0, -t21 * t35, -t23 * t35; 0, 0, 0, 0, 0, 0, 0, 0, t21 * t26 - t34, t23 * t26 - t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, t3 * qJD(3) + t31, t4 * qJD(3) - t32; 0, 0, 0, 0, 0, 0, 0, 0, t21 * t36 + t34, t23 * t36 + t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t13;
