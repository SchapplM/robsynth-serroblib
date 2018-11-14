% Calculate inertial parameters regressor of coriolis joint torque vector for
% S4PRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta3]';
% 
% Output:
% taug_reg [4x(4*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 14:12
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc_reg = S4PRPR3_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR3_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR3_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR3_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:11:25
% EndTime: 2018-11-14 14:11:26
% DurationCPUTime: 0.20s
% Computational Cost: add. (269->47), mult. (674->79), div. (0->0), fcn. (528->6), ass. (0->37)
t28 = sin(pkin(6));
t29 = cos(pkin(6));
t31 = sin(qJ(2));
t33 = cos(qJ(2));
t38 = t28 * t31 - t29 * t33;
t21 = t38 * qJD(1);
t25 = qJD(2) * pkin(2) + t33 * qJD(1);
t43 = qJD(1) * t31;
t12 = t29 * t25 - t28 * t43;
t10 = qJD(2) * pkin(3) + t12;
t23 = t28 * t33 + t29 * t31;
t18 = t23 * qJD(2);
t16 = qJD(1) * t18;
t20 = t38 * qJD(2);
t17 = qJD(1) * t20;
t30 = sin(qJ(4));
t32 = cos(qJ(4));
t13 = t28 * t25 + t29 * t43;
t44 = t30 * t13;
t1 = (qJD(4) * t10 - t17) * t32 - qJD(4) * t44 - t30 * t16;
t47 = pkin(2) * t28;
t19 = t23 * qJD(1);
t26 = t29 * pkin(2) + pkin(3);
t36 = t32 * t26 - t30 * t47;
t46 = t36 * qJD(4) + t30 * t19 + t32 * t21;
t37 = t30 * t26 + t32 * t47;
t45 = -t37 * qJD(4) + t32 * t19 - t30 * t21;
t6 = t30 * t10 + t32 * t13;
t40 = -t30 * t23 - t32 * t38;
t39 = t32 * t23 - t30 * t38;
t2 = -t6 * qJD(4) - t32 * t16 + t30 * t17;
t34 = qJD(2) ^ 2;
t27 = qJD(2) + qJD(4);
t5 = t32 * t10 - t44;
t4 = -t39 * qJD(4) - t32 * t18 + t30 * t20;
t3 = t40 * qJD(4) - t30 * t18 - t32 * t20;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34 * t31, -t34 * t33, 0, 0, 0, 0, 0, 0, 0, 0, -t18 * qJD(2), t20 * qJD(2), 0, -t12 * t18 - t13 * t20 + t16 * t38 - t17 * t23, 0, 0, 0, 0, 0, 0, t4 * t27, -t3 * t27, 0, t1 * t39 + t2 * t40 + t6 * t3 + t5 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12 * t19 + t13 * t21 + (-t16 * t29 - t17 * t28) * pkin(2), 0, 0, 0, 0, 0, 0, t45 * t27 + t2, -t46 * t27 - t1, 0, t1 * t37 + t2 * t36 + t45 * t5 + t46 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6 * t27 + t2, t5 * t27 - t1, 0, 0;];
tauc_reg  = t7;
