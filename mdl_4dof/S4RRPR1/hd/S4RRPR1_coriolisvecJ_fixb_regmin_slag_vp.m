% Calculate minimal parameter regressor of coriolis joint torque vector for
% S4RRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% taug_reg [4x10]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:54
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc_reg = S4RRPR1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:53:29
% EndTime: 2018-11-14 13:53:29
% DurationCPUTime: 0.19s
% Computational Cost: add. (188->53), mult. (482->85), div. (0->0), fcn. (274->6), ass. (0->37)
t20 = qJD(1) + qJD(2);
t19 = qJD(4) + t20;
t39 = qJD(4) - t19;
t21 = sin(pkin(7));
t38 = pkin(2) * t21;
t24 = sin(qJ(2));
t37 = t21 * t24;
t22 = cos(pkin(7));
t36 = t22 * t24;
t35 = pkin(1) * qJD(1);
t34 = t24 * t35;
t23 = sin(qJ(4));
t25 = cos(qJ(4));
t26 = cos(qJ(2));
t28 = pkin(1) * (-t21 * t26 - t36);
t10 = qJD(2) * t28;
t6 = qJD(1) * t10;
t27 = pkin(1) * (t22 * t26 - t37);
t12 = qJD(2) * t27;
t7 = qJD(1) * t12;
t33 = -t23 * t7 + t25 * t6;
t18 = t26 * pkin(1) + pkin(2);
t32 = -pkin(1) * t37 + t22 * t18;
t15 = t20 * pkin(2) + t26 * t35;
t3 = t22 * t15 - t21 * t34;
t1 = t20 * pkin(3) + t3;
t4 = t21 * t15 + t22 * t34;
t31 = -t23 * t1 - t25 * t4;
t30 = (-qJD(2) + t20) * t35;
t29 = pkin(1) * qJD(2) * (-qJD(1) - t20);
t17 = t22 * pkin(2) + pkin(3);
t13 = pkin(1) * t36 + t21 * t18;
t11 = qJD(1) * t27;
t9 = qJD(1) * t28;
t8 = pkin(3) + t32;
t2 = qJD(4) * t23 * t4;
t5 = [0, 0, 0, 0, t24 * t29, t26 * t29, t3 * t10 + t4 * t12 + t7 * t13 + t6 * t32, 0 (t25 * t10 - t23 * t12) * t19 + ((-t13 * t25 - t23 * t8) * t19 + t31) * qJD(4) + t33, t2 + (-(-qJD(4) * t13 + t10) * t19 - t6) * t23 + (-(qJD(4) * t8 + t12) * t19 - t7 - qJD(4) * t1) * t25; 0, 0, 0, 0, t24 * t30, t26 * t30, -t4 * t11 - t3 * t9 + (t21 * t7 + t22 * t6) * pkin(2), 0 -(-t23 * t11 + t25 * t9) * t19 + ((-t17 * t23 - t25 * t38) * t19 + t31) * qJD(4) + t33, t2 - t25 * t7 - t23 * t6 + (t25 * t11 + t23 * t9) * t19 + (-(t17 * t25 - t23 * t38) * t19 - t25 * t1) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, t39 * t31 + t33, t2 + (-t4 * t19 - t6) * t23 + (-t1 * t39 - t7) * t25;];
tauc_reg  = t5;
