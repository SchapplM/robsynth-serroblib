% Calculate minimal parameter regressor of coriolis joint torque vector for
% S4RRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% taug_reg [4x10]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:55
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc_reg = S4RRRP1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:54:31
% EndTime: 2018-11-14 13:54:32
% DurationCPUTime: 0.19s
% Computational Cost: add. (171->43), mult. (410->67), div. (0->0), fcn. (206->4), ass. (0->36)
t16 = qJD(1) + qJD(2);
t20 = cos(qJ(2));
t36 = pkin(1) * qJD(1);
t32 = t20 * t36;
t11 = t16 * pkin(2) + t32;
t19 = cos(qJ(3));
t30 = qJD(2) * t32;
t17 = sin(qJ(3));
t18 = sin(qJ(2));
t33 = t18 * t36;
t31 = t17 * t33;
t37 = (qJD(2) + qJD(3)) * t31;
t1 = (qJD(3) * t11 + t30) * t19 - t37;
t38 = t18 * t19;
t27 = -t17 * t20 - t38;
t34 = qJD(3) * t19;
t22 = (t27 * qJD(2) - t18 * t34) * pkin(1);
t21 = qJD(1) * t22;
t35 = qJD(3) * t17;
t2 = -t11 * t35 + t21;
t40 = t2 * pkin(3);
t39 = t17 * t18;
t29 = (-qJD(2) + t16) * t36;
t28 = pkin(1) * qJD(2) * (-qJD(1) - t16);
t26 = t19 * t20 - t39;
t15 = qJD(3) + t16;
t25 = (-pkin(2) * t15 - t11) * qJD(3);
t6 = t19 * t11 - t31;
t14 = t20 * pkin(1) + pkin(2);
t9 = t26 * t36;
t8 = t27 * t36;
t7 = t17 * t11 + t19 * t33;
t5 = t15 * pkin(3) + t6;
t4 = -t14 * t35 + t22;
t3 = t14 * t34 + (t26 * qJD(2) - t18 * t35) * pkin(1);
t10 = [0, 0, 0, 0, t18 * t28, t20 * t28, 0, t4 * t15 + t2, -t3 * t15 - t1, t1 * (pkin(1) * t38 + t17 * t14) + t7 * t3 + t2 * (-pkin(1) * t39 + t19 * t14 + pkin(3)) + t5 * t4; 0, 0, 0, 0, t18 * t29, t20 * t29, 0, -t8 * t15 + t17 * t25 + t21, t9 * t15 + (t25 - t30) * t19 + t37, t40 - t5 * t8 - t7 * t9 + (t1 * t17 + t2 * t19 + (-t17 * t5 + t19 * t7) * qJD(3)) * pkin(2); 0, 0, 0, 0, 0, 0, 0, t7 * t15 + t2, t6 * t15 - t1, t40 + (t5 - t6) * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
tauc_reg  = t10;
