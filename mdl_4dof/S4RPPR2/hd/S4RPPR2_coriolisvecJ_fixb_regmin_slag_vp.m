% Calculate minimal parameter regressor of coriolis joint torque vector for
% S4RPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta3]';
% 
% Output:
% taug_reg [4x12]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:48
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc_reg = S4RPPR2_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR2_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR2_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR2_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:47:32
% EndTime: 2018-11-14 13:47:32
% DurationCPUTime: 0.14s
% Computational Cost: add. (114->34), mult. (214->57), div. (0->0), fcn. (110->4), ass. (0->27)
t9 = -qJD(1) + qJD(4);
t28 = qJD(4) - t9;
t14 = -pkin(1) - pkin(2);
t10 = sin(pkin(6));
t27 = t10 * qJ(2);
t26 = qJD(1) * qJ(2);
t25 = qJD(1) * qJD(2);
t24 = 0.2e1 * t25;
t11 = cos(pkin(6));
t23 = t11 * t14 - t27;
t7 = t14 * qJD(1) + qJD(2);
t6 = t11 * t7;
t1 = t6 + (-pkin(3) - t27) * qJD(1);
t12 = sin(qJ(4));
t13 = cos(qJ(4));
t3 = t10 * t7 + t11 * t26;
t22 = t13 * t1 - t12 * t3;
t21 = -t12 * t1 - t13 * t3;
t20 = (-t10 * t26 + t6) * t10 - t3 * t11;
t19 = t10 * t13 + t11 * t12;
t18 = t10 * t12 - t11 * t13;
t17 = t19 * qJD(1);
t16 = t18 * qJD(1);
t15 = qJD(1) ^ 2;
t5 = t11 * qJ(2) + t10 * t14;
t4 = -pkin(3) + t23;
t2 = [0, 0, 0, 0, t24, qJ(2) * t24, t10 * t24, t11 * t24 ((-t10 * t23 + t11 * t5) * qJD(1) - t20) * qJD(2), 0 ((-t12 * t4 - t13 * t5) * t9 - t21) * qJD(4) + (-t19 * t9 + t17) * qJD(2) (-(-t12 * t5 + t13 * t4) * t9 + t22) * qJD(4) + (t18 * t9 - t16) * qJD(2); 0, 0, 0, 0, -t15, -t15 * qJ(2), -t10 * t15, -t11 * t15, t20 * qJD(1), 0 (-t19 * qJD(4) + t17) * t9 (t18 * qJD(4) - t16) * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19 * t25 + t28 * t21, t18 * t25 - t28 * t22;];
tauc_reg  = t2;
