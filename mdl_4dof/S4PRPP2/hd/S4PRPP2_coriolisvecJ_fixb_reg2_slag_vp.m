% Calculate inertial parameters regressor of coriolis joint torque vector for
% S4PRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,theta3]';
% 
% Output:
% taug_reg [4x(4*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 14:01
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc_reg = S4PRPP2_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP2_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPP2_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRPP2_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:00:26
% EndTime: 2018-11-14 14:00:26
% DurationCPUTime: 0.11s
% Computational Cost: add. (89->37), mult. (245->48), div. (0->0), fcn. (178->4), ass. (0->28)
t17 = sin(pkin(5));
t18 = cos(pkin(5));
t19 = sin(qJ(2));
t20 = cos(qJ(2));
t13 = t17 * t20 + t18 * t19;
t9 = t13 * qJD(1);
t12 = t17 * t19 - t18 * t20;
t8 = t13 * qJD(2);
t6 = qJD(1) * t8;
t29 = t6 * t12;
t28 = t8 * qJD(2);
t27 = qJD(1) * t19;
t10 = t12 * qJD(2);
t26 = t10 * qJD(2);
t25 = t20 * qJD(1);
t24 = t17 * t27;
t11 = t12 * qJD(1);
t23 = -t11 + t24;
t16 = qJD(2) * pkin(2) + t25;
t5 = t17 * t16 + t18 * t27;
t4 = t18 * t16 - t24;
t21 = qJD(2) ^ 2;
t15 = t18 * qJD(2) * t25;
t7 = -qJD(2) * t24 + t15;
t3 = t15 + (qJD(4) - t24) * qJD(2);
t2 = qJD(2) * qJ(4) + t5;
t1 = -qJD(2) * pkin(3) + qJD(4) - t4;
t14 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21 * t19, -t21 * t20, 0, 0, 0, 0, 0, 0, 0, 0, -t28, t26, 0, -t5 * t10 + t7 * t13 - t4 * t8 + t29, 0, 0, 0, 0, 0, 0, -t28, 0, -t26, t1 * t8 - t2 * t10 + t3 * t13 + t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23 * qJD(2) - t15, 0, t5 * t11 + t4 * t9 + (t17 * t7 - t6 * t18) * pkin(2), 0, 0, 0, 0, 0, 0, 0, 0, t15 + (0.2e1 * qJD(4) - t23) * qJD(2), t3 * (t17 * pkin(2) + qJ(4)) + t6 * (-t18 * pkin(2) - pkin(3)) - t1 * t9 + (qJD(4) + t11) * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21 (-t2 + t9) * qJD(2);];
tauc_reg  = t14;
