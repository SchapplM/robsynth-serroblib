% Calculate inertial parameters regressor of coriolis joint torque vector for
% S3RRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2]';
% 
% Output:
% taug_reg [3x(3*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 10:15
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc_reg = S3RRP1_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RRP1_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RRP1_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3RRP1_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:15:09
% EndTime: 2018-11-14 10:15:09
% DurationCPUTime: 0.08s
% Computational Cost: add. (53->23), mult. (107->32), div. (0->0), fcn. (33->2), ass. (0->20)
t11 = cos(qJ(2));
t18 = pkin(1) * qJD(2);
t13 = qJD(1) * t18;
t7 = t11 * t13;
t9 = qJD(1) + qJD(2);
t8 = t9 * qJD(3);
t3 = t7 + t8;
t19 = pkin(1) * qJD(1);
t10 = sin(qJ(2));
t17 = t10 * t18;
t16 = t11 * t18;
t15 = t10 * t19;
t14 = t11 * t19;
t12 = t9 * t14 - t7;
t6 = qJD(3) + t16;
t5 = t9 * qJ(3) + t15;
t4 = -t9 * pkin(2) + qJD(3) - t14;
t2 = (-qJD(1) - t9) * t17;
t1 = (-qJD(2) + t9) * t15;
t20 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, -t9 * t16 - t7, 0, 0, 0, 0, 0, 0, 0, 0, t2, 0, t6 * t9 + t3, t3 * (t10 * pkin(1) + qJ(3)) + t5 * t6 + (t4 + (-t11 * pkin(1) - pkin(2)) * qJD(1)) * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t12, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t12 + 0.2e1 * t8, t3 * qJ(3) + t5 * qJD(3) + (-t11 * t5 + (-pkin(2) * qJD(2) - t4) * t10) * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9 ^ 2, t10 * t13 - t5 * t9;];
tauc_reg  = t20;
