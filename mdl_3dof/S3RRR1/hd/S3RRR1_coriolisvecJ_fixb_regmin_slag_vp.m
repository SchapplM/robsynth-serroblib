% Calculate minimal parameter regressor of coriolis joint torque vector for
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
% taug_reg [3x9]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 10:16
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc_reg = S3RRR1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RRR1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RRR1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S3RRR1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:15:59
% EndTime: 2018-11-14 10:15:59
% DurationCPUTime: 0.15s
% Computational Cost: add. (78->30), mult. (188->50), div. (0->0), fcn. (90->4), ass. (0->26)
t21 = qJD(2) + qJD(3);
t7 = sin(qJ(3));
t8 = sin(qJ(2));
t28 = t7 * t8;
t9 = cos(qJ(3));
t27 = t8 * t9;
t25 = pkin(1) * qJD(1);
t18 = t25 * t28;
t26 = t21 * t18;
t24 = qJD(3) * t8;
t10 = cos(qJ(2));
t23 = qJD(2) * t10;
t22 = t10 * qJD(1);
t20 = pkin(1) * t22;
t6 = qJD(1) + qJD(2);
t1 = t6 * pkin(2) + t20;
t5 = qJD(1) + t21;
t19 = (-qJD(3) + t5) * t1;
t17 = -t10 * t7 - t27;
t16 = (-qJD(2) + t6) * t25;
t15 = pkin(1) * qJD(2) * (-qJD(1) - t6);
t14 = qJD(3) * (-pkin(2) * t5 - t1);
t13 = qJD(3) * (-(t10 * pkin(1) + pkin(2)) * t5 - t1);
t12 = (t10 * t9 - t28) * t5;
t11 = t17 * qJD(2) - t9 * t24;
t2 = [0, 0, 0, 0, t8 * t15, t10 * t15, 0, t7 * t13 + (qJD(1) + t5) * pkin(1) * t11, t9 * t13 + (t7 * t5 * t24 + (-t9 * t22 - t12) * qJD(2)) * pkin(1) + t26; 0, 0, 0, 0, t8 * t16, t10 * t16, 0, t7 * t14 + (-t17 * t5 + t11) * t25, t9 * t14 + (-t9 * t23 + t12) * t25 + t26; 0, 0, 0, 0, 0, 0, 0, t7 * t19 + (-t7 * t23 + (t5 - t21) * t27) * t25, -t5 * t18 + (-qJD(2) * t20 + t19) * t9 + t26;];
tauc_reg  = t2;
