% Calculate minimal parameter regressor of coriolis joint torque vector for
% S3PRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d2,d3]';
% 
% Output:
% taug_reg [3x7]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 10:13
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc_reg = S3PRR2_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PRR2_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3PRR2_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3PRR2_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:12:50
% EndTime: 2018-11-14 10:12:51
% DurationCPUTime: 0.12s
% Computational Cost: add. (45->18), mult. (102->30), div. (0->0), fcn. (68->4), ass. (0->18)
t4 = qJD(2) + qJD(3);
t5 = sin(qJ(3));
t6 = sin(qJ(2));
t21 = t5 * t6;
t8 = cos(qJ(2));
t20 = t5 * t8;
t7 = cos(qJ(3));
t19 = t6 * t7;
t18 = t7 * t8;
t16 = qJD(1) * qJD(2);
t3 = qJD(2) * pkin(2) + qJD(1) * t8;
t14 = (-qJD(3) + t4) * t3;
t13 = -t19 - t20;
t12 = t18 - t21;
t11 = qJD(3) * (-pkin(2) * t4 - t3);
t10 = t13 * qJD(2);
t9 = qJD(2) ^ 2;
t1 = [0, 0, -t9 * t6, -t9 * t8, 0 (t13 * qJD(3) + t10) * t4, -t4 ^ 2 * t12; 0, 0, 0, 0, 0, t5 * t11 + (-qJD(3) * t19 - t13 * t4 + t10) * qJD(1), t7 * t11 + (-qJD(2) * t18 + (t12 + t21) * t4) * qJD(1); 0, 0, 0, 0, 0, t5 * t14 - t16 * t20 (-t8 * t16 + t14) * t7;];
tauc_reg  = t1;
