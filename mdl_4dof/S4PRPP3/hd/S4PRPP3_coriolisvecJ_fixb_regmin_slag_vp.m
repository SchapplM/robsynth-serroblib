% Calculate minimal parameter regressor of coriolis joint torque vector for
% S4PRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2]';
% 
% Output:
% taug_reg [4x10]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 14:02
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc_reg = S4PRPP3_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP3_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPP3_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S4PRPP3_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:01:25
% EndTime: 2018-11-14 14:01:25
% DurationCPUTime: 0.08s
% Computational Cost: add. (45->25), mult. (96->23), div. (0->0), fcn. (40->2), ass. (0->21)
t12 = sin(qJ(2));
t13 = cos(qJ(2));
t18 = t12 * qJD(1);
t9 = qJD(2) * qJ(3) + t18;
t22 = t9 * t13;
t17 = t13 * qJD(1);
t7 = (qJD(3) + t17) * qJD(2);
t24 = qJD(2) * t22 + t7 * t12;
t23 = t7 * qJ(3) + t9 * qJD(3);
t14 = qJD(2) ^ 2;
t21 = t14 * t12;
t20 = qJD(2) * pkin(2);
t19 = qJD(2) * t12;
t16 = (-pkin(2) - pkin(3)) * qJD(2);
t15 = qJD(3) - t17;
t11 = 0.2e1 * qJD(2) * qJD(3);
t10 = t14 * t13;
t8 = t15 - t20;
t4 = t16 + t15;
t1 = (-t9 + t18) * qJD(2);
t2 = [0, 0, -t21, -t10, -t21, t10 (t8 - t17) * t19 + t24, -t21, t10 (t4 - t17) * t19 + t24; 0, 0, 0, 0, 0, t11 (-t22 + (-t8 - t20) * t12) * qJD(1) + t23, 0, t11 (-t22 + (t16 - t4) * t12) * qJD(1) + t23; 0, 0, 0, 0, 0, -t14, t1, 0, -t14, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
tauc_reg  = t2;
