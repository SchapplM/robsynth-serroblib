% Calculate inertial parameters regressor of coriolis joint torque vector for
% S3RPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d3]';
% 
% Output:
% taug_reg [3x(3*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 10:14
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc_reg = S3RPR1_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RPR1_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RPR1_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3RPR1_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:14:28
% EndTime: 2018-11-14 10:14:29
% DurationCPUTime: 0.10s
% Computational Cost: add. (109->26), mult. (192->39), div. (0->0), fcn. (74->2), ass. (0->25)
t12 = cos(qJ(3));
t11 = sin(qJ(3));
t24 = qJ(2) * qJD(1);
t19 = t11 * t24;
t23 = qJD(1) * qJD(2);
t26 = qJD(3) * t12;
t13 = -pkin(1) - pkin(2);
t8 = t13 * qJD(1) + qJD(2);
t1 = -qJD(3) * t19 + t12 * t23 + t8 * t26;
t22 = qJD(1) - qJD(3);
t5 = t12 * t8 - t19;
t29 = t22 * t5 + t1;
t25 = t11 * qJD(2);
t27 = t11 * t8;
t2 = -(qJ(2) * t26 + t25) * qJD(1) - qJD(3) * t27;
t6 = t12 * t24 + t27;
t28 = -t6 * t22 + t2;
t20 = 0.2e1 * t23;
t18 = t22 ^ 2;
t17 = t12 * qJ(2) + t11 * t13;
t16 = -t11 * qJ(2) + t12 * t13;
t14 = qJD(1) ^ 2;
t4 = -t17 * qJD(3) - t25;
t3 = t12 * qJD(2) + t16 * qJD(3);
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, qJ(2) * t20, 0, 0, 0, 0, 0, 0, -t22 * t4 - t2, t22 * t3 + t1, 0, t1 * t17 + t2 * t16 + t6 * t3 + t5 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, -t14 * qJ(2), 0, 0, 0, 0, 0, 0, -t11 * t18, -t12 * t18, 0, t29 * t11 + t28 * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, -t29, 0, 0;];
tauc_reg  = t7;
