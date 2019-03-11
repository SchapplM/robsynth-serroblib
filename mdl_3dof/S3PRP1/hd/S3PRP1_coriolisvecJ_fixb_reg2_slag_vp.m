% Calculate inertial parameters regressor of coriolis joint torque vector for
% S3PRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d2]';
% 
% Output:
% tauc_reg [3x(3*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S3PRP1_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PRP1_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3PRP1_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S3PRP1_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [3x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:03:01
% EndTime: 2019-03-08 18:03:01
% DurationCPUTime: 0.06s
% Computational Cost: add. (22->15), mult. (49->19), div. (0->0), fcn. (21->2), ass. (0->13)
t4 = sin(qJ(2));
t8 = t4 * qJD(1);
t3 = qJD(2) * qJ(3) + t8;
t5 = cos(qJ(2));
t12 = t3 * t5;
t6 = qJD(2) ^ 2;
t11 = t6 * t4;
t10 = t6 * t5;
t9 = qJD(2) * pkin(2);
t7 = t5 * qJD(1);
t2 = qJD(3) - t7 - t9;
t1 = (qJD(3) + t7) * qJD(2);
t13 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, -t10, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, t10, t1 * t4 + (t12 + (t2 - t7) * t4) * qJD(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * qJD(2) * qJD(3), t1 * qJ(3) + t3 * qJD(3) + (-t12 + (-t2 - t9) * t4) * qJD(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6 (-t3 + t8) * qJD(2);];
tauc_reg  = t13;
