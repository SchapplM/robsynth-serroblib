% Calculate inertial parameters regressor of coriolis joint torque vector for
% S3RPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1]';
% 
% Output:
% tauc_reg [3x(3*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S3RPP1_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RPP1_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RPP1_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S3RPP1_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [3x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:05:01
% EndTime: 2019-03-08 18:05:02
% DurationCPUTime: 0.05s
% Computational Cost: add. (16->12), mult. (30->15), div. (0->0), fcn. (0->0), ass. (0->7)
t6 = 2 * qJD(1);
t5 = (-pkin(1) - qJ(3));
t3 = qJD(2) * t6;
t4 = qJD(1) ^ 2;
t2 = qJ(2) * qJD(1) + qJD(3);
t1 = qJD(1) * t5 + qJD(2);
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, qJ(2) * t3, 0, 0, 0, 0, 0, 0, 0, t3, qJD(3) * t6, t2 * qJD(2) - t1 * qJD(3) + (qJ(2) * qJD(2) - t5 * qJD(3)) * qJD(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4, -t4 * qJ(2), 0, 0, 0, 0, 0, 0, 0, -t4, 0 (-qJD(3) - t2) * qJD(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4 (qJD(2) + t1) * qJD(1);];
tauc_reg  = t7;
