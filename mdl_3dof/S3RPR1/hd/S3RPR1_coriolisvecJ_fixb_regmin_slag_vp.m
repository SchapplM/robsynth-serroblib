% Calculate minimal parameter regressor of coriolis joint torque vector for
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
% tauc_reg [3x9]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S3RPR1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RPR1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RPR1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3RPR1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:05:57
% EndTime: 2019-03-08 18:05:57
% DurationCPUTime: 0.11s
% Computational Cost: add. (44->21), mult. (76->33), div. (0->0), fcn. (26->2), ass. (0->15)
t5 = -pkin(1) - pkin(2);
t1 = t5 * qJD(1) + qJD(2);
t3 = sin(qJ(3));
t15 = t3 * t1;
t14 = qJ(2) * t3;
t11 = qJD(1) - qJD(3);
t13 = qJD(3) + t11;
t12 = qJD(1) * qJ(2);
t10 = 0.2e1 * qJD(1) * qJD(2);
t4 = cos(qJ(3));
t9 = t13 * t4;
t8 = qJD(2) * (qJD(1) + t11);
t7 = t11 ^ 2;
t6 = qJD(1) ^ 2;
t2 = [0, 0, 0, 0, t10, qJ(2) * t10, 0, t3 * t8 + (-(-qJ(2) * t4 - t3 * t5) * t11 + t4 * t12 + t15) * qJD(3), t4 * t8 + ((t4 * t5 - t14) * t11 - t3 * t12 + t4 * t1) * qJD(3); 0, 0, 0, 0, -t6, -t6 * qJ(2), 0, -t3 * t7, -t4 * t7; 0, 0, 0, 0, 0, 0, 0, -t13 * t15 + (-qJ(2) * t9 - t3 * qJD(2)) * qJD(1), -t1 * t9 + (-qJD(2) * t4 + t13 * t14) * qJD(1);];
tauc_reg  = t2;
