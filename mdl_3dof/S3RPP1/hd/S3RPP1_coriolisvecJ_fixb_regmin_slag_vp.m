% Calculate minimal parameter regressor of coriolis joint torque vector for
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
% taug_reg [3x9]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 10:14
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc_reg = S3RPP1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RPP1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RPP1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S3RPP1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [3x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:13:44
% EndTime: 2018-11-14 10:13:44
% DurationCPUTime: 0.06s
% Computational Cost: add. (16->12), mult. (30->15), div. (0->0), fcn. (0->0), ass. (0->7)
t6 = 2 * qJD(1);
t5 = (-pkin(1) - qJ(3));
t3 = qJD(2) * t6;
t4 = qJD(1) ^ 2;
t2 = qJD(1) * qJ(2) + qJD(3);
t1 = t5 * qJD(1) + qJD(2);
t7 = [0, 0, 0, 0, t3, qJ(2) * t3, t3, qJD(3) * t6, t2 * qJD(2) - t1 * qJD(3) + (qJ(2) * qJD(2) - t5 * qJD(3)) * qJD(1); 0, 0, 0, 0, -t4, -t4 * qJ(2), -t4, 0 (-qJD(3) - t2) * qJD(1); 0, 0, 0, 0, 0, 0, 0, -t4 (qJD(2) + t1) * qJD(1);];
tauc_reg  = t7;
