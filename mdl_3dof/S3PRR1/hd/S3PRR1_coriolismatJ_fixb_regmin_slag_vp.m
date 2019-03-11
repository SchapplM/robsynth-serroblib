% Calculate minimal parameter regressor of coriolis matrix for
% S3PRR1
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
% cmat_reg [(3*%NQJ)%x7]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S3PRR1_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PRR1_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3PRR1_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3PRR1_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:04:01
% EndTime: 2019-03-08 18:04:01
% DurationCPUTime: 0.05s
% Computational Cost: add. (18->9), mult. (42->17), div. (0->0), fcn. (42->4), ass. (0->11)
t12 = pkin(2) * qJD(3);
t11 = qJD(2) * pkin(2);
t10 = qJD(2) + qJD(3);
t9 = pkin(2) * t10;
t8 = cos(qJ(2));
t7 = cos(qJ(3));
t6 = sin(qJ(2));
t5 = sin(qJ(3));
t2 = t10 * (t5 * t6 - t7 * t8);
t1 = t10 * (-t5 * t8 - t7 * t6);
t3 = [0, 0, 0, 0, 0, 0, 0; 0, 0, -t6 * qJD(2), -t8 * qJD(2), 0, t1, t2; 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -t5 * t12, -t7 * t12; 0, 0, 0, 0, 0, -t5 * t9, -t7 * t9; 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, t5 * t11, t7 * t11; 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t3;
