% Calculate minimal parameter regressor of coriolis matrix for
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
% cmat_reg [(3*%NQJ)%x7]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S3PRP1_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PRP1_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3PRP1_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S3PRP1_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [3x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:03:01
% EndTime: 2019-03-08 18:03:01
% DurationCPUTime: 0.05s
% Computational Cost: add. (7->7), mult. (12->8), div. (0->0), fcn. (8->2), ass. (0->6)
t1 = sin(qJ(2));
t5 = t1 * qJD(2);
t2 = cos(qJ(2));
t4 = t2 * qJD(2);
t3 = qJD(2) * qJ(3);
t6 = [0, 0, 0, 0, 0, 0, 0; 0, 0, -t5, -t4, -t5, t4 (-t1 * pkin(2) + t2 * qJ(3)) * qJD(2) + t1 * qJD(3); 0, 0, 0, 0, 0, 0, t5; 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, qJD(3), qJ(3) * qJD(3); 0, 0, 0, 0, 0, qJD(2), t3; 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -qJD(2), -t3; 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t6;
