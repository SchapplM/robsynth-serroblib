% Calculate minimal parameter regressor of coriolis matrix for
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
% cmat_reg [(3*%NQJ)%x9]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S3RPR1_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RPR1_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RPR1_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3RPR1_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:05:57
% EndTime: 2019-03-08 18:05:57
% DurationCPUTime: 0.06s
% Computational Cost: add. (30->15), mult. (43->18), div. (0->0), fcn. (32->2), ass. (0->14)
t14 = -pkin(1) - pkin(2);
t6 = sin(qJ(3));
t13 = t6 * qJD(1);
t12 = t6 * qJD(2);
t7 = cos(qJ(3));
t11 = t7 * qJD(1);
t10 = t7 * qJD(2);
t9 = qJ(2) * qJD(1);
t8 = qJD(1) - qJD(3);
t4 = t8 * t7;
t3 = t8 * t6;
t2 = t7 * qJ(2) + t6 * t14;
t1 = t6 * qJ(2) - t7 * t14;
t5 = [0, 0, 0, 0, qJD(2), qJ(2) * qJD(2), 0, t2 * qJD(3) + t12, -t1 * qJD(3) + t10; 0, 0, 0, 0, qJD(1), t9, 0, t13, t11; 0, 0, 0, 0, 0, 0, 0, t8 * t2, -t8 * t1; 0, 0, 0, 0, -qJD(1), -t9, 0, -t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, t3, t4; 0, 0, 0, 0, 0, 0, 0, -t2 * qJD(1) - t12, t1 * qJD(1) - t10; 0, 0, 0, 0, 0, 0, 0, -t13, -t11; 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t5;
