% Calculate inertial parameters regressor of coriolis matrix for
% S3RRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2]';
% 
% Output:
% cmat_reg [(3*3)x(3*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 10:15
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function cmat_reg = S3RRP1_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RRP1_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RRP1_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3RRP1_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:15:09
% EndTime: 2018-11-14 10:15:10
% DurationCPUTime: 0.09s
% Computational Cost: add. (37->21), mult. (75->22), div. (0->0), fcn. (35->2), ass. (0->18)
t17 = pkin(1) * qJD(1);
t16 = pkin(1) * qJD(2);
t11 = sin(qJ(2));
t12 = cos(qJ(2));
t5 = t11 * pkin(1) + qJ(3);
t1 = (t5 * t12 + (-t12 * pkin(1) - pkin(2)) * t11) * pkin(1);
t15 = t1 * qJD(1);
t8 = t12 * t16;
t14 = t8 + qJD(3);
t9 = qJD(1) + qJD(2);
t13 = t11 * t16;
t10 = qJ(3) * qJD(3);
t7 = t12 * t17;
t6 = t11 * t17;
t4 = t9 * qJ(3);
t3 = t5 * qJD(3);
t2 = -t6 - t13;
t18 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, -t8, 0, 0, 0, 0, 0, 0, 0, 0, -t13, 0, t14, t1 * qJD(2) + t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, -t7 - t8, 0, 0, 0, 0, 0, 0, 0, 0, t2, 0, t7 + t14, t15 + t3 + (-pkin(2) * t11 + qJ(3) * t12) * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, t9 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, t7, 0, 0, 0, 0, 0, 0, 0, 0, t6, 0, -t7 + qJD(3), t10 - t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3), t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, -qJ(3) * qJD(2) - t5 * qJD(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t18;
