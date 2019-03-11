% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S3RRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2,d3]';
% 
% Output:
% MMD_reg [((3+1)*3/2)x(3*10)]
%   inertial parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S3RRR1_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RRR1_inertiaDJ_reg2_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RRR1_inertiaDJ_reg2_slag_vp: qJD has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S3RRR1_inertiaDJ_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:08:08
% EndTime: 2019-03-08 18:08:08
% DurationCPUTime: 0.09s
% Computational Cost: add. (58->20), mult. (182->38), div. (0->0), fcn. (110->4), ass. (0->20)
t10 = cos(qJ(3));
t9 = sin(qJ(2));
t21 = t10 * t9;
t20 = pkin(1) * qJD(2);
t8 = sin(qJ(3));
t19 = qJD(3) * t8;
t18 = qJD(3) * t10;
t17 = t8 * t9 * pkin(1);
t16 = t9 * t20;
t15 = pkin(2) * t19;
t11 = cos(qJ(2));
t14 = t11 * t20;
t13 = pkin(2) * t18;
t7 = t11 * pkin(1) + pkin(2);
t1 = -t7 * t18 - t10 * t14 + (qJD(2) + qJD(3)) * t17;
t12 = (-t9 * t18 + (-t11 * t8 - t21) * qJD(2)) * pkin(1);
t4 = pkin(1) * t21 + t8 * t7;
t3 = t10 * t7 - t17;
t2 = -t7 * t19 + t12;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t16, -0.2e1 * t14, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t2, 0.2e1 * t1, 0, -0.2e1 * t4 * t1 + 0.2e1 * t3 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, -t14, 0, 0, 0, 0, 0, 0, 0, 0 (-pkin(2) - t7) * t19 + t12, t1 - t13, 0 (-t1 * t8 + t10 * t2 + (t10 * t4 - t3 * t8) * qJD(3)) * pkin(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t15, -0.2e1 * t13, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, -t13, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t5;
