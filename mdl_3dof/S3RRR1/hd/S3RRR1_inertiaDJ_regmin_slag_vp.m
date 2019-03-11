% Calculate minimal parameter regressor of joint inertia matrix time derivative for
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
% MMD_reg [((3+1)*3/2)x9]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S3RRR1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RRR1_inertiaDJ_regmin_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RRR1_inertiaDJ_regmin_slag_vp: qJD has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S3RRR1_inertiaDJ_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:08:08
% EndTime: 2019-03-08 18:08:08
% DurationCPUTime: 0.08s
% Computational Cost: add. (30->15), mult. (104->28), div. (0->0), fcn. (58->4), ass. (0->18)
t18 = pkin(1) * qJD(2);
t7 = sin(qJ(2));
t15 = t7 * t18;
t6 = sin(qJ(3));
t17 = qJD(3) * t6;
t19 = t7 * pkin(1) * t17 + t6 * t15;
t8 = cos(qJ(3));
t16 = qJD(3) * t8;
t9 = cos(qJ(2));
t14 = t9 * t18;
t13 = pkin(2) * t17;
t12 = pkin(2) * t16;
t5 = t9 * pkin(1) + pkin(2);
t11 = (-pkin(2) - t5) * qJD(3);
t10 = (-t7 * t16 + (-t6 * t9 - t7 * t8) * qJD(2)) * pkin(1);
t2 = -t5 * t17 + t10;
t1 = (-qJD(3) * t5 - t14) * t8 + t19;
t3 = [0, 0, 0, 0, -0.2e1 * t15, -0.2e1 * t14, 0, 0.2e1 * t2, 0.2e1 * t1; 0, 0, 0, 0, -t15, -t14, 0, t6 * t11 + t10 (t11 - t14) * t8 + t19; 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t13, -0.2e1 * t12; 0, 0, 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, -t13, -t12; 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t3;
