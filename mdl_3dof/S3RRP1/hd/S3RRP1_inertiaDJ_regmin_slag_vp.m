% Calculate minimal parameter regressor of joint inertia matrix time derivative for
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
% MMD_reg [((3+1)*3/2)x9]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S3RRP1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RRP1_inertiaDJ_regmin_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RRP1_inertiaDJ_regmin_slag_vp: qJD has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3RRP1_inertiaDJ_regmin_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:06:52
% EndTime: 2019-03-08 18:06:52
% DurationCPUTime: 0.06s
% Computational Cost: add. (13->10), mult. (45->17), div. (0->0), fcn. (16->2), ass. (0->10)
t7 = 2 * qJD(3);
t9 = pkin(1) * qJD(2);
t5 = sin(qJ(2));
t8 = t5 * t9;
t6 = cos(qJ(2));
t4 = t6 * t9;
t3 = t5 * pkin(1) + qJ(3);
t2 = -0.2e1 * t8;
t1 = t4 + qJD(3);
t10 = [0, 0, 0, 0, t2, -0.2e1 * t4, t2, 0.2e1 * t1, 0.2e1 * t3 * t1 + 0.2e1 * (-t6 * pkin(1) - pkin(2)) * t8; 0, 0, 0, 0, -t8, -t4, -t8, t7 + t4, -pkin(2) * t8 + t1 * qJ(3) + t3 * qJD(3); 0, 0, 0, 0, 0, 0, 0, t7, qJ(3) * t7; 0, 0, 0, 0, 0, 0, 0, 0, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t10;
