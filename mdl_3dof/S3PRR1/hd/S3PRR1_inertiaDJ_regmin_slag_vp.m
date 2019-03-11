% Calculate minimal parameter regressor of joint inertia matrix time derivative for
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
% MMD_reg [((3+1)*3/2)x7]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S3PRR1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PRR1_inertiaDJ_regmin_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3PRR1_inertiaDJ_regmin_slag_vp: qJD has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3PRR1_inertiaDJ_regmin_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:04:00
% EndTime: 2019-03-08 18:04:01
% DurationCPUTime: 0.06s
% Computational Cost: add. (16->7), mult. (44->13), div. (0->0), fcn. (38->4), ass. (0->11)
t12 = qJD(2) + qJD(3);
t11 = pkin(2) * qJD(3);
t3 = sin(qJ(3));
t10 = t3 * t11;
t5 = cos(qJ(3));
t9 = t5 * t11;
t6 = cos(qJ(2));
t4 = sin(qJ(2));
t2 = t12 * (-t3 * t6 - t4 * t5);
t1 = t12 * (t3 * t4 - t5 * t6);
t7 = [0, 0, 0, 0, 0, 0, 0; 0, 0, -t4 * qJD(2), -t6 * qJD(2), 0, t2, t1; 0, 0, 0, 0, 0, -0.2e1 * t10, -0.2e1 * t9; 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, -t10, -t9; 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t7;
