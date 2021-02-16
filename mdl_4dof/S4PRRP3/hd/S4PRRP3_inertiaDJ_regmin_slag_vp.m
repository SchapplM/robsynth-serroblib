% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S4PRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x15]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-14 22:27
% Revision: beb2ba9bd8c5bd556f42a244985f3dab86917626 (2021-01-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4PRRP3_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP3_inertiaDJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP3_inertiaDJ_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP3_inertiaDJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-14 22:27:05
% EndTime: 2021-01-14 22:27:05
% DurationCPUTime: 0.10s
% Computational Cost: add. (46->23), mult. (132->46), div. (0->0), fcn. (75->2), ass. (0->17)
t17 = 2 * qJD(3);
t10 = cos(qJ(3));
t16 = t10 * pkin(3);
t15 = -qJ(4) - pkin(5);
t9 = sin(qJ(3));
t14 = qJD(3) * t9;
t7 = qJD(3) * t10;
t13 = -2 * pkin(2) * qJD(3);
t12 = 0.2e1 * t14;
t11 = pkin(3) * t14;
t8 = t9 ^ 2;
t5 = -pkin(2) - t16;
t4 = t15 * t10;
t3 = t15 * t9;
t2 = -t9 * qJD(4) + t15 * t7;
t1 = -t10 * qJD(4) - t15 * t14;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9 * t1 + t10 * t2 + (-t10 * t4 - t3 * t9) * qJD(3); 0, 0, 0, 0, t10 * t12, (t10 ^ 2 - t8) * t17, 0, 0, 0, t9 * t13, t10 * t13, (t5 - t16) * t12, (pkin(3) * t8 + t10 * t5) * t17, -0.2e1 * t1 * t10 - 0.2e1 * t2 * t9 + 0.2e1 * (-t10 * t3 + t4 * t9) * qJD(3), 0.2e1 * t4 * t1 + 0.2e1 * t5 * t11 + 0.2e1 * t3 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, -t7, -t14, -t7, 0, -t11; 0, 0, 0, 0, 0, 0, t7, -t14, 0, -pkin(5) * t7, pkin(5) * t14, t2, t1, -pkin(3) * t7, t2 * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, t7, 0, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t6;
