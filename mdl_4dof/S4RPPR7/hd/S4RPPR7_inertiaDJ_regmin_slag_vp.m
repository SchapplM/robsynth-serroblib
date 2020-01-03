% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S4RPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta3]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x17]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4RPPR7_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR7_inertiaDJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR7_inertiaDJ_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR7_inertiaDJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:41:41
% EndTime: 2019-12-31 16:41:41
% DurationCPUTime: 0.10s
% Computational Cost: add. (41->20), mult. (121->43), div. (0->0), fcn. (94->4), ass. (0->17)
t12 = sin(pkin(6));
t13 = cos(pkin(6));
t15 = sin(qJ(4));
t16 = cos(qJ(4));
t24 = -t15 * t12 + t16 * t13;
t7 = (t12 ^ 2 + t13 ^ 2) * qJD(3);
t23 = 2 * qJD(2);
t14 = -pkin(1) - qJ(3);
t22 = -pkin(5) + t14;
t18 = qJ(2) * qJD(2);
t3 = t16 * t12 + t15 * t13;
t9 = t12 * pkin(3) + qJ(2);
t6 = t22 * t13;
t5 = t22 * t12;
t2 = t24 * qJD(4);
t1 = t3 * qJD(4);
t4 = [0, 0, 0, 0, t23, 2 * t18, t12 * t23, t13 * t23, 0.2e1 * t7, -0.2e1 * t14 * t7 + (2 * t18), -0.2e1 * t24 * t1, 0.2e1 * t1 * t3 - 0.2e1 * t2 * t24, 0, 0, 0, 0.2e1 * qJD(2) * t3 + 0.2e1 * t9 * t2, 0.2e1 * qJD(2) * t24 - 0.2e1 * t9 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2), 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, -t2, 0, (-t15 * t6 - t16 * t5) * qJD(4) - t24 * qJD(3), (t15 * t5 - t16 * t6) * qJD(4) + t3 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t4;
