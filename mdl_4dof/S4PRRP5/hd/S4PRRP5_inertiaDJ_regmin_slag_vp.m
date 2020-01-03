% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S4PRRP5
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
% MMD_reg [((4+1)*4/2)x13]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:29
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4PRRP5_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP5_inertiaDJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP5_inertiaDJ_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP5_inertiaDJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:29:22
% EndTime: 2019-12-31 16:29:22
% DurationCPUTime: 0.12s
% Computational Cost: add. (50->25), mult. (162->57), div. (0->0), fcn. (106->4), ass. (0->23)
t8 = sin(qJ(3));
t6 = t8 ^ 2;
t10 = cos(qJ(3));
t7 = t10 ^ 2;
t22 = t6 + t7;
t21 = -qJ(4) - pkin(5);
t20 = t8 * qJD(3);
t9 = sin(qJ(2));
t19 = t9 * qJD(2);
t18 = t10 * qJD(3);
t11 = cos(qJ(2));
t17 = t11 * qJD(2);
t16 = -0.2e1 * pkin(2) * qJD(3);
t15 = pkin(3) * t20;
t14 = qJD(3) * t21;
t13 = -t8 * t17 - t9 * t18;
t1 = t10 * qJD(4) + t8 * t14;
t2 = -t8 * qJD(4) + t10 * t14;
t3 = t21 * t8;
t4 = t21 * t10;
t12 = t1 * t10 - t2 * t8 + (-t10 * t3 + t4 * t8) * qJD(3);
t5 = -t10 * pkin(3) - pkin(2);
t23 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * (-0.1e1 + t22) * t9 * t17; 0, 0, -t19, -t17, 0, 0, 0, 0, 0, -t10 * t19 - t11 * t20, -t11 * t18 + t8 * t19, t22 * t17, (-t15 + (-t10 * t4 - t3 * t8) * qJD(2)) * t11 + (qJD(2) * t5 + t12) * t9; 0, 0, 0, 0, 0.2e1 * t8 * t18, 0.2e1 * (-t6 + t7) * qJD(3), 0, 0, 0, t8 * t16, t10 * t16, 0.2e1 * t12, -0.2e1 * t4 * t1 + 0.2e1 * t5 * t15 + 0.2e1 * t3 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, -t10 * t17 + t9 * t20, 0, t13 * pkin(3); 0, 0, 0, 0, 0, 0, t18, -t20, 0, -pkin(5) * t18, pkin(5) * t20, -pkin(3) * t18, t2 * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t23;
