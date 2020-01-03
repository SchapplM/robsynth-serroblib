% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S4PRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x18]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4PRRR4_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR4_inertiaDJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR4_inertiaDJ_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR4_inertiaDJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:32:39
% EndTime: 2019-12-31 16:32:40
% DurationCPUTime: 0.12s
% Computational Cost: add. (84->29), mult. (239->57), div. (0->0), fcn. (192->4), ass. (0->27)
t28 = qJD(3) + qJD(4);
t27 = pkin(5) + pkin(6);
t14 = sin(qJ(4));
t15 = sin(qJ(3));
t26 = t14 * t15;
t25 = qJD(3) * t15;
t17 = cos(qJ(3));
t24 = qJD(3) * t17;
t16 = cos(qJ(4));
t23 = qJD(4) * t16;
t22 = -0.2e1 * pkin(2) * qJD(3);
t21 = pkin(3) * t25;
t20 = qJD(4) * t14 * pkin(3);
t19 = pkin(3) * t23;
t18 = qJD(3) * t27;
t6 = t14 * t17 + t16 * t15;
t13 = -t17 * pkin(3) - pkin(2);
t10 = t27 * t17;
t9 = t27 * t15;
t8 = t17 * t18;
t7 = t15 * t18;
t5 = -t16 * t17 + t26;
t4 = t28 * t6;
t3 = -t16 * t24 - t17 * t23 + t28 * t26;
t2 = t14 * t7 - t16 * t8 + (-t10 * t16 + t14 * t9) * qJD(4);
t1 = t14 * t8 + t16 * t7 + (t10 * t14 + t16 * t9) * qJD(4);
t11 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0.2e1 * t15 * t24, 0.2e1 * (-t15 ^ 2 + t17 ^ 2) * qJD(3), 0, 0, 0, t15 * t22, t17 * t22, -0.2e1 * t6 * t3, 0.2e1 * t3 * t5 - 0.2e1 * t6 * t4, 0, 0, 0, 0.2e1 * t13 * t4 + 0.2e1 * t5 * t21, -0.2e1 * t13 * t3 + 0.2e1 * t6 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t25, -t24, 0, 0, 0, 0, 0, -t4, t3; 0, 0, 0, 0, 0, 0, t24, -t25, 0, -pkin(5) * t24, pkin(5) * t25, 0, 0, -t3, -t4, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t20, -0.2e1 * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, -t4, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20, -t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t11;
