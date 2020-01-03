% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S4RPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x16]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4RPRR5_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR5_inertiaDJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR5_inertiaDJ_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR5_inertiaDJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:51:40
% EndTime: 2019-12-31 16:51:40
% DurationCPUTime: 0.12s
% Computational Cost: add. (70->28), mult. (166->53), div. (0->0), fcn. (104->4), ass. (0->26)
t27 = 2 * qJD(2);
t12 = sin(qJ(4));
t13 = sin(qJ(3));
t15 = cos(qJ(3));
t16 = -pkin(1) - pkin(2);
t17 = t15 * qJ(2) + t13 * t16;
t2 = t13 * qJD(2) + t17 * qJD(3);
t26 = t2 * t12;
t14 = cos(qJ(4));
t25 = t2 * t14;
t24 = qJD(3) * t13;
t23 = qJD(3) * t15;
t22 = t12 * qJD(4);
t21 = t14 * qJD(4);
t20 = -0.2e1 * pkin(3) * qJD(4);
t19 = t12 * t21;
t6 = t13 * qJ(2) - t15 * t16 + pkin(3);
t18 = qJD(4) * (pkin(3) + t6);
t9 = 0.2e1 * t19;
t8 = (-t12 ^ 2 + t14 ^ 2) * qJD(4);
t7 = -pkin(6) + t17;
t5 = 0.2e1 * t8;
t4 = t14 * t24 + t15 * t22;
t3 = t12 * t24 - t15 * t21;
t1 = qJ(2) * t24 - t15 * qJD(2) - t16 * t23;
t10 = [0, 0, 0, 0, t27, qJ(2) * t27, 0, 0.2e1 * t2, -0.2e1 * t1, t9, t5, 0, 0, 0, -0.2e1 * t6 * t22 + 0.2e1 * t25, -0.2e1 * t6 * t21 - 0.2e1 * t26; 0, 0, 0, 0, 0, 0, 0, t24, t23, 0, 0, 0, 0, 0, t4, -t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, -t2, t1, -0.2e1 * t19, -0.2e1 * t8, 0, 0, 0, t12 * t18 - t25, t14 * t18 + t26; 0, 0, 0, 0, 0, 0, 0, -t24, -t23, 0, 0, 0, 0, 0, -t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, t5, 0, 0, 0, t12 * t20, t14 * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, t22, 0, t12 * t1 - t7 * t21, t14 * t1 + t7 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12 * t23 - t13 * t21, t13 * t22 - t14 * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, -t22, 0, -pkin(6) * t21, pkin(6) * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t10;
