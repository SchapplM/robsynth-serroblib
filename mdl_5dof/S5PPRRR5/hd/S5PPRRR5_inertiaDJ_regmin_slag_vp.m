% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5PPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x15]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PPRRR5_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR5_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR5_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRR5_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:35:45
% EndTime: 2019-12-31 17:35:46
% DurationCPUTime: 0.15s
% Computational Cost: add. (66->27), mult. (209->48), div. (0->0), fcn. (170->6), ass. (0->30)
t34 = qJD(3) + qJD(4);
t18 = sin(qJ(4));
t33 = t18 * pkin(3);
t21 = cos(qJ(4));
t15 = -t21 * pkin(3) - pkin(4);
t20 = cos(qJ(5));
t16 = qJD(5) * t20;
t17 = sin(qJ(5));
t25 = qJD(4) * t33;
t32 = t15 * t16 + t17 * t25;
t19 = sin(qJ(3));
t31 = t18 * t19;
t22 = cos(qJ(3));
t30 = qJD(3) * t22;
t29 = qJD(4) * t21;
t28 = qJD(5) * t17;
t27 = pkin(4) * t28;
t26 = pkin(4) * t16;
t24 = pkin(3) * t29;
t7 = t18 * t22 + t21 * t19;
t23 = t15 * t28 - t20 * t25;
t14 = pkin(7) + t33;
t11 = 0.2e1 * t17 * t16;
t6 = -t21 * t22 + t31;
t5 = 0.2e1 * (-t17 ^ 2 + t20 ^ 2) * qJD(5);
t4 = t34 * t7;
t3 = -t21 * t30 - t22 * t29 + t34 * t31;
t2 = -t4 * t20 + t6 * t28;
t1 = t6 * t16 + t4 * t17;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, -qJD(3) * t19, -t30, 0, -t4, t3, 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, 0, -0.2e1 * t25, -0.2e1 * t24, t11, t5, 0, 0, 0, 0.2e1 * t23, 0.2e1 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, -t4, t3, 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, 0, -t25, -t24, t11, t5, 0, 0, 0, t23 - t27, -t26 + t32; 0, 0, 0, 0, 0, 0, 0, 0, t11, t5, 0, 0, 0, -0.2e1 * t27, -0.2e1 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7 * t16 + t17 * t3, t20 * t3 + t7 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, -t28, 0, -t14 * t16 - t17 * t24, t14 * t28 - t20 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, -t28, 0, -pkin(7) * t16, pkin(7) * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t8;
