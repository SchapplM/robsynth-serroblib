% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5PRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x19]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRRPP4_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP4_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP4_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRRPP4_inertiaDJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:41:12
% EndTime: 2019-12-31 17:41:13
% DurationCPUTime: 0.18s
% Computational Cost: add. (100->42), mult. (247->76), div. (0->0), fcn. (143->2), ass. (0->28)
t14 = cos(qJ(3));
t13 = sin(qJ(3));
t25 = t13 * qJ(4);
t28 = pkin(3) + pkin(4);
t29 = t28 * t14 + t25;
t16 = 2 * qJD(4);
t27 = pkin(6) - qJ(5);
t26 = qJ(4) * t14;
t10 = t13 * qJD(3);
t24 = t13 * qJD(4);
t11 = t14 * qJD(3);
t23 = t14 * qJD(4);
t22 = -0.2e1 * pkin(2) * qJD(3);
t21 = pkin(6) * t10;
t20 = pkin(6) * t11;
t8 = t27 * t14;
t19 = -pkin(3) * t10 + t24;
t18 = -t14 * pkin(3) - t25;
t17 = t18 * qJD(3) + t23;
t12 = qJ(4) * t16;
t7 = t27 * t13;
t6 = -pkin(2) + t18;
t5 = pkin(2) + t29;
t4 = qJ(4) * t11 + t19;
t3 = qJD(3) * t8 - t13 * qJD(5);
t2 = -t14 * qJD(5) - t27 * t10;
t1 = (-pkin(4) * t13 + t26) * qJD(3) + t19;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13 * t2 - t14 * t3 + (t13 * t7 + t14 * t8) * qJD(3); 0, 0, 0, 0, 0.2e1 * t13 * t11, 0.2e1 * (-t13 ^ 2 + t14 ^ 2) * qJD(3), 0, 0, 0, t13 * t22, t14 * t22, 0.2e1 * t6 * t10 + 0.2e1 * t4 * t14, 0, -0.2e1 * t6 * t11 + 0.2e1 * t4 * t13, -0.2e1 * t6 * t4, 0.2e1 * t1 * t14 - 0.2e1 * t5 * t10, 0.2e1 * t1 * t13 + 0.2e1 * t5 * t11, -0.2e1 * t3 * t13 - 0.2e1 * t2 * t14 + 0.2e1 * (t13 * t8 - t14 * t7) * qJD(3), 0.2e1 * t5 * t1 + 0.2e1 * t8 * t2 + 0.2e1 * t7 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, -t11, -t10, 0, t11, t4, -t10, t11, 0, t24 + (-t13 * t28 + t26) * qJD(3); 0, 0, 0, 0, 0, 0, t11, -t10, 0, -t20, t21, -t20, t17, -t21, t17 * pkin(6), -t3, t2, t29 * qJD(3) - t23, t2 * qJ(4) + t8 * qJD(4) - t28 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, t12, 0, t16, 0, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, 0, 0, 0, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, 0, t20, 0, 0, -t11, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, t11, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t9;
