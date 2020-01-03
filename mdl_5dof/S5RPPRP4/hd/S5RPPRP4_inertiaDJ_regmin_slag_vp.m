% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RPPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x18]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPPRP4_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP4_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP4_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP4_inertiaDJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:52:20
% EndTime: 2019-12-31 17:52:20
% DurationCPUTime: 0.16s
% Computational Cost: add. (133->32), mult. (232->68), div. (0->0), fcn. (169->4), ass. (0->26)
t13 = sin(qJ(4));
t14 = cos(qJ(4));
t12 = cos(pkin(7));
t25 = t12 * qJD(2);
t18 = -qJD(5) + t25;
t11 = sin(pkin(7));
t15 = -pkin(1) - pkin(2);
t28 = t12 * qJ(2) + t11 * t15;
t6 = -pkin(6) + t28;
t27 = qJ(5) - t6;
t19 = qJD(4) * t27;
t1 = t13 * t19 + t18 * t14;
t2 = -t18 * t13 + t14 * t19;
t3 = t27 * t13;
t4 = t27 * t14;
t30 = (-t13 * t4 + t14 * t3) * qJD(4) - t1 * t14 + t2 * t13;
t29 = 0.2e1 * qJD(2);
t26 = t11 * qJD(2);
t24 = t13 * qJD(4);
t23 = t14 * qJD(4);
t22 = pkin(4) * t24;
t21 = t11 * t23;
t20 = -t11 * qJ(2) + t12 * t15;
t5 = pkin(3) - t20;
t7 = -t22 + t26;
t8 = [0, 0, 0, 0, t29, qJ(2) * t29, 0.2e1 * t26, 0.2e1 * t25, (-t20 * t11 + t28 * t12) * t29, 0.2e1 * t13 * t23, 0.2e1 * (-t13 ^ 2 + t14 ^ 2) * qJD(4), 0, 0, 0, 0.2e1 * t14 * t26 - 0.2e1 * t5 * t24, -0.2e1 * t13 * t26 - 0.2e1 * t5 * t23, 0.2e1 * t30, -0.2e1 * t4 * t1 + 0.2e1 * t3 * t2 + 0.2e1 * (t14 * pkin(4) + t5) * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12 * t24, t12 * t23, 0, -t30 * t11 - t7 * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1 * t13 + t2 * t14 + (-t13 * t3 - t14 * t4) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, t24, 0, -t13 * t25 - t6 * t23, -t14 * t25 + t6 * t24, pkin(4) * t23, t2 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, t11 * t24, 0, -pkin(4) * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, -t23, 0, -t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t8;
