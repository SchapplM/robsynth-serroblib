% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RPPRP6
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
% MMD_reg [((5+1)*5/2)x21]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPPRP6_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP6_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP6_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP6_inertiaDJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:55:16
% EndTime: 2019-12-31 17:55:17
% DurationCPUTime: 0.20s
% Computational Cost: add. (248->45), mult. (530->86), div. (0->0), fcn. (465->4), ass. (0->33)
t23 = sin(pkin(7));
t24 = cos(pkin(7));
t13 = (t23 ^ 2 + t24 ^ 2) * qJD(3);
t42 = 2 * qJD(2);
t41 = 2 * qJD(5);
t37 = sin(qJ(4));
t33 = t37 * t23;
t38 = cos(qJ(4));
t34 = t38 * t24;
t10 = t33 - t34;
t31 = qJD(4) * t37;
t32 = qJD(4) * t38;
t8 = -t23 * t32 - t24 * t31;
t40 = t10 * t8;
t25 = -pkin(1) - qJ(3);
t39 = -pkin(6) + t25;
t18 = t23 * pkin(3) + qJ(2);
t35 = qJ(2) * qJD(2);
t30 = t39 * t37;
t11 = t38 * t23 + t37 * t24;
t9 = -t23 * t31 + t24 * t32;
t29 = -t11 * t9 + t40;
t28 = t39 * t34;
t27 = pkin(4) * t8 + t9 * qJ(5) + t11 * qJD(5);
t12 = t39 * t23;
t2 = t11 * qJD(3) - qJD(4) * t28 + t12 * t31;
t3 = t12 * t32 - qJD(3) * t33 + (t38 * qJD(3) + qJD(4) * t30) * t24;
t5 = t37 * t12 - t28;
t6 = t38 * t12 + t24 * t30;
t26 = t3 * t10 - t2 * t11 - t5 * t8 + t6 * t9;
t4 = t11 * pkin(4) + t10 * qJ(5) + t18;
t1 = t9 * pkin(4) - t8 * qJ(5) + t10 * qJD(5) + qJD(2);
t7 = [0, 0, 0, 0, t42, 0.2e1 * t35, t23 * t42, t24 * t42, 0.2e1 * t13, -0.2e1 * t25 * t13 + 0.2e1 * t35, -0.2e1 * t40, 0.2e1 * t10 * t9 - 0.2e1 * t8 * t11, 0, 0, 0, 0.2e1 * qJD(2) * t11 + 0.2e1 * t18 * t9, -0.2e1 * qJD(2) * t10 + 0.2e1 * t18 * t8, 0.2e1 * t1 * t11 + 0.2e1 * t4 * t9, -0.2e1 * t26, 0.2e1 * t1 * t10 - 0.2e1 * t4 * t8, 0.2e1 * t4 * t1 - 0.2e1 * t6 * t2 + 0.2e1 * t5 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t29, 0, t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2), 0, 0, 0, 0, 0, t9, t8, t9, 0, -t8, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t9, 0, -t3, t2, -t3, -t27, -t2, -t3 * pkin(4) - t2 * qJ(5) + t6 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t9, t8, 0, t9, t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, qJ(5) * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t7;
