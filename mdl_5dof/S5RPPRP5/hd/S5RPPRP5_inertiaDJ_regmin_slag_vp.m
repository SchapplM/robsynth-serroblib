% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RPPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x22]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPPRP5_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP5_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP5_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP5_inertiaDJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:53:43
% EndTime: 2019-12-31 17:53:44
% DurationCPUTime: 0.24s
% Computational Cost: add. (202->57), mult. (501->104), div. (0->0), fcn. (425->4), ass. (0->31)
t26 = sin(pkin(7));
t24 = t26 ^ 2;
t27 = cos(pkin(7));
t42 = (t27 ^ 2 + t24) * qJD(2);
t41 = 2 * qJD(5);
t40 = -pkin(6) + qJ(2);
t39 = qJ(2) * t42;
t38 = t26 * qJD(2);
t37 = t26 * qJD(3);
t36 = t27 * qJD(2);
t28 = sin(qJ(4));
t35 = t28 * qJD(4);
t29 = cos(qJ(4));
t34 = t29 * qJD(4);
t32 = -t27 * pkin(2) - t26 * qJ(3) - pkin(1);
t31 = t40 * t26;
t9 = t27 * pkin(3) - t32;
t30 = t29 * t31;
t12 = t26 * t28 + t27 * t29;
t14 = t40 * t27;
t6 = t29 * t14 + t28 * t31;
t13 = t26 * t29 - t27 * t28;
t11 = 0.2e1 * t42;
t8 = t26 * t34 - t27 * t35;
t7 = qJD(4) * t12;
t5 = t28 * t14 - t30;
t4 = t12 * pkin(4) - t13 * qJ(5) + t9;
t3 = qJD(4) * t6 + t28 * t36 - t29 * t38;
t2 = -qJD(4) * t30 + t14 * t35 - t28 * t38 - t29 * t36;
t1 = t8 * pkin(4) + t7 * qJ(5) - t13 * qJD(5) + t37;
t10 = [0, 0, 0, 0, 0, t11, 0.2e1 * t39, 0.2e1 * t27 * t37, t11, 0.2e1 * t24 * qJD(3), -0.2e1 * t32 * t37 + 0.2e1 * t39, -0.2e1 * t13 * t7, 0.2e1 * t7 * t12 - 0.2e1 * t13 * t8, 0, 0, 0, 0.2e1 * t12 * t37 + 0.2e1 * t9 * t8, 0.2e1 * t13 * t37 - 0.2e1 * t9 * t7, 0.2e1 * t1 * t12 + 0.2e1 * t4 * t8, 0.2e1 * t2 * t12 + 0.2e1 * t3 * t13 - 0.2e1 * t5 * t7 - 0.2e1 * t6 * t8, -0.2e1 * t1 * t13 + 0.2e1 * t4 * t7, 0.2e1 * t4 * t1 - 0.2e1 * t6 * t2 + 0.2e1 * t5 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, 0, 0, 0, 0, 0, -t8, t7, -t8, 0, -t7, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, 0, 0, 0, 0, 0, 0, 0, 0, -t28 * t8 + t29 * t7 + (-t12 * t29 + t13 * t28) * qJD(4), 0, -t2 * t28 - t3 * t29 + (t28 * t5 + t29 * t6) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, -t8, 0, -t3, t2, -t3, t7 * pkin(4) - t8 * qJ(5) - t12 * qJD(5), -t2, -t3 * pkin(4) - t2 * qJ(5) + t6 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35, -t34, -t35, 0, t34, t28 * qJD(5) + (-pkin(4) * t28 + qJ(5) * t29) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, qJ(5) * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t10;
