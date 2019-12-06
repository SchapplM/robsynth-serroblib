% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5PRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x20]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:44
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRRRP3_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP3_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP3_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP3_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:44:13
% EndTime: 2019-12-05 16:44:14
% DurationCPUTime: 0.27s
% Computational Cost: add. (341->59), mult. (835->110), div. (0->0), fcn. (712->4), ass. (0->39)
t45 = qJD(3) + qJD(4);
t25 = sin(qJ(3));
t44 = pkin(6) + pkin(7);
t32 = qJD(3) * t44;
t15 = t25 * t32;
t26 = cos(qJ(3));
t16 = t26 * t32;
t24 = sin(qJ(4));
t17 = t44 * t25;
t18 = t44 * t26;
t42 = cos(qJ(4));
t28 = -t42 * t17 - t24 * t18;
t3 = -t28 * qJD(4) + t42 * t15 + t24 * t16;
t14 = t24 * t26 + t42 * t25;
t11 = t45 * t14;
t43 = t11 * pkin(4);
t31 = t42 * qJD(4);
t33 = t42 * t26;
t40 = t24 * t25;
t10 = -qJD(3) * t33 - t26 * t31 + t45 * t40;
t41 = t14 * t10;
t39 = t25 * qJD(3);
t38 = t26 * qJD(3);
t37 = -0.2e1 * pkin(2) * qJD(3);
t36 = pkin(3) * t39;
t35 = qJD(4) * t24 * pkin(3);
t34 = t42 * pkin(3);
t23 = -t26 * pkin(3) - pkin(2);
t30 = pkin(3) * t31;
t27 = t24 * t17 - t42 * t18;
t4 = t27 * qJD(4) + t24 * t15 - t42 * t16;
t22 = t34 + pkin(4);
t13 = -t33 + t40;
t9 = t36 + t43;
t8 = -t13 * qJ(5) - t27;
t7 = -t14 * qJ(5) + t28;
t2 = t10 * qJ(5) - t14 * qJD(5) + t4;
t1 = -t11 * qJ(5) - t13 * qJD(5) - t3;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t13 * t11 - 0.2e1 * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14 * t1 - t10 * t8 - t11 * t7 - t13 * t2; 0, 0, 0, 0, 0.2e1 * t25 * t38, 0.2e1 * (-t25 ^ 2 + t26 ^ 2) * qJD(3), 0, 0, 0, t25 * t37, t26 * t37, -0.2e1 * t41, 0.2e1 * t13 * t10 - 0.2e1 * t14 * t11, 0, 0, 0, 0.2e1 * t23 * t11 + 0.2e1 * t13 * t36, -0.2e1 * t23 * t10 + 0.2e1 * t14 * t36, -0.2e1 * t1 * t13 + 0.2e1 * t7 * t10 - 0.2e1 * t8 * t11 - 0.2e1 * t2 * t14, 0.2e1 * t8 * t1 + 0.2e1 * t7 * t2 + 0.2e1 * (t13 * pkin(4) + t23) * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t39, -t38, 0, 0, 0, 0, 0, -t11, t10, 0, -t11 * t22 + (-t10 * t24 + (t13 * t24 + t42 * t14) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, t38, -t39, 0, -pkin(6) * t38, pkin(6) * t39, 0, 0, -t10, -t11, 0, t4, t3, t22 * t10 + (-t11 * t24 + (-t42 * t13 + t14 * t24) * qJD(4)) * pkin(3), t2 * t22 + (t1 * t24 + (-t24 * t7 + t42 * t8) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t35, -0.2e1 * t30, 0, 0.2e1 * (t34 - t22) * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, t10, 0, -t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, -t11, 0, t4, t3, pkin(4) * t10, t2 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35, -t30, 0, -pkin(4) * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t5;
