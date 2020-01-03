% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S4RPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x(4*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4RPPR6_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR6_inertiaDJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR6_inertiaDJ_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR6_inertiaDJ_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:40:44
% EndTime: 2019-12-31 16:40:45
% DurationCPUTime: 0.18s
% Computational Cost: add. (116->37), mult. (306->79), div. (0->0), fcn. (251->4), ass. (0->27)
t22 = sin(pkin(6));
t20 = t22 ^ 2;
t23 = cos(pkin(6));
t35 = (t23 ^ 2 + t20) * qJD(2);
t25 = cos(qJ(4));
t34 = t23 * t25;
t33 = -pkin(5) + qJ(2);
t32 = qJ(2) * t35;
t31 = t22 * qJD(2);
t30 = t22 * qJD(3);
t24 = sin(qJ(4));
t29 = t24 * qJD(4);
t28 = t25 * qJD(4);
t26 = t22 * qJ(3) + pkin(1);
t12 = t33 * t22;
t13 = t33 * t23;
t4 = t24 * t12 + t25 * t13;
t11 = t22 * t25 - t23 * t24;
t10 = t22 * t24 + t34;
t9 = 0.2e1 * t35;
t7 = (pkin(2) + pkin(3)) * t23 + t26;
t6 = t22 * t28 - t23 * t29;
t5 = qJD(4) * t10;
t3 = t25 * t12 - t24 * t13;
t2 = t11 * qJD(2) - qJD(4) * t4;
t1 = -qJD(2) * t34 - t12 * t28 + t13 * t29 - t24 * t31;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, 0.2e1 * t32, 0, 0, 0, 0, 0, 0, 0.2e1 * t23 * t30, t9, 0.2e1 * t20 * qJD(3), -0.2e1 * (-t23 * pkin(2) - t26) * t30 + 0.2e1 * t32, -0.2e1 * t11 * t5, 0.2e1 * t5 * t10 - 0.2e1 * t11 * t6, 0, 0.2e1 * t10 * t6, 0, 0, 0.2e1 * t10 * t30 + 0.2e1 * t7 * t6, 0.2e1 * t11 * t30 - 0.2e1 * t7 * t5, 0.2e1 * t1 * t10 - 0.2e1 * t2 * t11 + 0.2e1 * t3 * t5 - 0.2e1 * t4 * t6, -0.2e1 * t4 * t1 + 0.2e1 * t3 * t2 + 0.2e1 * t7 * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30, 0, 0, 0, 0, 0, 0, -t6, t5, 0, -t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, 0, 0, 0, 0, 0, 0, 0, 0, -t24 * t6 + t25 * t5 + (-t10 * t25 + t11 * t24) * qJD(4), -t1 * t24 + t2 * t25 + (-t24 * t3 + t25 * t4) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, -t6, 0, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29, -t28, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t8;
