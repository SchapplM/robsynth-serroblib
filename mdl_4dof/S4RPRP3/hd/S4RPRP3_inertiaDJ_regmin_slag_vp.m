% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S4RPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x15]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 10:20
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4RPRP3_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP3_inertiaDJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP3_inertiaDJ_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP3_inertiaDJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 10:20:36
% EndTime: 2021-01-15 10:20:37
% DurationCPUTime: 0.10s
% Computational Cost: add. (68->25), mult. (154->48), div. (0->0), fcn. (97->4), ass. (0->19)
t20 = 2 * qJD(3);
t13 = cos(qJ(3));
t19 = t13 * pkin(3);
t6 = sin(pkin(6)) * pkin(1) + pkin(5);
t18 = qJ(4) + t6;
t12 = sin(qJ(3));
t17 = qJD(3) * t12;
t9 = qJD(3) * t13;
t16 = 0.2e1 * t17;
t15 = 0.2e1 * t9;
t14 = pkin(3) * t17;
t7 = -cos(pkin(6)) * pkin(1) - pkin(2);
t10 = t12 ^ 2;
t5 = t7 - t19;
t4 = t18 * t13;
t3 = t18 * t12;
t2 = -t12 * qJD(4) - t18 * t9;
t1 = -t13 * qJD(4) + t18 * t17;
t8 = [0, 0, 0, 0, t12 * t15, (t13 ^ 2 - t10) * t20, 0, 0, 0, t7 * t16, t7 * t15, (t5 - t19) * t16, (pkin(3) * t10 + t13 * t5) * t20, -0.2e1 * t1 * t13 - 0.2e1 * t2 * t12 + 0.2e1 * (-t12 * t4 + t13 * t3) * qJD(3), -0.2e1 * t4 * t1 + 0.2e1 * t5 * t14 - 0.2e1 * t3 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1 * t12 + t2 * t13 + (t12 * t3 + t13 * t4) * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t9, -t17, 0, -t6 * t9, t6 * t17, t2, t1, -pkin(3) * t9, t2 * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, -t9, -t17, -t9, 0, -t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, t9, 0, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t8;
