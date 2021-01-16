% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S4RPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x17]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 10:27
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4RPRP6_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP6_inertiaDJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP6_inertiaDJ_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP6_inertiaDJ_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 10:27:40
% EndTime: 2021-01-15 10:27:41
% DurationCPUTime: 0.10s
% Computational Cost: add. (71->27), mult. (139->52), div. (0->0), fcn. (78->2), ass. (0->18)
t17 = 2 * qJD(2);
t8 = sin(qJ(3));
t16 = qJD(3) * t8;
t9 = cos(qJ(3));
t15 = qJD(3) * t9;
t10 = -pkin(1) - pkin(5);
t14 = qJ(4) - t10;
t13 = qJD(3) * t10;
t12 = qJ(2) * qJD(3);
t11 = pkin(3) * t16;
t5 = t14 * t9;
t7 = t8 * pkin(3) + qJ(2);
t6 = pkin(3) * t15 + qJD(2);
t4 = t14 * t8;
t3 = -qJD(3) * t5 - t8 * qJD(4);
t2 = -t9 * qJD(4) + t14 * t16;
t1 = t2 * t9 + t3 * t8 + (-t4 * t9 + t5 * t8) * qJD(3);
t18 = [0, 0, 0, 0, t17, qJ(2) * t17, -0.2e1 * t8 * t15, 0.2e1 * (t8 ^ 2 - t9 ^ 2) * qJD(3), 0, 0, 0, 0.2e1 * qJD(2) * t8 + 0.2e1 * t9 * t12, 0.2e1 * qJD(2) * t9 - 0.2e1 * t8 * t12, 0.2e1 * t7 * t15 + 0.2e1 * t6 * t8, -0.2e1 * t7 * t16 + 0.2e1 * t6 * t9, -0.2e1 * t1, -0.2e1 * t5 * t2 - 0.2e1 * t4 * t3 + 0.2e1 * t7 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, -t16, -t15, 0, -t8 * t13, -t9 * t13, t2, -t3, t11, t2 * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, -t15, -t16, -t15, 0, -t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, -t16, 0, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t18;
