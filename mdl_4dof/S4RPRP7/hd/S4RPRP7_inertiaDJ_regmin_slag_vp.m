% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S4RPRP7
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
% Datum: 2019-12-31 16:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4RPRP7_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP7_inertiaDJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP7_inertiaDJ_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP7_inertiaDJ_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:47:17
% EndTime: 2019-12-31 16:47:17
% DurationCPUTime: 0.09s
% Computational Cost: add. (43->22), mult. (94->43), div. (0->0), fcn. (49->2), ass. (0->16)
t15 = 2 * qJD(2);
t14 = 2 * qJD(4);
t4 = sin(qJ(3));
t13 = qJD(3) * t4;
t5 = cos(qJ(3));
t12 = qJD(3) * t5;
t6 = -pkin(1) - pkin(5);
t11 = qJD(3) * t6;
t10 = qJ(2) * qJD(3);
t9 = t4 * t11;
t8 = t5 * t11;
t7 = -t4 * pkin(3) + t5 * qJ(4);
t2 = t7 * qJD(3) + t4 * qJD(4);
t3 = qJ(2) - t7;
t1 = -t5 * qJD(4) + qJD(2) + (pkin(3) * t5 + qJ(4) * t4) * qJD(3);
t16 = [0, 0, 0, 0, t15, qJ(2) * t15, -0.2e1 * t4 * t12, 0.2e1 * (t4 ^ 2 - t5 ^ 2) * qJD(3), 0, 0, 0, 0.2e1 * qJD(2) * t4 + 0.2e1 * t5 * t10, 0.2e1 * qJD(2) * t5 - 0.2e1 * t4 * t10, 0.2e1 * t1 * t4 + 0.2e1 * t3 * t12, 0, -0.2e1 * t1 * t5 + 0.2e1 * t3 * t13, 0.2e1 * t3 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, -t13, -t12, 0, -t9, -t8, -t9, -t2, t8, t2 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, -t12, -t13, 0, t12, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, qJ(4) * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, 0, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t16;
