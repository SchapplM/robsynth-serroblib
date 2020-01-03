% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S4PRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x13]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:27
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4PRRP3_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP3_inertiaDJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP3_inertiaDJ_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP3_inertiaDJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:26:56
% EndTime: 2019-12-31 16:26:56
% DurationCPUTime: 0.10s
% Computational Cost: add. (36->19), mult. (102->43), div. (0->0), fcn. (58->2), ass. (0->13)
t12 = -qJ(4) - pkin(5);
t5 = sin(qJ(3));
t11 = t5 * qJD(3);
t6 = cos(qJ(3));
t10 = t6 * qJD(3);
t9 = -0.2e1 * pkin(2) * qJD(3);
t8 = pkin(3) * t11;
t7 = qJD(3) * t12;
t4 = t12 * t6;
t3 = t12 * t5;
t2 = -t5 * qJD(4) + t6 * t7;
t1 = t6 * qJD(4) + t5 * t7;
t13 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5 * t1 + t6 * t2 + (-t3 * t5 - t4 * t6) * qJD(3); 0, 0, 0, 0, 0.2e1 * t5 * t10, 0.2e1 * (-t5 ^ 2 + t6 ^ 2) * qJD(3), 0, 0, 0, t5 * t9, t6 * t9, 0.2e1 * t1 * t6 - 0.2e1 * t2 * t5 + 0.2e1 * (-t3 * t6 + t4 * t5) * qJD(3), -0.2e1 * t4 * t1 + 0.2e1 * t3 * t2 + 0.2e1 * (-t6 * pkin(3) - pkin(2)) * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, -t10, 0, -t8; 0, 0, 0, 0, 0, 0, t10, -t11, 0, -pkin(5) * t10, pkin(5) * t11, -pkin(3) * t10, t2 * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t13;
