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
% MMD_reg [((4+1)*4/2)x15]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2019-12-31 16:46:14
% EndTime: 2019-12-31 16:46:14
% DurationCPUTime: 0.11s
% Computational Cost: add. (56->21), mult. (109->43), div. (0->0), fcn. (60->2), ass. (0->16)
t15 = 2 * qJD(2);
t9 = -pkin(1) - pkin(5);
t14 = qJ(4) - t9;
t7 = sin(qJ(3));
t13 = t7 * qJD(3);
t8 = cos(qJ(3));
t12 = t8 * qJD(3);
t11 = qJ(2) * qJD(3);
t10 = pkin(3) * t13;
t5 = t14 * t8;
t6 = pkin(3) * t12 + qJD(2);
t4 = t14 * t7;
t3 = -qJD(3) * t5 - t7 * qJD(4);
t2 = -t8 * qJD(4) + t14 * t13;
t1 = t2 * t8 + t3 * t7 + (-t4 * t8 + t5 * t7) * qJD(3);
t16 = [0, 0, 0, 0, t15, qJ(2) * t15, -0.2e1 * t7 * t12, 0.2e1 * (t7 ^ 2 - t8 ^ 2) * qJD(3), 0, 0, 0, 0.2e1 * qJD(2) * t7 + 0.2e1 * t8 * t11, 0.2e1 * qJD(2) * t8 - 0.2e1 * t7 * t11, -0.2e1 * t1, -0.2e1 * t4 * t3 - 0.2e1 * t5 * t2 + 0.2e1 * (t7 * pkin(3) + qJ(2)) * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, -t13, -t12, 0, -t9 * t13, -t9 * t12, t10, t2 * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, -t12, 0, -t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t16;
