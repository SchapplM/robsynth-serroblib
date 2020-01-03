% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S4PRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x12]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4PRPR5_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR5_inertiaDJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR5_inertiaDJ_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR5_inertiaDJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:23:17
% EndTime: 2019-12-31 16:23:18
% DurationCPUTime: 0.12s
% Computational Cost: add. (31->17), mult. (93->38), div. (0->0), fcn. (84->6), ass. (0->16)
t9 = sin(qJ(4));
t15 = t9 * qJD(4);
t11 = cos(qJ(4));
t14 = t11 * qJD(4);
t13 = 0.2e1 * t14;
t10 = sin(qJ(2));
t12 = cos(qJ(2));
t7 = sin(pkin(7));
t8 = cos(pkin(7));
t4 = t8 * t10 + t7 * t12;
t3 = t7 * t10 - t8 * t12;
t6 = -t8 * pkin(2) - pkin(3);
t5 = t7 * pkin(2) + pkin(5);
t2 = t3 * qJD(2);
t1 = t4 * qJD(2);
t16 = [0, 0, 0, 0, 0.2e1 * t3 * t1 - 0.2e1 * t4 * t2, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t10 * qJD(2), -t12 * qJD(2), (-t1 * t8 - t2 * t7) * pkin(2), 0, 0, 0, 0, 0, -t1 * t11 + t3 * t15, t1 * t9 + t3 * t14; 0, 0, 0, 0, 0, t9 * t13, 0.2e1 * (t11 ^ 2 - t9 ^ 2) * qJD(4), 0, 0, 0, 0.2e1 * t6 * t15, t6 * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4 * t14 + t9 * t2, t11 * t2 + t4 * t15; 0, 0, 0, 0, 0, 0, 0, t14, -t15, 0, -t5 * t14, t5 * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, -t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t16;
