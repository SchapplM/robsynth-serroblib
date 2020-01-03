% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S4RPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2,theta3]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x15]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:38
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4RPPR3_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR3_inertiaDJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR3_inertiaDJ_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPPR3_inertiaDJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:37:56
% EndTime: 2019-12-31 16:37:56
% DurationCPUTime: 0.11s
% Computational Cost: add. (37->16), mult. (108->33), div. (0->0), fcn. (90->6), ass. (0->16)
t11 = sin(pkin(7));
t12 = cos(pkin(7));
t13 = sin(qJ(4));
t14 = cos(qJ(4));
t15 = t13 * t11 - t14 * t12;
t1 = t15 * qJD(4);
t20 = -0.2e1 * t1;
t8 = sin(pkin(6)) * pkin(1) + qJ(3);
t19 = pkin(5) + t8;
t16 = 0.2e1 * (t11 ^ 2 + t12 ^ 2) * qJD(3);
t5 = t14 * t11 + t13 * t12;
t6 = -cos(pkin(6)) * pkin(1) - pkin(2) - t12 * pkin(3);
t4 = t19 * t12;
t3 = t19 * t11;
t2 = t5 * qJD(4);
t7 = [0, 0, 0, 0, 0, 0, t16, t8 * t16, t5 * t20, 0.2e1 * t1 * t15 - 0.2e1 * t5 * t2, 0, 0, 0, 0.2e1 * t6 * t2, t6 * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, -t2, 0, (t13 * t3 - t14 * t4) * qJD(4) - t5 * qJD(3), (t13 * t4 + t14 * t3) * qJD(4) + t15 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t7;
