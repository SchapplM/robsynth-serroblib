% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RPPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x18]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPPPR4_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR4_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR4_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR4_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:45:15
% EndTime: 2019-12-31 17:45:16
% DurationCPUTime: 0.12s
% Computational Cost: add. (53->23), mult. (139->45), div. (0->0), fcn. (112->6), ass. (0->18)
t14 = sin(pkin(8));
t16 = cos(pkin(8));
t17 = sin(qJ(5));
t18 = cos(qJ(5));
t26 = -t17 * t14 + t18 * t16;
t8 = (t14 ^ 2 + t16 ^ 2) * qJD(4);
t25 = 2 * qJD(3);
t10 = -cos(pkin(7)) * pkin(1) - pkin(2) - qJ(4);
t24 = -pkin(6) + t10;
t11 = sin(pkin(7)) * pkin(1) + qJ(3);
t20 = t11 * qJD(3);
t5 = t18 * t14 + t17 * t16;
t7 = t14 * pkin(4) + t11;
t4 = t26 * qJD(5);
t3 = t5 * qJD(5);
t2 = t24 * t16;
t1 = t24 * t14;
t6 = [0, 0, 0, 0, 0, t25, 0.2e1 * t20, t14 * t25, t16 * t25, 0.2e1 * t8, -0.2e1 * t10 * t8 + 0.2e1 * t20, -0.2e1 * t26 * t3, -0.2e1 * t26 * t4 + 0.2e1 * t3 * t5, 0, 0, 0, 0.2e1 * qJD(3) * t5 + 0.2e1 * t7 * t4, 0.2e1 * qJD(3) * t26 - 0.2e1 * t7 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3), 0, 0, 0, 0, 0, t4, -t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, -t4, 0, (-t1 * t18 - t17 * t2) * qJD(5) - t26 * qJD(4), (t1 * t17 - t18 * t2) * qJD(5) + t5 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t6;
