% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RPPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x16]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPPRP3_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP3_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP3_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP3_inertiaDJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:51:02
% EndTime: 2019-12-31 17:51:02
% DurationCPUTime: 0.11s
% Computational Cost: add. (94->28), mult. (156->49), div. (0->0), fcn. (100->4), ass. (0->17)
t17 = 2 * qJD(3);
t7 = -cos(pkin(7)) * pkin(1) - pkin(2) - pkin(6);
t16 = qJ(5) - t7;
t10 = sin(qJ(4));
t15 = qJD(4) * t10;
t11 = cos(qJ(4));
t14 = qJD(4) * t11;
t13 = pkin(4) * t15;
t12 = pkin(4) * t14;
t8 = sin(pkin(7)) * pkin(1) + qJ(3);
t5 = t16 * t11;
t6 = qJD(3) + t12;
t4 = t16 * t10;
t3 = -qJD(4) * t5 - t10 * qJD(5);
t2 = -t11 * qJD(5) + t16 * t15;
t1 = t3 * t10 + t2 * t11 + (t10 * t5 - t11 * t4) * qJD(4);
t9 = [0, 0, 0, 0, 0, t17, t8 * t17, -0.2e1 * t10 * t14, 0.2e1 * (t10 ^ 2 - t11 ^ 2) * qJD(4), 0, 0, 0, 0.2e1 * qJD(3) * t10 + 0.2e1 * t8 * t14, 0.2e1 * qJD(3) * t11 - 0.2e1 * t8 * t15, -0.2e1 * t1, -0.2e1 * t4 * t3 - 0.2e1 * t5 * t2 + 0.2e1 * (t10 * pkin(4) + t8) * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2 * t10 + t3 * t11 + (t10 * t4 + t11 * t5) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, -t14, 0, -t7 * t15, -t7 * t14, t13, t2 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, t15, 0, -t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, -t14, 0, -t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t9;
