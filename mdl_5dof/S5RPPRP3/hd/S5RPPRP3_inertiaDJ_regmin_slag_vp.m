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
% MMD_reg [((5+1)*5/2)x18]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 17:04
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
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
% StartTime: 2021-01-15 17:04:41
% EndTime: 2021-01-15 17:04:42
% DurationCPUTime: 0.14s
% Computational Cost: add. (114->35), mult. (192->57), div. (0->0), fcn. (124->4), ass. (0->18)
t18 = 2 * qJD(3);
t8 = -cos(pkin(7)) * pkin(1) - pkin(2) - pkin(6);
t17 = qJ(5) - t8;
t12 = sin(qJ(4));
t10 = t12 * qJD(4);
t13 = cos(qJ(4));
t16 = t13 * qJD(4);
t15 = pkin(4) * t10;
t14 = pkin(4) * t16;
t5 = t17 * t13;
t9 = sin(pkin(7)) * pkin(1) + qJ(3);
t7 = qJD(3) + t14;
t6 = t12 * pkin(4) + t9;
t4 = t17 * t12;
t3 = -qJD(4) * t5 - t12 * qJD(5);
t2 = -t13 * qJD(5) + t17 * t10;
t1 = t3 * t12 + t2 * t13 + (t12 * t5 - t13 * t4) * qJD(4);
t11 = [0, 0, 0, 0, 0, t18, t9 * t18, -0.2e1 * t12 * t16, 0.2e1 * (t12 ^ 2 - t13 ^ 2) * qJD(4), 0, 0, 0, 0.2e1 * qJD(3) * t12 + 0.2e1 * t9 * t16, 0.2e1 * qJD(3) * t13 - 0.2e1 * t9 * t10, 0.2e1 * t7 * t12 + 0.2e1 * t6 * t16, -0.2e1 * t6 * t10 + 0.2e1 * t7 * t13, -0.2e1 * t1, -0.2e1 * t5 * t2 - 0.2e1 * t4 * t3 + 0.2e1 * t6 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2 * t12 + t3 * t13 + (t12 * t4 + t13 * t5) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, -t16, 0, -t8 * t10, -t8 * t16, t2, -t3, t15, t2 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, t10, -t16, t10, 0, -t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, -t16, -t10, -t16, 0, -t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, -t10, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t11;
