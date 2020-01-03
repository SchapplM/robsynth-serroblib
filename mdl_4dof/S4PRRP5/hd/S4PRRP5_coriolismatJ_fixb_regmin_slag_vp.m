% Calculate minimal parameter regressor of coriolis matrix for
% S4PRRP5
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
% cmat_reg [(4*%NQJ)%x13]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:29
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4PRRP5_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP5_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP5_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP5_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:29:23
% EndTime: 2019-12-31 16:29:24
% DurationCPUTime: 0.17s
% Computational Cost: add. (119->36), mult. (319->66), div. (0->0), fcn. (241->4), ass. (0->35)
t21 = sin(qJ(3));
t19 = t21 ^ 2;
t23 = cos(qJ(3));
t20 = t23 ^ 2;
t13 = t20 + t19;
t42 = pkin(3) * t21;
t41 = -qJ(4) - pkin(5);
t17 = -t23 * pkin(3) - pkin(2);
t1 = t17 * t42;
t40 = t1 * qJD(2);
t22 = sin(qJ(2));
t24 = cos(qJ(2));
t6 = (-0.1e1 + t13) * t24 * t22;
t39 = t6 * qJD(1);
t38 = qJD(2) * t21;
t37 = qJD(2) * t23;
t36 = t13 * qJD(2);
t14 = t20 - t19;
t35 = t14 * qJD(2);
t34 = t21 * qJD(3);
t33 = t22 * qJD(2);
t32 = t23 * qJD(3);
t31 = t24 * qJD(2);
t30 = pkin(2) * t38;
t29 = pkin(2) * t37;
t28 = pkin(3) * t38;
t27 = t21 * t37;
t26 = t22 * t32;
t12 = t41 * t23;
t5 = -t12 * t23 - t19 * t41;
t7 = (0.1e1 / 0.2e1 - t20 / 0.2e1 - t19 / 0.2e1) * t22;
t25 = -t7 * qJD(1) + t5 * qJD(2);
t8 = (0.1e1 + t13) * t22 / 0.2e1;
t2 = t24 * t42;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6 * qJD(2); 0, 0, -t33, -t31, 0, 0, 0, 0, 0, -t23 * t33 - t24 * t34, t21 * t33 - t24 * t32, t13 * t31, t39 + (t22 * t17 + t5 * t24) * qJD(2) - t2 * qJD(3) + t8 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21 * t31 - t26, t22 * t34 - t23 * t31, 0, -pkin(3) * t26 - t2 * qJD(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8 * qJD(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7 * qJD(4) - t39; 0, 0, 0, 0, t21 * t32, t14 * qJD(3), 0, 0, 0, -pkin(2) * t34, -pkin(2) * t32, t13 * qJD(4), t1 * qJD(3) + t5 * qJD(4); 0, 0, 0, 0, t27, t35, t32, -t34, 0, -pkin(5) * t32 - t30, pkin(5) * t34 - t29, -pkin(3) * t32, t12 * pkin(3) * qJD(3) + t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, -t27, -t35, 0, 0, 0, t30, t29, 0, -qJD(4) * t42 - t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7 * qJD(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36, pkin(3) * t34 - t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t3;
