% Calculate minimal parameter regressor of coriolis matrix for
% S5PRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta4]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x17]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:38
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5PRPPR5_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR5_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR5_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR5_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:38:39
% EndTime: 2019-12-31 17:38:40
% DurationCPUTime: 0.23s
% Computational Cost: add. (131->45), mult. (264->70), div. (0->0), fcn. (282->6), ass. (0->36)
t44 = cos(qJ(2));
t43 = sin(qJ(2));
t19 = sin(pkin(8));
t20 = cos(pkin(8));
t23 = -pkin(2) - pkin(3);
t12 = t20 * qJ(3) + t19 * t23;
t11 = -t19 * qJ(3) + t20 * t23;
t4 = -t11 * t19 + t12 * t20;
t42 = t4 * qJD(2);
t21 = sin(qJ(5));
t41 = qJD(2) * t21;
t22 = cos(qJ(5));
t40 = qJD(2) * t22;
t39 = qJD(5) * t21;
t38 = qJD(5) * t22;
t13 = -t21 ^ 2 + t22 ^ 2;
t37 = t13 * qJD(2);
t36 = t19 * qJD(2);
t35 = t19 * qJD(3);
t34 = t20 * qJD(2);
t33 = t20 * qJD(3);
t32 = qJD(2) * qJ(3);
t31 = t22 * t34;
t30 = t21 * t36;
t29 = t21 * t34;
t28 = t21 * t40;
t27 = t22 * t36;
t26 = qJD(2) * t44;
t25 = qJD(2) * t43;
t9 = pkin(4) - t11;
t24 = qJD(2) * t9 - t33;
t10 = -pkin(6) + t12;
t8 = -t44 * t19 + t43 * t20;
t7 = -t43 * t19 - t44 * t20;
t2 = -t7 * t19 + t8 * t20;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t25, -t26, -t25, t26, (-t43 * pkin(2) + t44 * qJ(3)) * qJD(2) + t43 * qJD(3), -t8 * qJD(2), -t7 * qJD(2), (t8 * t11 - t7 * t12) * qJD(2) + t2 * qJD(3), 0, 0, 0, 0, 0, t7 * t39 - t8 * t40, t7 * t38 + t8 * t41; 0, 0, 0, 0, 0, 0, t25, 0, 0, t2 * qJD(2), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8 * t38 + t7 * t41, t8 * t39 + t7 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, qJD(3), qJ(3) * qJD(3), t35, t33, t4 * qJD(3), t21 * t38, t13 * qJD(5), 0, 0, 0, t22 * t35 - t9 * t39, -t21 * t35 - t9 * t38; 0, 0, 0, 0, 0, qJD(2), t32, t36, t34, t42, 0, 0, 0, 0, 0, t27, -t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, t37, -t38, t39, 0, -t10 * t38 - t9 * t41, t10 * t39 - t9 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -qJD(2), -t32, -t36, -t34, -t42, 0, 0, 0, 0, 0, t20 * t39 - t27, t20 * t38 + t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19 * t38 + t29, t19 * t39 + t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t39, -t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28, -t37, 0, 0, 0, t24 * t21, t24 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29, -t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t1;
