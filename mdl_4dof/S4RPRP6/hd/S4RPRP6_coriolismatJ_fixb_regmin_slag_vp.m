% Calculate minimal parameter regressor of coriolis matrix for
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
% cmat_reg [(4*%NQJ)%x17]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 10:27
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4RPRP6_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP6_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP6_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP6_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 10:27:43
% EndTime: 2021-01-15 10:27:44
% DurationCPUTime: 0.23s
% Computational Cost: add. (179->53), mult. (295->59), div. (0->0), fcn. (202->2), ass. (0->44)
t28 = -pkin(1) - pkin(5);
t51 = -qJ(4) + t28;
t26 = sin(qJ(3));
t50 = t26 * pkin(3);
t16 = qJ(2) + t50;
t27 = cos(qJ(3));
t1 = t16 * t27 * pkin(3);
t49 = t1 * qJD(1);
t10 = t51 * t26;
t11 = t51 * t27;
t4 = t10 * t26 + t11 * t27;
t48 = t4 * qJD(1);
t6 = (t16 + t50) * t27;
t47 = t6 * qJD(1);
t25 = t27 ^ 2;
t9 = t25 * pkin(3) - t16 * t26;
t46 = t9 * qJD(1);
t45 = qJD(3) * t26;
t44 = qJD(3) * t27;
t43 = qJD(3) * t28;
t42 = t10 * qJD(3);
t24 = t26 ^ 2;
t31 = -t24 / 0.2e1 - t25 / 0.2e1;
t13 = -0.1e1 / 0.2e1 + t31;
t41 = t13 * qJD(1);
t14 = t24 + t25;
t40 = t14 * qJD(1);
t15 = t24 - t25;
t39 = t15 * qJD(1);
t38 = t16 * qJD(1);
t20 = t27 * qJD(1);
t37 = t27 * qJD(4);
t36 = qJ(2) * qJD(3);
t35 = qJD(1) * qJ(2);
t34 = pkin(3) * t45;
t33 = pkin(3) * t20;
t32 = t26 * t20;
t30 = t26 * t35;
t29 = t27 * t35;
t22 = qJD(2) * t27;
t21 = qJD(2) * t26;
t19 = t26 * qJD(1);
t12 = 0.1e1 / 0.2e1 + t31;
t2 = [0, 0, 0, 0, qJD(2), qJ(2) * qJD(2), -t26 * t44, t15 * qJD(3), 0, 0, 0, t27 * t36 + t21, -t26 * t36 + t22, t6 * qJD(3) + t21, t9 * qJD(3) + t22, t14 * qJD(4), t16 * qJD(2) + t1 * qJD(3) - t4 * qJD(4); 0, 0, 0, 0, qJD(1), t35, 0, 0, 0, 0, 0, t19, t20, t19, t20, 0, t12 * qJD(4) + t38; 0, 0, 0, 0, 0, 0, -t32, t39, -t45, -t44, 0, -t26 * t43 + t29, -t27 * t43 - t30, -t42 + t47, -t11 * qJD(3) + t46, t34, -pkin(3) * t42 + t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, t12 * qJD(2) - t48; 0, 0, 0, 0, -qJD(1), -t35, 0, 0, 0, 0, 0, -t19, -t20, -t19, -t20, 0, t13 * qJD(4) - t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t45, -t44, -t45, -t44, 0, -t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41; 0, 0, 0, 0, 0, 0, t32, -t39, 0, 0, 0, -t29, t30, -t37 - t47, t26 * qJD(4) - t46, 0, -pkin(3) * t37 - t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20, t19, 0, -t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, -t45, -t40, pkin(3) * t44 - t13 * qJD(2) + t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, -t19, 0, t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t2;
