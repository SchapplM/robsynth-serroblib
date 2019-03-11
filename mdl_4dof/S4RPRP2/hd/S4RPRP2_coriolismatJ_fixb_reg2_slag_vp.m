% Calculate inertial parameters regressor of coriolis matrix for
% S4RPRP2
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
% cmat_reg [(4*4)x(4*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4RPRP2_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP2_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP2_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP2_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:30:55
% EndTime: 2019-03-08 18:30:55
% DurationCPUTime: 0.26s
% Computational Cost: add. (270->38), mult. (365->42), div. (0->0), fcn. (265->2), ass. (0->36)
t48 = -pkin(3) / 0.2e1;
t47 = -pkin(1) - pkin(2);
t28 = sin(qJ(3));
t29 = cos(qJ(3));
t20 = t28 * qJ(2) - t29 * t47;
t19 = -pkin(3) - t20;
t46 = t19 * t28;
t45 = t20 * t28;
t21 = t29 * qJ(2) + t28 * t47;
t44 = t21 * t29;
t1 = (t48 + t20 / 0.2e1 + t19 / 0.2e1) * t28;
t43 = t1 * qJD(1);
t3 = (-t19 - t20) * t21;
t42 = t3 * qJD(1);
t7 = t44 - t46;
t39 = t7 * qJD(1);
t8 = t44 + t45;
t38 = t8 * qJD(1);
t37 = qJD(3) * t28;
t36 = t21 * qJD(3);
t35 = t28 * qJD(1);
t34 = t28 * qJD(2);
t33 = t29 * qJD(1);
t32 = t29 * qJD(2);
t31 = qJ(2) * qJD(1);
t30 = qJD(1) - qJD(3);
t23 = -qJD(3) * t29 + t33;
t22 = t35 - t37;
t14 = -t21 * qJD(1) - t34;
t13 = t20 * qJD(1) - t32;
t12 = -t20 * qJD(3) + t32;
t11 = t34 + t36;
t10 = t30 * t21;
t9 = t30 * t20;
t2 = -t45 / 0.2e1 - t46 / 0.2e1 + t28 * t48;
t4 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2), qJ(2) * qJD(2), 0, 0, 0, 0, 0, 0, t11, t12, 0, t8 * qJD(2), 0, 0, 0, 0, 0, 0, t11, t12, 0, t7 * qJD(2) + t3 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(1), t31, 0, 0, 0, 0, 0, 0, t35, t33, 0, t38, 0, 0, 0, 0, 0, 0, t35, t33, 0, t2 * qJD(3) + t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t9, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t9, 0, -pkin(3) * t36 + t2 * qJD(2) + t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(1), -t31, 0, 0, 0, 0, 0, 0, -t22, -t23, 0, -t38, 0, 0, 0, 0, 0, 0, -t22, -t23, 0, -t1 * qJD(3) - t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, t23, 0, 0, 0, 0, 0, 0, 0, 0, t22, t23, 0, -pkin(3) * t37 - t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, t13, 0, 0, 0, 0, 0, 0, 0, 0, t14, t13, 0, t1 * qJD(2) - t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35, -t33, 0, 0, 0, 0, 0, 0, 0, 0, -t35, -t33, 0, t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t4;
