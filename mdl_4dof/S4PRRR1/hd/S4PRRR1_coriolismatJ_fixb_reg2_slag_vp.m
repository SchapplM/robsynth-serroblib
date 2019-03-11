% Calculate inertial parameters regressor of coriolis matrix for
% S4PRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% 
% Output:
% cmat_reg [(4*4)x(4*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4PRRR1_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR1_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR1_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR1_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:25:22
% EndTime: 2019-03-08 18:25:22
% DurationCPUTime: 0.20s
% Computational Cost: add. (177->41), mult. (445->68), div. (0->0), fcn. (314->4), ass. (0->39)
t26 = sin(qJ(4));
t7 = t26 * pkin(3);
t29 = cos(qJ(3));
t49 = t29 * pkin(2);
t27 = sin(qJ(3));
t48 = t26 * t27;
t47 = t26 * t29;
t28 = cos(qJ(4));
t46 = t28 * t27;
t45 = t28 * t29;
t44 = pkin(2) * qJD(2);
t43 = pkin(2) * qJD(3);
t42 = pkin(3) * qJD(3);
t41 = pkin(3) * qJD(4);
t33 = pkin(3) + t49;
t21 = t28 * t33;
t13 = pkin(2) * t48 - t21;
t30 = t26 * t33;
t14 = pkin(2) * t46 + t30;
t17 = (t46 + t47) * pkin(2);
t18 = (t45 - t48) * pkin(2);
t2 = t13 * t17 + t14 * t18;
t40 = t2 * qJD(2);
t39 = t7 * qJD(2);
t8 = t21 / 0.2e1 + (-t49 / 0.2e1 + pkin(3) / 0.2e1) * t28;
t38 = t8 * qJD(2);
t37 = t13 * qJD(2);
t36 = t14 * qJD(2);
t35 = t17 * qJD(2);
t34 = t18 * qJD(2);
t32 = pkin(2) * (-qJD(2) - qJD(3));
t31 = pkin(3) * (-qJD(3) - qJD(4));
t16 = t18 * qJD(3);
t15 = t17 * qJD(3);
t12 = t14 * qJD(4);
t11 = t13 * qJD(4);
t6 = -t28 * pkin(3) / 0.2e1 - t21 / 0.2e1 + (t48 - t45 / 0.2e1) * pkin(2);
t5 = -t7 / 0.2e1 - t30 / 0.2e1 + (-t46 - t47 / 0.2e1) * pkin(2);
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t27 * t43, -t29 * t43, 0, 0, 0, 0, 0, 0, 0, 0, -t15 - t12, -t16 + t11, 0, t2 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27 * t32, t29 * t32, 0, 0, 0, 0, 0, 0, 0, 0, t5 * qJD(4) - t15 - t35, t6 * qJD(4) - t16 - t34, 0, t40 + (-t17 * t28 + t18 * t26) * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5 * qJD(3) - t12 - t36, t6 * qJD(3) + t11 + t37, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27 * t44, t29 * t44, 0, 0, 0, 0, 0, 0, 0, 0, -t7 * qJD(4) + t35, -t8 * qJD(4) + t34, 0, -t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26 * t41, -t28 * t41, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26 * t31 - t39, t28 * t31 - t38, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7 * qJD(3) + t36, t8 * qJD(3) - t37, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26 * t42 + t39, t28 * t42 + t38, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t1;
