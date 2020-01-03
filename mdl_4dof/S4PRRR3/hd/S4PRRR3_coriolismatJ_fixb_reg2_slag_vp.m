% Calculate inertial parameters regressor of coriolis matrix for
% S4PRRR3
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
% Datum: 2019-12-31 16:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4PRRR3_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR3_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR3_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR3_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:31:41
% EndTime: 2019-12-31 16:31:42
% DurationCPUTime: 0.33s
% Computational Cost: add. (125->51), mult. (359->67), div. (0->0), fcn. (231->4), ass. (0->46)
t29 = cos(qJ(3));
t50 = t29 * pkin(2);
t20 = -pkin(3) - t50;
t55 = pkin(3) / 0.2e1 - t20 / 0.2e1;
t26 = sin(qJ(4));
t24 = t26 ^ 2;
t28 = cos(qJ(4));
t25 = t28 ^ 2;
t17 = t25 - t24;
t42 = qJD(2) + qJD(3);
t54 = t42 * t17;
t39 = -t50 / 0.2e1;
t53 = t39 - t55;
t49 = pkin(2) * qJD(2);
t48 = pkin(2) * qJD(3);
t47 = pkin(3) * qJD(3);
t27 = sin(qJ(3));
t19 = t27 * pkin(2) + pkin(6);
t36 = (t24 + t25) * t29;
t1 = (t19 * t36 + t20 * t27) * pkin(2);
t46 = t1 * qJD(2);
t45 = qJD(2) * t20;
t10 = pkin(2) * t36;
t44 = t10 * qJD(2);
t43 = t26 * qJD(4);
t23 = t28 * qJD(4);
t41 = t27 * t48;
t40 = t27 * t49;
t38 = t26 * t45;
t37 = t28 * t45;
t35 = pkin(2) * t42;
t34 = t26 * t40;
t33 = t27 * t35;
t32 = t39 + t55;
t2 = t32 * t26;
t31 = t2 * qJD(2) + t26 * t47;
t3 = t32 * t28;
t30 = t3 * qJD(2) + t28 * t47;
t18 = t26 * t23;
t16 = t26 * t41;
t13 = t17 * qJD(4);
t9 = t42 * t28 * t26;
t8 = t10 * qJD(3);
t5 = t53 * t28;
t4 = t53 * t26;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43, -t23, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41, -t29 * t48, 0, 0, t18, t13, 0, -t18, 0, 0, t20 * t43 - t28 * t41, t20 * t23 + t16, t8, t1 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33, -t29 * t35, 0, 0, t18, t13, 0, -t18, 0, 0, t4 * qJD(4) - t28 * t33, t5 * qJD(4) + t16 + t34, t8 + t44, t46 + (-pkin(3) * t27 + pkin(6) * t36) * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, t54, t23, -t9, -t43, 0, t4 * qJD(3) - t19 * t23 + t38, t5 * qJD(3) + t19 * t43 + t37, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, t29 * t49, 0, 0, t18, t13, 0, -t18, 0, 0, -t2 * qJD(4) + t28 * t40, -t3 * qJD(4) - t34, -t44, -t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, t13, 0, -t18, 0, 0, -pkin(3) * t43, -pkin(3) * t23, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, t54, t23, -t9, -t43, 0, -pkin(6) * t23 - t31, pkin(6) * t43 - t30, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, -t54, 0, t9, 0, 0, t2 * qJD(3) - t38, t3 * qJD(3) - t37, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, -t54, 0, t9, 0, 0, t31, t30, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t6;
