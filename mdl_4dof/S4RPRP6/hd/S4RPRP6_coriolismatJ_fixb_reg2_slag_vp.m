% Calculate inertial parameters regressor of coriolis matrix for
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
% cmat_reg [(4*4)x(4*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4RPRP6_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP6_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP6_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP6_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:46:16
% EndTime: 2019-12-31 16:46:17
% DurationCPUTime: 0.36s
% Computational Cost: add. (189->61), mult. (327->58), div. (0->0), fcn. (228->2), ass. (0->46)
t35 = -pkin(1) - pkin(5);
t55 = -qJ(4) + t35;
t33 = sin(qJ(3));
t54 = t33 * pkin(3);
t21 = qJ(2) + t54;
t34 = cos(qJ(3));
t1 = t21 * t34 * pkin(3);
t53 = t1 * qJD(1);
t10 = t55 * t33;
t11 = t55 * t34;
t4 = t10 * t33 + t11 * t34;
t52 = t4 * qJD(1);
t6 = (t21 + t54) * t34;
t51 = t6 * qJD(1);
t32 = t34 ^ 2;
t9 = t32 * pkin(3) - t21 * t33;
t50 = t9 * qJD(1);
t49 = t10 * qJD(3);
t31 = t33 ^ 2;
t38 = -t31 / 0.2e1 - t32 / 0.2e1;
t13 = -0.1e1 / 0.2e1 + t38;
t48 = t13 * qJD(1);
t17 = t31 + t32;
t47 = t17 * qJD(1);
t18 = t31 - t32;
t46 = t18 * qJD(1);
t45 = t21 * qJD(1);
t44 = t33 * qJD(3);
t25 = t34 * qJD(1);
t43 = t34 * qJD(3);
t42 = t34 * qJD(4);
t41 = qJ(2) * qJD(3);
t29 = qJD(1) * qJ(2);
t40 = pkin(3) * t44;
t39 = pkin(3) * t25;
t37 = t33 * t29;
t36 = t34 * t29;
t30 = qJ(2) * qJD(2);
t27 = qJD(2) * t34;
t26 = qJD(2) * t33;
t24 = t33 * qJD(1);
t20 = t33 * t43;
t19 = t33 * t25;
t14 = t18 * qJD(3);
t12 = 0.1e1 / 0.2e1 + t38;
t2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2), t30, -t20, t14, 0, t20, 0, 0, t34 * t41 + t26, -t33 * t41 + t27, 0, t30, -t20, t14, 0, t20, 0, 0, t6 * qJD(3) + t26, t9 * qJD(3) + t27, t17 * qJD(4), t21 * qJD(2) + t1 * qJD(3) - t4 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(1), t29, 0, 0, 0, 0, 0, 0, t24, t25, 0, t29, 0, 0, 0, 0, 0, 0, t24, t25, 0, t12 * qJD(4) + t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, t46, -t44, t19, -t43, 0, -t35 * t44 + t36, -t35 * t43 - t37, 0, 0, -t19, t46, -t44, t19, -t43, 0, -t49 + t51, -t11 * qJD(3) + t50, t40, -pkin(3) * t49 + t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, t12 * qJD(2) - t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(1), -t29, 0, 0, 0, 0, 0, 0, -t24, -t25, 0, -t29, 0, 0, 0, 0, 0, 0, -t24, -t25, 0, t13 * qJD(4) - t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44, -t43, 0, 0, 0, 0, 0, 0, 0, 0, -t44, -t43, 0, -t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, -t46, 0, -t19, 0, 0, -t36, t37, 0, 0, t19, -t46, 0, -t19, 0, 0, -t42 - t51, t33 * qJD(4) - t50, 0, -pkin(3) * t42 - t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t25, t24, 0, -t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, -t44, -t47, pkin(3) * t43 - t13 * qJD(2) + t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, -t24, 0, t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t2;
