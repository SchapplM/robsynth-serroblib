% Calculate minimal parameter regressor of coriolis matrix for
% S4PRRR4
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
% cmat_reg [(4*%NQJ)%x18]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4PRRR4_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR4_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR4_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR4_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:32:41
% EndTime: 2019-12-31 16:32:41
% DurationCPUTime: 0.22s
% Computational Cost: add. (189->44), mult. (421->60), div. (0->0), fcn. (414->4), ass. (0->39)
t36 = qJD(3) + qJD(4);
t26 = sin(qJ(4));
t27 = sin(qJ(3));
t28 = cos(qJ(3));
t50 = cos(qJ(4));
t18 = t26 * t27 - t50 * t28;
t56 = t36 * t18;
t53 = pkin(5) + pkin(6);
t21 = t53 * t27;
t22 = t53 * t28;
t55 = t36 * (t26 * t21 - t50 * t22);
t54 = t36 * (t50 * t21 + t26 * t22);
t52 = pkin(3) * t26;
t51 = pkin(3) * t27;
t20 = t26 * t28 + t50 * t27;
t3 = t18 ^ 2 - t20 ^ 2;
t49 = t3 * qJD(2);
t25 = -t28 * pkin(3) - pkin(2);
t6 = t18 * t51 + t25 * t20;
t44 = t6 * qJD(2);
t7 = -t25 * t18 + t20 * t51;
t43 = t7 * qJD(2);
t42 = qJD(2) * t25;
t41 = qJD(2) * t28;
t40 = qJD(3) * t27;
t39 = qJD(3) * t28;
t38 = qJD(4) * t25;
t23 = -t27 ^ 2 + t28 ^ 2;
t37 = t23 * qJD(2);
t35 = pkin(2) * t27 * qJD(2);
t34 = pkin(2) * t41;
t33 = t18 * t42;
t32 = t20 * t42;
t31 = t27 * t41;
t30 = t50 * qJD(3);
t29 = t50 * qJD(4);
t10 = t36 * t20;
t11 = t20 * t18 * qJD(2);
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t40, -t39, 0, 0, 0, 0, 0, -t10, t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, t27 * t39, t23 * qJD(3), 0, 0, 0, -pkin(2) * t40, -pkin(2) * t39, -t18 * t10, t36 * t3, 0, 0, 0, t6 * qJD(3) + t20 * t38, t7 * qJD(3) - t18 * t38; 0, 0, 0, 0, t31, t37, t39, -t40, 0, -pkin(5) * t39 - t35, pkin(5) * t40 - t34, -t11, t49, -t56, -t10, 0, t44 + t55, t43 + t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, t49, -t56, -t10, 0, t32 + t55, -t33 + t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, -t31, -t37, 0, 0, 0, t35, t34, t11, -t49, 0, 0, 0, -t44, -t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(4) * t52, -pkin(3) * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36 * t52, (-t30 - t29) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, -t49, 0, 0, 0, -t32, t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3) * t52, pkin(3) * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t1;
