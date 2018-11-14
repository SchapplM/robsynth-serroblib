% Calculate minimal parameter regressor of coriolis matrix for
% S4RRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% cmat_reg [(4*%NQJ)%x10]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:55
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function cmat_reg = S4RRRP1_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP1_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP1_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP1_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:54:31
% EndTime: 2018-11-14 13:54:31
% DurationCPUTime: 0.18s
% Computational Cost: add. (256->57), mult. (595->92), div. (0->0), fcn. (415->4), ass. (0->51)
t33 = sin(qJ(2));
t34 = cos(qJ(3));
t56 = t34 * t33;
t32 = sin(qJ(3));
t35 = cos(qJ(2));
t57 = t32 * t35;
t22 = (t56 + t57) * pkin(1);
t61 = t22 * pkin(3);
t9 = t32 * pkin(2);
t60 = t34 * pkin(2);
t59 = t35 * pkin(1);
t58 = t32 * t33;
t55 = t34 * t35;
t54 = pkin(1) * qJD(1);
t53 = pkin(1) * qJD(2);
t52 = pkin(2) * qJD(2);
t51 = pkin(2) * qJD(3);
t41 = pkin(2) + t59;
t26 = t34 * t41;
t18 = pkin(1) * t58 - t26;
t15 = pkin(3) - t18;
t38 = t32 * t41;
t19 = pkin(1) * t56 + t38;
t3 = (-t15 - t18) * t19;
t50 = t3 * qJD(1);
t23 = (t55 - t58) * pkin(1);
t4 = -t15 * t22 + t19 * t23;
t49 = t4 * qJD(1);
t48 = t9 * qJD(1);
t10 = t26 / 0.2e1 + (-t59 / 0.2e1 + pkin(2) / 0.2e1) * t34;
t47 = t10 * qJD(1);
t46 = t18 * qJD(1);
t45 = t19 * qJD(1);
t14 = t19 * qJD(3);
t44 = t22 * qJD(1);
t43 = t23 * qJD(1);
t42 = t32 * t51;
t40 = pkin(1) * (-qJD(1) - qJD(2));
t39 = pkin(2) * (-qJD(2) - qJD(3));
t29 = pkin(3) + t60;
t16 = (-t29 + t60) * t9;
t36 = (t19 * t34 / 0.2e1 + (-t18 / 0.2e1 - t15 / 0.2e1) * t32) * pkin(2) - t19 * t29 / 0.2e1;
t2 = t61 / 0.2e1 + t36;
t37 = -t2 * qJD(1) - t16 * qJD(2);
t21 = t23 * qJD(2);
t20 = t22 * qJD(2);
t13 = t18 * qJD(3);
t6 = -t60 / 0.2e1 - t26 / 0.2e1 + (t58 - t55 / 0.2e1) * pkin(1);
t5 = -t9 / 0.2e1 - t38 / 0.2e1 + (-t56 - t57 / 0.2e1) * pkin(1);
t1 = -t61 / 0.2e1 + t36;
t7 = [0, 0, 0, 0, -t33 * t53, -t35 * t53, 0, -t20 - t14, -t21 + t13, t4 * qJD(2) + t3 * qJD(3); 0, 0, 0, 0, t33 * t40, t35 * t40, 0, t5 * qJD(3) - t20 - t44, t6 * qJD(3) - t21 - t43, t49 + (-t22 * t29 + t23 * t9) * qJD(2) + t1 * qJD(3); 0, 0, 0, 0, 0, 0, 0, t5 * qJD(2) - t14 - t45, t6 * qJD(2) + t13 + t46, -pkin(3) * t14 + t1 * qJD(2) + t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, t33 * t54, t35 * t54, 0, -t9 * qJD(3) + t44, -t10 * qJD(3) + t43, t2 * qJD(3) - t49; 0, 0, 0, 0, 0, 0, 0, -t42, -t34 * t51, t16 * qJD(3); 0, 0, 0, 0, 0, 0, 0, t32 * t39 - t48, t34 * t39 - t47, -pkin(3) * t42 - t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, t9 * qJD(2) + t45, t10 * qJD(2) - t46, -t2 * qJD(2) - t50; 0, 0, 0, 0, 0, 0, 0, t32 * t52 + t48, t34 * t52 + t47, t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t7;
