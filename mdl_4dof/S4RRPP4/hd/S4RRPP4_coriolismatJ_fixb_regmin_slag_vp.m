% Calculate minimal parameter regressor of coriolis matrix for
% S4RRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
% 
% Output:
% cmat_reg [(4*%NQJ)%x18]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4RRPP4_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP4_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP4_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP4_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:59:19
% EndTime: 2019-12-31 16:59:20
% DurationCPUTime: 0.29s
% Computational Cost: add. (224->69), mult. (409->86), div. (0->0), fcn. (301->2), ass. (0->57)
t34 = sin(qJ(2));
t29 = t34 * qJ(3);
t35 = cos(qJ(2));
t61 = pkin(2) + pkin(3);
t64 = t61 * t35 + t29;
t37 = -t35 * pkin(2) - t29;
t45 = t35 * qJD(3);
t63 = qJD(2) * t37 + t45;
t62 = pkin(5) - qJ(4);
t10 = pkin(1) + t64;
t56 = t35 * qJ(3);
t11 = -t61 * t34 + t56;
t1 = t10 * t11;
t59 = t1 * qJD(1);
t2 = t10 * t35 + t11 * t34;
t58 = t2 * qJD(1);
t3 = -t10 * t34 + t11 * t35;
t57 = t3 * qJD(1);
t13 = -pkin(1) + t37;
t15 = t34 * pkin(2) - t56;
t4 = t13 * t35 + t15 * t34;
t55 = t4 * qJD(1);
t5 = -t13 * t34 + t15 * t35;
t54 = t5 * qJD(1);
t14 = t62 * t34;
t16 = t62 * t35;
t7 = t14 * t34 + t16 * t35;
t53 = t7 * qJD(1);
t44 = -pkin(2) / 0.2e1 - pkin(3) / 0.2e1;
t8 = t56 + (-t61 / 0.2e1 + t44) * t34;
t52 = t8 * qJD(1);
t51 = qJD(2) * t34;
t27 = qJD(2) * t35;
t50 = t16 * qJD(2);
t32 = t34 ^ 2;
t33 = t35 ^ 2;
t17 = t33 + t32;
t49 = t17 * qJD(1);
t18 = t33 - t32;
t48 = t18 * qJD(1);
t26 = t34 * qJD(1);
t47 = t34 * qJD(3);
t46 = t35 * qJD(1);
t43 = pkin(1) * t26;
t42 = pkin(1) * t46;
t41 = pkin(5) * t51;
t40 = pkin(5) * t27;
t39 = t13 * t15 * qJD(1);
t38 = t13 * t26;
t31 = qJ(3) * qJD(3);
t30 = qJD(2) * qJ(3);
t25 = t32 * qJD(1);
t24 = t32 * qJD(3);
t20 = t34 * t46;
t19 = t34 * t45;
t9 = (t61 / 0.2e1 + t44) * t34;
t6 = [0, 0, 0, t34 * t27, t18 * qJD(2), 0, 0, 0, -pkin(1) * t51, -pkin(1) * t27, -t5 * qJD(2) + t19, 0, -t4 * qJD(2) + t24, (qJD(2) * t15 - t47) * t13, t3 * qJD(2) + t19, t2 * qJD(2) + t24, t17 * qJD(4), t1 * qJD(2) - t7 * qJD(4) + t10 * t47; 0, 0, 0, t20, t48, t27, -t51, 0, -t40 - t43, t41 - t42, -t40 - t54, t63, -t41 - t55, t63 * pkin(5) + t39, -t50 + t57, -t14 * qJD(2) + t58, t64 * qJD(2) - t45, t59 + (-t14 * qJ(3) - t16 * t61) * qJD(2) + t16 * qJD(3) + t9 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, t27, t25, -t38 + t40, t20, t25, -t27, t10 * t26 + t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, t9 * qJD(2) - t53; 0, 0, 0, -t20, -t48, 0, 0, 0, t43, t42, t54, 0, t55, -t39, t34 * qJD(4) - t57, -t35 * qJD(4) - t58, 0, -t8 * qJD(4) - t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3), t31, 0, qJD(3), 0, t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2), t30, 0, qJD(2), 0, t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, -t46, 0, -t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20, 0, -t25, t38, -t20, -t25, 0, (-qJD(1) * t10 - qJD(4)) * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(2), -t30, 0, -qJD(2), 0, -t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t51, t27, -t49, t8 * qJD(2) + t47 + t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, t46, 0, t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t6;
