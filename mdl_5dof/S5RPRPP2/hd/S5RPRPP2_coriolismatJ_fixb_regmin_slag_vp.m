% Calculate minimal parameter regressor of coriolis matrix for
% S5RPRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x19]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RPRPP2_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP2_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP2_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP2_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:11:14
% EndTime: 2019-12-31 18:11:15
% DurationCPUTime: 0.40s
% Computational Cost: add. (334->79), mult. (559->90), div. (0->0), fcn. (431->4), ass. (0->64)
t45 = cos(qJ(3));
t71 = pkin(3) + pkin(4);
t76 = t71 * t45;
t44 = sin(qJ(3));
t38 = t44 * pkin(3);
t64 = t45 * qJ(4);
t20 = -t38 + t64;
t32 = t44 * qJD(4);
t75 = t20 * qJD(3) + t32;
t26 = sin(pkin(7)) * pkin(1) + pkin(6);
t74 = -qJ(5) + t26;
t37 = t44 * qJ(4);
t54 = t45 * qJD(4);
t69 = t45 * pkin(3);
t73 = (-t37 - t69) * qJD(3) + t54;
t72 = -t38 / 0.2e1;
t70 = t44 * pkin(4);
t68 = t71 * t44;
t27 = -cos(pkin(7)) * pkin(1) - pkin(2);
t47 = -t27 + t37;
t13 = t47 + t76;
t19 = t20 - t70;
t1 = t13 * t19;
t67 = t1 * qJD(1);
t3 = t13 * t45 + t19 * t44;
t66 = t3 * qJD(1);
t4 = -t13 * t44 + t19 * t45;
t65 = t4 * qJD(1);
t17 = t74 * t44;
t18 = t74 * t45;
t7 = t17 * t44 + t18 * t45;
t63 = t7 * qJD(1);
t16 = -t47 - t69;
t8 = t16 * t45 - t20 * t44;
t62 = t8 * qJD(1);
t9 = -t16 * t44 - t20 * t45;
t61 = t9 * qJD(1);
t10 = t64 + t72 + (-pkin(4) / 0.2e1 - t71 / 0.2e1) * t44;
t60 = t10 * qJD(1);
t59 = t18 * qJD(3);
t41 = t44 ^ 2;
t42 = t45 ^ 2;
t21 = t41 + t42;
t57 = t21 * qJD(1);
t22 = t42 - t41;
t56 = t22 * qJD(1);
t34 = t44 * qJD(1);
t33 = t44 * qJD(3);
t55 = t45 * qJD(1);
t35 = t45 * qJD(3);
t53 = t16 * t20 * qJD(1);
t52 = t16 * t34;
t51 = t27 * t34;
t50 = t27 * t55;
t49 = t26 * t33;
t48 = t26 * t35;
t40 = qJ(4) * qJD(4);
t39 = qJD(3) * qJ(4);
t31 = t41 * qJD(1);
t30 = t41 * qJD(4);
t24 = t44 * t55;
t23 = t44 * t54;
t15 = t68 / 0.2e1 + t72 - t70 / 0.2e1;
t2 = [0, 0, 0, 0, t44 * t35, t22 * qJD(3), 0, 0, 0, t27 * t33, t27 * t35, -t9 * qJD(3) + t23, 0, -t8 * qJD(3) + t30, -t75 * t16, t4 * qJD(3) + t23, t3 * qJD(3) + t30, t21 * qJD(5), t1 * qJD(3) - t7 * qJD(5) + t13 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, t24, t56, t35, -t33, 0, -t48 + t51, t49 + t50, -t48 - t61, t73, -t49 - t62, t73 * t26 - t53, -t59 + t65, -t17 * qJD(3) + t66, (t37 + t76) * qJD(3) - t54, t67 + (-t17 * qJ(4) - t18 * t71) * qJD(3) + t18 * qJD(4) + t15 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, t35, t31, t48 - t52, t24, t31, -t35, t13 * t34 + t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, t15 * qJD(3) - t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33, -t35, -t33, 0, t35, t75, -t33, t35, 0, (t64 - t68) * qJD(3) + t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, 0, 0, 0, t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, -t24, -t56, 0, 0, 0, -t51, -t50, t61, 0, t62, t53, t44 * qJD(5) - t65, -t45 * qJD(5) - t66, 0, -t10 * qJD(5) - t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4), t40, 0, qJD(4), 0, t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3), t39, 0, qJD(3), 0, t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, -t55, 0, -t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, 0, -t31, t52, -t24, -t31, 0, (-qJD(1) * t13 - qJD(5)) * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3), -t39, 0, -qJD(3), 0, -t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33, t35, -t57, t10 * qJD(3) + t32 + t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, t55, 0, t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t2;
