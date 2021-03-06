% Calculate inertial parameters regressor of coriolis matrix for
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
% cmat_reg [(5*5)x(5*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RPRPP2_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP2_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP2_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP2_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:11:17
% EndTime: 2019-12-31 18:11:18
% DurationCPUTime: 0.84s
% Computational Cost: add. (350->93), mult. (611->90), div. (0->0), fcn. (477->4), ass. (0->66)
t49 = cos(qJ(3));
t75 = pkin(3) + pkin(4);
t80 = t75 * t49;
t48 = sin(qJ(3));
t42 = t48 * pkin(3);
t68 = t49 * qJ(4);
t20 = -t42 + t68;
t36 = t48 * qJD(4);
t79 = t20 * qJD(3) + t36;
t30 = sin(pkin(7)) * pkin(1) + pkin(6);
t78 = -qJ(5) + t30;
t41 = t48 * qJ(4);
t58 = t49 * qJD(4);
t73 = t49 * pkin(3);
t77 = (-t41 - t73) * qJD(3) + t58;
t76 = -t42 / 0.2e1;
t74 = t48 * pkin(4);
t72 = t75 * t48;
t31 = -cos(pkin(7)) * pkin(1) - pkin(2);
t51 = -t31 + t41;
t13 = t51 + t80;
t19 = t20 - t74;
t1 = t13 * t19;
t71 = t1 * qJD(1);
t3 = t13 * t49 + t19 * t48;
t70 = t3 * qJD(1);
t4 = -t13 * t48 + t19 * t49;
t69 = t4 * qJD(1);
t17 = t78 * t48;
t18 = t78 * t49;
t7 = t17 * t48 + t18 * t49;
t67 = t7 * qJD(1);
t16 = -t51 - t73;
t8 = t16 * t49 - t20 * t48;
t66 = t8 * qJD(1);
t9 = -t16 * t48 - t20 * t49;
t65 = t9 * qJD(1);
t10 = t68 + t76 + (-pkin(4) / 0.2e1 - t75 / 0.2e1) * t48;
t64 = t10 * qJD(1);
t63 = t18 * qJD(3);
t45 = t48 ^ 2;
t46 = t49 ^ 2;
t24 = t45 + t46;
t61 = t24 * qJD(1);
t25 = t46 - t45;
t60 = t25 * qJD(1);
t22 = t25 * qJD(3);
t38 = t48 * qJD(1);
t37 = t48 * qJD(3);
t59 = t49 * qJD(1);
t39 = t49 * qJD(3);
t57 = t16 * t20 * qJD(1);
t56 = t16 * t38;
t55 = t31 * t38;
t54 = t31 * t59;
t53 = t30 * t37;
t52 = t30 * t39;
t44 = qJ(4) * qJD(4);
t43 = qJD(3) * qJ(4);
t35 = t45 * qJD(1);
t34 = t45 * qJD(4);
t28 = t48 * t39;
t27 = t48 * t59;
t26 = t48 * t58;
t15 = t72 / 0.2e1 + t76 - t74 / 0.2e1;
t2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, t22, 0, -t28, 0, 0, t31 * t37, t31 * t39, 0, 0, t28, 0, -t22, 0, 0, -t28, -t9 * qJD(3) + t26, 0, -t8 * qJD(3) + t34, -t79 * t16, t28, -t22, 0, -t28, 0, 0, t4 * qJD(3) + t26, t3 * qJD(3) + t34, t24 * qJD(5), t1 * qJD(3) - t7 * qJD(5) + t13 * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, t60, t39, -t27, -t37, 0, -t52 + t55, t53 + t54, 0, 0, t27, t39, -t60, 0, t37, -t27, -t52 - t65, t77, -t53 - t66, t30 * t77 - t57, t27, -t60, -t39, -t27, -t37, 0, -t63 + t69, -t17 * qJD(3) + t70, (t41 + t80) * qJD(3) - t58, t71 + (-t17 * qJ(4) - t18 * t75) * qJD(3) + t18 * qJD(4) + t15 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, t39, t35, t52 - t56, 0, 0, 0, 0, 0, 0, t27, t35, -t39, t13 * t38 + t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, t15 * qJD(3) - t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, -t39, 0, 0, 0, 0, 0, 0, 0, 0, -t37, 0, t39, t79, 0, 0, 0, 0, 0, 0, -t37, t39, 0, (t68 - t72) * qJD(3) + t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t27, -t60, 0, t27, 0, 0, -t55, -t54, 0, 0, -t27, 0, t60, 0, 0, t27, t65, 0, t66, t57, -t27, t60, 0, t27, 0, 0, t48 * qJD(5) - t69, -t49 * qJD(5) - t70, 0, -t10 * qJD(5) - t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4), t44, 0, 0, 0, 0, 0, 0, 0, qJD(4), 0, t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3), t43, 0, 0, 0, 0, 0, 0, 0, qJD(3), 0, t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, -t59, 0, -t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t27, 0, -t35, t56, 0, 0, 0, 0, 0, 0, -t27, -t35, 0, (-qJD(1) * t13 - qJD(5)) * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3), -t43, 0, 0, 0, 0, 0, 0, 0, -qJD(3), 0, -t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, t39, -t61, t10 * qJD(3) + t36 + t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38, t59, 0, t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t2;
