% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5PPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
% 
% Output:
% tauc_reg [5x15]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PPRRR1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:12:51
% EndTime: 2019-12-05 15:12:52
% DurationCPUTime: 0.32s
% Computational Cost: add. (337->65), mult. (849->93), div. (0->0), fcn. (680->8), ass. (0->58)
t37 = sin(pkin(9));
t38 = cos(pkin(9));
t41 = sin(qJ(3));
t44 = cos(qJ(3));
t77 = -t41 * t37 + t44 * t38;
t22 = t77 * qJD(1);
t27 = t44 * t37 + t41 * t38;
t23 = t27 * qJD(1);
t40 = sin(qJ(4));
t43 = cos(qJ(4));
t65 = t43 * t23;
t12 = t40 * t22 + t65;
t34 = qJD(3) + qJD(4);
t45 = qJD(5) ^ 2;
t60 = qJD(4) * t40;
t76 = (pkin(3) * t60 - t12) * t34 + (t40 * pkin(3) + pkin(7)) * t45;
t19 = qJD(3) * pkin(3) + t22;
t24 = t77 * qJD(3);
t20 = qJD(1) * t24;
t25 = t27 * qJD(3);
t21 = qJD(1) * t25;
t55 = -t40 * t21 - t23 * t60;
t75 = -(qJD(4) * t19 + t20) * t43 - t55;
t11 = t40 * t19 + t65;
t56 = t40 * t20 + t43 * t21;
t3 = t11 * qJD(4) + t56;
t39 = sin(qJ(5));
t42 = cos(qJ(5));
t67 = t40 * t23;
t10 = t43 * t19 - t67;
t73 = t34 * pkin(4);
t8 = -t10 - t73;
t61 = qJD(5) * t8;
t74 = t3 * t39 + t42 * t61;
t72 = t43 * pkin(3);
t15 = t43 * t27 + t40 * t77;
t71 = (t15 * qJD(4) + t40 * t24 + t43 * t25) * t34;
t70 = t11 * t34;
t68 = t39 * t42;
t63 = t45 * t39;
t62 = t39 ^ 2 - t42 ^ 2;
t59 = 0.2e1 * qJD(5) * t34;
t58 = -t8 * t34 + t75;
t57 = -pkin(3) * t34 - t19;
t52 = pkin(7) * t45 - t70;
t51 = t15 * t45 + t71;
t50 = -t40 * t27 + t43 * t77;
t48 = qJD(5) * (t10 - t73);
t4 = t50 * qJD(4) + t43 * t24 - t40 * t25;
t47 = qJD(5) * (-t34 * t50 - t4);
t13 = t43 * t22 - t67;
t46 = qJD(5) * (-qJD(4) * t72 + (-pkin(4) - t72) * t34 + t13);
t33 = t34 ^ 2;
t32 = t45 * t42;
t28 = t59 * t68;
t18 = t62 * t59;
t6 = t39 * t61;
t1 = [0, 0, 0, -t25 * qJD(3), -t24 * qJD(3), 0, -t71, -t4 * t34, 0, 0, 0, 0, 0, t39 * t47 - t51 * t42, t51 * t39 + t42 * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t63, -t32; 0, 0, 0, 0, 0, 0, t12 * t34 + (t57 * t40 - t65) * qJD(4) - t56, t13 * t34 + (t57 * qJD(4) - t20) * t43 - t55, t28, -t18, t32, -t63, 0, t6 + t39 * t46 + (-t3 - t76) * t42, t76 * t39 + t42 * t46 + t74; 0, 0, 0, 0, 0, 0, -t3 + t70, t10 * t34 + t75, t28, -t18, t32, -t63, 0, t6 + t39 * t48 + (-t3 - t52) * t42, t52 * t39 + t42 * t48 + t74; 0, 0, 0, 0, 0, 0, 0, 0, -t33 * t68, t62 * t33, 0, 0, 0, t58 * t39, t58 * t42;];
tauc_reg = t1;
