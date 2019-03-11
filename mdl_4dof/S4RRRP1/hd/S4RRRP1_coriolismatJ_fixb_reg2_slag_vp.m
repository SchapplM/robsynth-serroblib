% Calculate inertial parameters regressor of coriolis matrix for
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
% cmat_reg [(4*4)x(4*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4RRRP1_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP1_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP1_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP1_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:36:05
% EndTime: 2019-03-08 18:36:06
% DurationCPUTime: 0.33s
% Computational Cost: add. (429->63), mult. (1024->97), div. (0->0), fcn. (721->4), ass. (0->70)
t53 = sin(qJ(2));
t54 = cos(qJ(3));
t78 = t54 * t53;
t52 = sin(qJ(3));
t55 = cos(qJ(2));
t79 = t52 * t55;
t42 = (t78 + t79) * pkin(1);
t84 = t42 * pkin(3);
t27 = t52 * pkin(2);
t83 = t54 * pkin(2);
t82 = t55 * pkin(1);
t61 = pkin(2) + t82;
t58 = t52 * t61;
t39 = pkin(1) * t78 + t58;
t77 = t54 * t55;
t80 = t52 * t53;
t43 = (t77 - t80) * pkin(1);
t81 = t39 * t43;
t76 = pkin(1) * qJD(1);
t75 = pkin(1) * qJD(2);
t74 = pkin(2) * qJD(2);
t73 = pkin(2) * qJD(3);
t46 = t54 * t61;
t38 = pkin(1) * t80 - t46;
t35 = pkin(3) - t38;
t5 = (-t35 - t38) * t39;
t72 = t5 * qJD(1);
t9 = -t35 * t42 + t81;
t71 = t9 * qJD(1);
t10 = t38 * t42 + t81;
t70 = t10 * qJD(1);
t69 = t27 * qJD(1);
t28 = t46 / 0.2e1 + (-t82 / 0.2e1 + pkin(2) / 0.2e1) * t54;
t68 = t28 * qJD(1);
t67 = t38 * qJD(1);
t66 = t39 * qJD(1);
t33 = t39 * qJD(3);
t65 = t42 * qJD(1);
t64 = t43 * qJD(1);
t63 = t52 * t73;
t62 = t54 * t73;
t60 = pkin(1) * (-qJD(1) - qJD(2));
t59 = pkin(2) * (-qJD(2) - qJD(3));
t49 = pkin(3) + t83;
t56 = (t39 * t54 / 0.2e1 + (-t38 / 0.2e1 - t35 / 0.2e1) * t52) * pkin(2) - t39 * t49 / 0.2e1;
t2 = t84 / 0.2e1 + t56;
t36 = (-t49 + t83) * t27;
t57 = -t2 * qJD(1) - t36 * qJD(2);
t41 = t43 * qJD(2);
t40 = t42 * qJD(2);
t34 = t43 * t27;
t32 = t38 * qJD(3);
t24 = -t83 / 0.2e1 - t46 / 0.2e1 + (t80 - t77 / 0.2e1) * pkin(1);
t23 = -t27 / 0.2e1 - t58 / 0.2e1 + (-t78 - t79 / 0.2e1) * pkin(1);
t22 = t54 * t74 + t68;
t21 = t52 * t74 + t69;
t20 = -t41 + t32;
t19 = -t40 - t33;
t18 = t54 * t59 - t68;
t17 = t52 * t59 - t69;
t16 = -t28 * qJD(3) + t64;
t15 = -t27 * qJD(3) + t65;
t13 = t27 * qJD(2) + t66;
t12 = t28 * qJD(2) - t67;
t8 = t24 * qJD(3) - t41 - t64;
t7 = t23 * qJD(3) - t40 - t65;
t4 = t23 * qJD(2) - t33 - t66;
t3 = t24 * qJD(2) + t32 + t67;
t1 = -t84 / 0.2e1 + t56;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t53 * t75, -t55 * t75, 0, 0, 0, 0, 0, 0, 0, 0, t19, t20, 0, t10 * qJD(2), 0, 0, 0, 0, 0, 0, t19, t20, 0, t9 * qJD(2) + t5 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53 * t60, t55 * t60, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, t70 + (-t42 * t83 + t34) * qJD(2), 0, 0, 0, 0, 0, 0, t7, t8, 0, t71 + (-t42 * t49 + t34) * qJD(2) + t1 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t3, 0, 0, 0, 0, 0, 0, 0, 0, t4, t3, 0, -pkin(3) * t33 + t1 * qJD(2) + t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53 * t76, t55 * t76, 0, 0, 0, 0, 0, 0, 0, 0, t15, t16, 0, -t70, 0, 0, 0, 0, 0, 0, t15, t16, 0, t2 * qJD(3) - t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t63, -t62, 0, 0, 0, 0, 0, 0, 0, 0, -t63, -t62, 0, t36 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, t18, 0, 0, 0, 0, 0, 0, 0, 0, t17, t18, 0, -pkin(3) * t63 - t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, t12, 0, 0, 0, 0, 0, 0, 0, 0, t13, t12, 0, -t2 * qJD(2) - t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, t22, 0, 0, 0, 0, 0, 0, 0, 0, t21, t22, 0, t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t6;
