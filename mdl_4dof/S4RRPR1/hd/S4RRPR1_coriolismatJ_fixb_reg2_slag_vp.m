% Calculate inertial parameters regressor of coriolis matrix for
% S4RRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% cmat_reg [(4*4)x(4*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:54
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function cmat_reg = S4RRPR1_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR1_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR1_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR1_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:53:31
% EndTime: 2018-11-14 13:53:32
% DurationCPUTime: 0.36s
% Computational Cost: add. (691->65), mult. (1450->94), div. (0->0), fcn. (1230->6), ass. (0->65)
t50 = sin(pkin(7));
t84 = t50 * pkin(2);
t51 = cos(pkin(7));
t83 = t51 * pkin(2);
t54 = cos(qJ(2));
t82 = t54 * pkin(1);
t81 = cos(qJ(4));
t53 = sin(qJ(2));
t80 = t51 * t53;
t40 = (t50 * t54 + t80) * pkin(1);
t52 = sin(qJ(4));
t79 = t52 * t40;
t48 = t50 * t53 * pkin(1);
t41 = t51 * t82 - t48;
t78 = t52 * t41;
t49 = pkin(2) + t82;
t44 = t51 * t49;
t77 = -t48 + t44;
t76 = pkin(1) * qJD(1);
t75 = pkin(1) * qJD(2);
t66 = pkin(3) + t77;
t27 = t81 * t66;
t36 = pkin(1) * t80 + t50 * t49;
t17 = t52 * t36 - t27;
t63 = t81 * t36;
t18 = t52 * t66 + t63;
t34 = t81 * t40;
t21 = t34 + t78;
t62 = t81 * t41;
t22 = t62 - t79;
t4 = t17 * t21 + t18 * t22;
t74 = t4 * qJD(1);
t10 = t36 * t41 - t77 * t40;
t73 = t10 * qJD(1);
t72 = t17 * qJD(1);
t71 = t18 * qJD(1);
t70 = t21 * qJD(1);
t69 = t22 * qJD(1);
t68 = t40 * qJD(1);
t67 = t41 * qJD(1);
t64 = pkin(3) + t83;
t45 = t81 * t64;
t65 = -t27 / 0.2e1 - t45 / 0.2e1;
t61 = pkin(1) * (-qJD(1) - qJD(2));
t60 = t81 * t84;
t59 = t84 / 0.2e1 + t36 / 0.2e1;
t35 = t52 * t84 - t45;
t5 = t62 / 0.2e1 + (-t40 / 0.2e1 + t59) * t52 + t65;
t58 = -t5 * qJD(1) - t35 * qJD(2);
t37 = t52 * t64 + t60;
t55 = -t63 / 0.2e1 - t60 / 0.2e1;
t56 = -t83 / 0.2e1 - pkin(3) + t48 / 0.2e1 - t44 / 0.2e1;
t7 = t34 / 0.2e1 + (t41 / 0.2e1 + t56) * t52 + t55;
t57 = -t7 * qJD(1) + t37 * qJD(2);
t39 = t41 * qJD(2);
t38 = t40 * qJD(2);
t33 = t37 * qJD(4);
t32 = t35 * qJD(4);
t20 = t22 * qJD(2);
t19 = t21 * qJD(2);
t16 = t18 * qJD(4);
t15 = t17 * qJD(4);
t8 = -t78 / 0.2e1 - t34 / 0.2e1 + t56 * t52 + t55;
t6 = -t62 / 0.2e1 + t79 / 0.2e1 + t59 * t52 + t65;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t53 * t75, -t54 * t75, 0, 0, 0, 0, 0, 0, 0, 0, -t38, -t39, 0, t10 * qJD(2), 0, 0, 0, 0, 0, 0, -t19 - t16, -t20 + t15, 0, t4 * qJD(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53 * t61, t54 * t61, 0, 0, 0, 0, 0, 0, 0, 0, -t38 - t68, -t39 - t67, 0, t73 + (-t40 * t51 + t41 * t50) * qJD(2) * pkin(2), 0, 0, 0, 0, 0, 0, t8 * qJD(4) - t19 - t70, t6 * qJD(4) - t20 - t69, 0, t74 + (t21 * t35 + t22 * t37) * qJD(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8 * qJD(2) - t16 - t71, t6 * qJD(2) + t15 + t72, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53 * t76, t54 * t76, 0, 0, 0, 0, 0, 0, 0, 0, t68, t67, 0, -t73, 0, 0, 0, 0, 0, 0, t7 * qJD(4) + t70, t5 * qJD(4) + t69, 0, -t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33, t32, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33 - t57, t32 - t58, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7 * qJD(2) + t71, -t5 * qJD(2) - t72, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, t58, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t1;
