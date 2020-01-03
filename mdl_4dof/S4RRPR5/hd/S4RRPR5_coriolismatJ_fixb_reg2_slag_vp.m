% Calculate inertial parameters regressor of coriolis matrix for
% S4RRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
% 
% Output:
% cmat_reg [(4*4)x(4*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4RRPR5_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR5_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR5_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR5_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:03:36
% EndTime: 2019-12-31 17:03:37
% DurationCPUTime: 0.51s
% Computational Cost: add. (214->91), mult. (517->89), div. (0->0), fcn. (315->4), ass. (0->68)
t48 = sin(qJ(2));
t82 = t48 * pkin(1);
t34 = qJ(3) + t82;
t86 = -qJ(3) / 0.2e1 - t34 / 0.2e1;
t47 = sin(qJ(4));
t45 = t47 ^ 2;
t49 = cos(qJ(4));
t46 = t49 ^ 2;
t31 = t45 - t46;
t43 = qJD(1) + qJD(2);
t85 = t43 * t31;
t50 = cos(qJ(2));
t81 = t50 * pkin(1);
t79 = t34 * t50;
t75 = pkin(1) * qJD(2);
t39 = t50 * t75;
t41 = qJD(3) * t47;
t78 = t47 * t39 + t41;
t42 = qJD(3) * t49;
t77 = t49 * t39 + t42;
t76 = pkin(1) * qJD(1);
t59 = -pkin(2) - t81;
t33 = -pkin(6) + t59;
t58 = (t45 + t46) * t48;
t1 = (t33 * t58 + t79) * pkin(1);
t74 = t1 * qJD(1);
t8 = -pkin(1) * t79 - t59 * t82;
t73 = t8 * qJD(1);
t13 = pkin(1) * t58;
t72 = t13 * qJD(1);
t71 = t34 * qJD(1);
t70 = t47 * qJD(4);
t69 = t49 * qJD(4);
t68 = t39 + qJD(3);
t67 = qJ(3) * qJD(2);
t66 = qJ(3) * qJD(4);
t65 = t48 * t75;
t64 = t48 * t76;
t38 = t50 * t76;
t63 = t82 / 0.2e1;
t62 = t47 * t71;
t61 = t49 * t71;
t60 = t45 / 0.2e1 + t46 / 0.2e1;
t57 = t47 * t38;
t56 = t49 * t38;
t18 = t43 * t49;
t6 = -qJ(3) + (-0.1e1 / 0.2e1 + t60) * t82;
t55 = -t6 * qJD(1) + t67;
t54 = t63 + t86;
t2 = t54 * t47;
t53 = -t2 * qJD(1) + t47 * t67;
t3 = t54 * t49;
t52 = t3 * qJD(1) - t49 * t67;
t51 = -pkin(2) - pkin(6);
t44 = qJ(3) * qJD(3);
t40 = qJ(3) * t81;
t32 = t47 * t69;
t30 = t43 * qJ(3);
t27 = t34 * qJD(3);
t20 = t31 * qJD(4);
t17 = t43 * t47;
t14 = t43 * t82;
t12 = t47 * t18;
t11 = t13 * qJD(2);
t7 = t60 * t82 + qJ(3) + t63;
t5 = (t63 - t86) * t49;
t4 = (-t82 / 0.2e1 + t86) * t47;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t65, -t39, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, t68, -t8 * qJD(2) + t27, -t32, t20, 0, t32, 0, 0, t34 * t69 + t78, -t34 * t70 + t77, -t11, t1 * qJD(2) + t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, -t38 - t39, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, t38 + t68, -t73 + (-pkin(2) * t82 + t40) * qJD(2) + t27, -t32, t20, 0, t32, 0, 0, t5 * qJD(4) + t57 + t78, t4 * qJD(4) + t56 + t77, -t11 - t72, t74 + (t51 * t13 + t40) * qJD(2) + t7 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, t34 * qJD(2) + t71, 0, 0, 0, 0, 0, 0, t17, t18, 0, t7 * qJD(2) + t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, t85, -t70, t12, -t69, 0, t5 * qJD(2) - t33 * t70 + t61, t4 * qJD(2) - t33 * t69 - t62, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64, t38, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t64, -t38 + qJD(3), t44 + t73, -t32, t20, 0, t32, 0, 0, -t3 * qJD(4) + t41 - t57, t2 * qJD(4) + t42 - t56, t72, -t6 * qJD(3) - t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3), t44, -t32, t20, 0, t32, 0, 0, t49 * t66 + t41, -t47 * t66 + t42, 0, t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, t30, 0, 0, 0, 0, 0, 0, t17, t18, 0, t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, t85, -t70, t12, -t69, 0, -t51 * t70 - t52, -t51 * t69 - t53, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43, -t71 - t67, 0, 0, 0, 0, 0, 0, -t17, -t18, 0, t6 * qJD(2) - t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43, -t30, 0, 0, 0, 0, 0, 0, -t17, -t18, 0, -t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t70, -t69, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, -t85, 0, -t12, 0, 0, t3 * qJD(2) - t61, -t2 * qJD(2) + t62, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, -t85, 0, -t12, 0, 0, t52, t53, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t9;
