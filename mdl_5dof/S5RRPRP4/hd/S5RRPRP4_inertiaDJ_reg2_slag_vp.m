% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:53
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRPRP4_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP4_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP4_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP4_inertiaDJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:52:55
% EndTime: 2019-12-31 19:52:57
% DurationCPUTime: 0.54s
% Computational Cost: add. (205->74), mult. (513->94), div. (0->0), fcn. (274->4), ass. (0->53)
t45 = sin(qJ(2));
t64 = pkin(1) * qJD(2);
t36 = t45 * t64;
t44 = sin(qJ(4));
t42 = t44 ^ 2;
t46 = cos(qJ(4));
t43 = t46 ^ 2;
t77 = (t42 + t43) * t36;
t75 = -t44 * pkin(4) + t46 * qJ(5);
t74 = 0.2e1 * (t42 - t43) * qJD(4);
t49 = 2 * qJD(3);
t73 = 2 * qJD(5);
t22 = qJ(3) - t75;
t71 = t45 * pkin(1);
t14 = t22 + t71;
t47 = cos(qJ(2));
t58 = t47 * t64;
t37 = t46 * qJD(4);
t59 = qJD(4) * qJ(5);
t9 = pkin(4) * t37 - t46 * qJD(5) + t44 * t59 + qJD(3);
t2 = t9 + t58;
t72 = t14 * t37 + t2 * t44;
t70 = t22 * t37 + t9 * t44;
t29 = qJD(3) + t58;
t31 = qJ(3) + t71;
t69 = t31 * t29;
t54 = -t47 * pkin(1) - pkin(2);
t30 = -pkin(7) + t54;
t68 = t77 * t30;
t67 = t29 * t44 + t31 * t37;
t48 = -pkin(2) - pkin(7);
t66 = t77 * t48;
t60 = qJ(3) * qJD(4);
t65 = qJD(3) * t44 + t46 * t60;
t62 = t44 * qJD(4);
t61 = t44 * qJD(5);
t57 = t48 * t62;
t56 = t44 * t37;
t55 = t48 * t37;
t51 = t29 * qJ(3) + t31 * qJD(3);
t50 = t75 * qJD(4) + t61;
t41 = qJ(3) * t49;
t39 = qJD(3) * t46;
t28 = -0.2e1 * t56;
t27 = 0.2e1 * t56;
t20 = t29 * t46;
t12 = t22 * t62;
t11 = -pkin(4) * t62 + t46 * t59 + t61;
t7 = t14 * t62;
t6 = 0.2e1 * t77;
t5 = t30 * t37 + t44 * t36;
t4 = t30 * t62 - t46 * t36;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t36, -0.2e1 * t58, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t36, 0.2e1 * t29, 0.2e1 * t54 * t36 + 0.2e1 * t69, t28, t74, 0, t27, 0, 0, 0.2e1 * t67, -0.2e1 * t31 * t62 + 0.2e1 * t20, -t6, 0.2e1 * t68 + 0.2e1 * t69, t28, 0, -t74, 0, 0, t27, 0.2e1 * t72, -t6, -0.2e1 * t2 * t46 + 0.2e1 * t7, 0.2e1 * t14 * t2 + 0.2e1 * t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36, -t58, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, t49 + t58, -pkin(2) * t36 + t51, t28, t74, 0, t27, 0, 0, t65 + t67, t20 + t39 + (-qJ(3) - t31) * t62, -t77, t51 + t66, t28, 0, -t74, 0, 0, t27, t70 + t72, -t77, t12 + t7 + (-t2 - t9) * t46, t14 * t9 + t2 * t22 + t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, t41, t28, t74, 0, t27, 0, 0, 0.2e1 * t65, -0.2e1 * t44 * t60 + 0.2e1 * t39, 0, t41, t28, 0, -t74, 0, 0, t27, 0.2e1 * t70, 0, -0.2e1 * t9 * t46 + 0.2e1 * t12, 0.2e1 * t22 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t62, 0, -t37, 0, -t4, -t5, 0, 0, 0, -t62, 0, 0, t37, 0, -t4, -t11, t5, (pkin(4) * t46 + qJ(5) * t44) * t36 + t50 * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t62, 0, -t37, 0, -t57, -t55, 0, 0, 0, -t62, 0, 0, t37, 0, -t57, -t11, t55, t50 * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t62, -t37, 0, 0, 0, 0, 0, 0, 0, 0, -t62, 0, t37, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73, qJ(5) * t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t62, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t62, 0, t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t1;
