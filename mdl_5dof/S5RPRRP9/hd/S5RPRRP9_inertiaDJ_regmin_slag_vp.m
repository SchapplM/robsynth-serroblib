% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RPRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x25]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRRP9_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP9_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP9_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP9_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:49:11
% EndTime: 2019-12-31 18:49:13
% DurationCPUTime: 0.44s
% Computational Cost: add. (917->84), mult. (2065->139), div. (0->0), fcn. (2005->6), ass. (0->50)
t33 = sin(pkin(8));
t34 = cos(pkin(8));
t66 = sin(qJ(3));
t68 = cos(qJ(3));
t22 = t68 * t33 + t66 * t34;
t65 = pkin(6) + qJ(2);
t72 = t68 * t65;
t70 = t66 * qJD(2) + qJD(3) * t72;
t73 = t65 * t66;
t75 = -t68 * qJD(2) + qJD(3) * t73;
t74 = t33 * t75 - t70 * t34;
t69 = t70 * t33 + t75 * t34;
t36 = 2 * qJD(5);
t67 = cos(qJ(4));
t64 = t22 * qJD(3);
t35 = sin(qJ(4));
t63 = qJD(4) * t35;
t62 = pkin(3) * t63;
t28 = -t34 * pkin(2) - pkin(1);
t61 = t64 * pkin(3);
t58 = qJD(4) * t67;
t55 = 0.2e1 * (t33 ^ 2 + t34 ^ 2) * qJD(2);
t48 = t66 * t33 - t68 * t34;
t47 = t35 * t48;
t46 = t48 * qJD(3);
t45 = -0.2e1 * t46;
t44 = t67 * t48;
t16 = t48 * pkin(3) + t28;
t41 = -t22 * pkin(7) - t33 * t72 - t34 * t73;
t40 = t35 * t41;
t39 = t67 * t41;
t38 = -t64 * pkin(7) - t69;
t37 = -pkin(7) * t46 - t74;
t30 = pkin(3) * t58;
t29 = -t67 * pkin(3) - pkin(4);
t27 = pkin(3) * t35 + qJ(5);
t26 = -0.2e1 * t62;
t23 = t30 + qJD(5);
t15 = t67 * t22 - t47;
t14 = t35 * t22 + t44;
t13 = (t68 * pkin(7) + t72) * t34 + (-t66 * pkin(7) - t73) * t33;
t8 = -qJD(4) * t47 + t22 * t58 - t35 * t46 + t67 * t64;
t7 = t22 * t63 + t35 * t64 - (-qJD(3) - qJD(4)) * t44;
t6 = t67 * t13 + t40;
t5 = t35 * t13 - t39;
t4 = t14 * pkin(4) - t15 * qJ(5) + t16;
t3 = t8 * pkin(4) + t7 * qJ(5) - t15 * qJD(5) + t61;
t2 = qJD(4) * t40 + t13 * t58 + t35 * t38 + t67 * t37;
t1 = -qJD(4) * t39 + t13 * t63 + t35 * t37 - t67 * t38;
t9 = [0, 0, 0, 0, 0, t55, qJ(2) * t55, t22 * t45, 0.2e1 * t48 ^ 2 * qJD(3) - 0.2e1 * t22 * t64, 0, 0, 0, 0.2e1 * t28 * t64, t28 * t45, -0.2e1 * t15 * t7, 0.2e1 * t14 * t7 - 0.2e1 * t15 * t8, 0, 0, 0, 0.2e1 * t14 * t61 + 0.2e1 * t16 * t8, 0.2e1 * t15 * t61 - 0.2e1 * t16 * t7, 0.2e1 * t14 * t3 + 0.2e1 * t4 * t8, 0.2e1 * t1 * t14 + 0.2e1 * t15 * t2 - 0.2e1 * t5 * t7 - 0.2e1 * t6 * t8, -0.2e1 * t15 * t3 + 0.2e1 * t4 * t7, -0.2e1 * t1 * t6 + 0.2e1 * t2 * t5 + 0.2e1 * t3 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64, -t46, 0, 0, 0, 0, 0, t8, -t7, t8, 0, t7, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46, -t64, 0, t74, t69, 0, 0, -t7, -t8, 0, -t2, t1, -t2, -t14 * t23 + t15 * t62 - t27 * t8 - t29 * t7, -t1, -t1 * t27 + t2 * t29 + t23 * t6 + t5 * t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, -0.2e1 * t30, t26, 0, 0.2e1 * t23, 0.2e1 * t23 * t27 + 0.2e1 * t29 * t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, -t8, 0, -t2, t1, -t2, pkin(4) * t7 - qJ(5) * t8 - qJD(5) * t14, -t1, -pkin(4) * t2 - qJ(5) * t1 + qJD(5) * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t62, -t30, -t62, 0, t36 + t30, -pkin(4) * t62 + qJ(5) * t23 + qJD(5) * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, qJ(5) * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t9;
