% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RPRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x17]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:09
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRPP1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP1_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP1_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPP1_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:09:13
% EndTime: 2019-12-31 18:09:14
% DurationCPUTime: 0.20s
% Computational Cost: add. (353->49), mult. (770->96), div. (0->0), fcn. (649->6), ass. (0->38)
t51 = 2 * qJD(5);
t30 = sin(pkin(8));
t32 = sin(qJ(3));
t50 = t30 * t32;
t26 = sin(pkin(7)) * pkin(1) + pkin(6);
t49 = qJ(4) + t26;
t48 = cos(pkin(8));
t33 = cos(qJ(3));
t42 = t48 * t32;
t22 = t30 * t33 + t42;
t47 = t22 * qJD(5);
t46 = t32 * qJD(3);
t45 = t33 * qJD(3);
t44 = 0.2e1 * t45;
t29 = pkin(3) * t46;
t28 = -cos(pkin(7)) * pkin(1) - pkin(2);
t20 = t49 * t33;
t10 = t30 * t20 + t49 * t42;
t11 = t48 * t20 - t49 * t50;
t40 = qJD(3) * t49;
t16 = t33 * qJD(4) - t32 * t40;
t35 = -t32 * qJD(4) - t33 * t40;
t7 = t30 * t16 - t48 * t35;
t8 = t48 * t16 + t30 * t35;
t43 = t10 * t7 + t11 * t8;
t41 = t48 * t33;
t39 = -t33 * pkin(3) + t28;
t18 = t22 * qJD(3);
t19 = qJD(3) * t41 - t30 * t46;
t21 = -t41 + t50;
t38 = 0.2e1 * t21 * t18 + 0.2e1 * t22 * t19;
t37 = t10 * t18 + t11 * t19 + t7 * t21 + t8 * t22;
t34 = 0.2e1 * t10 * t19 - 0.2e1 * t11 * t18 - 0.2e1 * t8 * t21 + 0.2e1 * t7 * t22;
t27 = -t48 * pkin(3) - pkin(4);
t24 = t30 * pkin(3) + qJ(5);
t9 = t21 * pkin(4) - t22 * qJ(5) + t39;
t4 = t18 * pkin(4) - t19 * qJ(5) + t29 - t47;
t1 = [0, 0, 0, 0, t32 * t44, 0.2e1 * (-t32 ^ 2 + t33 ^ 2) * qJD(3), 0, 0, 0, 0.2e1 * t28 * t46, t28 * t44, t34, 0.2e1 * t39 * t29 + 0.2e1 * t43, 0.2e1 * t9 * t18 + 0.2e1 * t4 * t21, t34, -0.2e1 * t9 * t19 - 0.2e1 * t4 * t22, 0.2e1 * t9 * t4 + 0.2e1 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, 0, 0, 0, t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, 0, 0, 0, t38; 0, 0, 0, 0, 0, 0, t45, -t46, 0, -t26 * t45, t26 * t46, (-t18 * t30 - t48 * t19) * pkin(3), (t30 * t8 - t48 * t7) * pkin(3), -t7, -qJD(5) * t21 - t24 * t18 + t27 * t19, t8, t11 * qJD(5) + t8 * t24 + t7 * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46, -t45, 0, (-t48 * t18 + t19 * t30) * pkin(3), -t18, 0, t19, t18 * t27 + t19 * t24 + t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, t24 * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, t18, 0, -t19, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t1;
