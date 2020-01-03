% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
% 
% Output:
% tauc_reg [5x17]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRPR6_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR6_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR6_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR6_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:17:49
% EndTime: 2019-12-31 18:17:50
% DurationCPUTime: 0.20s
% Computational Cost: add. (290->60), mult. (614->75), div. (0->0), fcn. (304->6), ass. (0->50)
t27 = qJD(1) + qJD(3);
t35 = cos(qJ(3));
t63 = pkin(1) * sin(pkin(8));
t50 = qJD(1) * t63;
t24 = cos(pkin(8)) * pkin(1) + pkin(2);
t22 = t24 * qJD(1);
t33 = sin(qJ(3));
t60 = t33 * t22;
t10 = t35 * t50 + t60;
t55 = t27 * qJ(4);
t7 = t10 + t55;
t49 = qJD(3) * t63;
t46 = qJD(1) * t49;
t8 = qJD(3) * t60 + t35 * t46;
t65 = -t7 * t27 + t8;
t26 = t27 ^ 2;
t32 = sin(qJ(5));
t25 = t27 * qJD(4);
t53 = qJD(3) * t35;
t39 = -t22 * t53 + t33 * t46;
t5 = -t25 + t39;
t34 = cos(qJ(5));
t51 = t34 * qJD(5);
t64 = -t5 * t32 + t7 * t51;
t41 = t24 * t53 - t33 * t49;
t11 = -qJD(4) - t41;
t61 = t11 * t27;
t37 = qJD(5) ^ 2;
t59 = t37 * t32;
t58 = t37 * t34;
t57 = -t26 - t37;
t56 = t32 ^ 2 - t34 ^ 2;
t9 = -t35 * t22 + t33 * t50;
t54 = qJD(4) + t9;
t52 = t32 * qJD(5);
t40 = t33 * t24 + t35 * t63;
t12 = t40 * qJD(3);
t15 = qJ(4) + t40;
t47 = t15 * t27 + t12;
t45 = -t35 * t24 + t33 * t63 - pkin(3);
t44 = t10 * t27 - t8;
t43 = t12 * t27 + t8;
t42 = -(-pkin(7) + t45) * t37 - t61;
t38 = -t9 * t27 + t39;
t36 = -pkin(3) - pkin(7);
t19 = -0.2e1 * t27 * t32 * t51;
t13 = 0.2e1 * t56 * t27 * qJD(5);
t6 = -t27 * pkin(3) + t54;
t2 = t5 * t34;
t1 = [0, 0, 0, 0, 0, -t43, -t41 * t27 + t39, t43, -t5 - t61, -t7 * t11 + t6 * t12 - t5 * t15 + t8 * t45, t19, t13, -t59, -t58, 0, t42 * t32 + t47 * t51 + t64, -t2 + t42 * t34 + (-t47 - t7) * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t58, t59; 0, 0, 0, 0, 0, t44, t38, -t44, 0.2e1 * t25 - t38, -t8 * pkin(3) - t5 * qJ(4) - t6 * t10 + t54 * t7, t19, t13, -t59, -t58, 0, -t10 * t51 - t36 * t59 + (qJ(4) * t51 + t54 * t32) * t27 + t64, -t2 + (t54 * t27 - t36 * t37) * t34 + (t10 - t7 - t55) * t52; 0, 0, 0, 0, 0, 0, 0, 0, -t26, t65, 0, 0, 0, 0, 0, t57 * t32, t57 * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34 * t26 * t32, -t56 * t26, 0, 0, 0, t65 * t34, -t65 * t32;];
tauc_reg = t1;
