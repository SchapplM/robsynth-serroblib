% Calculate minimal parameter regressor of coriolis matrix for
% S5PPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x15]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5PPRRR5_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR5_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR5_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRR5_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:35:46
% EndTime: 2019-12-31 17:35:47
% DurationCPUTime: 0.28s
% Computational Cost: add. (210->46), mult. (536->74), div. (0->0), fcn. (486->6), ass. (0->49)
t46 = cos(qJ(4));
t66 = t46 * pkin(3);
t38 = -pkin(4) - t66;
t71 = pkin(4) / 0.2e1 - t38 / 0.2e1;
t42 = sin(qJ(5));
t45 = cos(qJ(5));
t34 = -t42 ^ 2 + t45 ^ 2;
t60 = qJD(3) + qJD(4);
t70 = t60 * t34;
t57 = -t66 / 0.2e1;
t69 = t57 - t71;
t65 = pkin(3) * qJD(4);
t64 = pkin(4) * qJD(4);
t63 = qJD(3) * pkin(3);
t62 = qJD(3) * t38;
t61 = qJD(5) * t42;
t41 = qJD(5) * t45;
t43 = sin(qJ(4));
t59 = t43 * t65;
t58 = t43 * t63;
t56 = t42 * t62;
t55 = t45 * t62;
t54 = pkin(3) * t60;
t53 = t42 * t58;
t44 = sin(qJ(3));
t47 = cos(qJ(3));
t27 = t43 * t47 + t46 * t44;
t14 = t60 * t27;
t52 = t60 * t42;
t51 = t43 * t54;
t50 = t57 + t71;
t15 = t50 * t42;
t49 = t15 * qJD(3) + t42 * t64;
t16 = t50 * t45;
t48 = t16 * qJD(3) + t45 * t64;
t37 = t43 * pkin(3) + pkin(7);
t35 = t42 * t41;
t33 = t42 * t59;
t30 = t34 * qJD(5);
t26 = t43 * t44 - t46 * t47;
t25 = t45 * t52;
t18 = t69 * t45;
t17 = t69 * t42;
t13 = t60 * t26;
t12 = t26 * t45;
t11 = t26 * t42;
t2 = t11 * qJD(5) - t45 * t14;
t1 = t12 * qJD(5) + t27 * t52;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, -qJD(3) * t44, -qJD(3) * t47, 0, -t14, t13, 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, 0, -t14, t13, 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60 * t11 - t27 * t41, t60 * t12 + t27 * t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, -t59, -t46 * t65, t35, t30, 0, 0, 0, t38 * t61 - t45 * t59, t38 * t41 + t33; 0, 0, 0, 0, 0, 0, -t51, -t46 * t54, t35, t30, 0, 0, 0, t17 * qJD(5) - t45 * t51, t18 * qJD(5) + t33 + t53; 0, 0, 0, 0, 0, 0, 0, 0, t25, t70, t41, -t61, 0, t17 * qJD(4) - t37 * t41 + t56, t18 * qJD(4) + t37 * t61 + t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t58, t46 * t63, t35, t30, 0, 0, 0, -t15 * qJD(5) + t45 * t58, -t16 * qJD(5) - t53; 0, 0, 0, 0, 0, 0, 0, 0, t35, t30, 0, 0, 0, -pkin(4) * t61, -pkin(4) * t41; 0, 0, 0, 0, 0, 0, 0, 0, t25, t70, t41, -t61, 0, -pkin(7) * t41 - t49, pkin(7) * t61 - t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, -t25, -t70, 0, 0, 0, t15 * qJD(4) - t56, t16 * qJD(4) - t55; 0, 0, 0, 0, 0, 0, 0, 0, -t25, -t70, 0, 0, 0, t49, t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t3;
