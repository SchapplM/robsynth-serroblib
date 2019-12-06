% Calculate minimal parameter regressor of coriolis matrix for
% S5PPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x14]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5PPRRP1_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP1_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP1_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP1_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:07:16
% EndTime: 2019-12-05 15:07:17
% DurationCPUTime: 0.25s
% Computational Cost: add. (262->43), mult. (659->72), div. (0->0), fcn. (628->6), ass. (0->46)
t40 = sin(qJ(4));
t37 = t40 ^ 2;
t68 = pkin(4) * t40;
t67 = cos(qJ(3));
t39 = sin(pkin(8));
t41 = sin(qJ(3));
t63 = cos(pkin(8));
t27 = t41 * t39 - t67 * t63;
t9 = t27 * t40;
t66 = -qJ(5) - pkin(6);
t42 = cos(qJ(4));
t38 = t42 ^ 2;
t33 = t37 + t38;
t47 = t41 * t63;
t48 = t67 * t39;
t28 = t48 + t47;
t3 = (0.1e1 - t33) * t28 * t27;
t65 = t3 * qJD(1);
t36 = -t42 * pkin(4) - pkin(3);
t4 = t36 * t68;
t64 = t4 * qJD(3);
t62 = qJD(3) * t40;
t61 = qJD(3) * t42;
t60 = t27 * qJD(3);
t59 = t28 * qJD(3);
t58 = t33 * qJD(3);
t34 = t38 - t37;
t57 = t34 * qJD(3);
t56 = t40 * qJD(4);
t55 = t42 * qJD(4);
t54 = pkin(3) * t62;
t53 = pkin(3) * t61;
t52 = pkin(4) * t56;
t51 = pkin(4) * t62;
t50 = t40 * t61;
t49 = t28 * t55;
t46 = (t38 / 0.2e1 + t37 / 0.2e1) * t28;
t30 = t66 * t42;
t45 = t30 * t42 + t37 * t66;
t43 = -t48 / 0.2e1 - t47 / 0.2e1;
t6 = t46 + t43;
t44 = t6 * qJD(1) - qJD(3) * t45;
t10 = t27 * t42;
t5 = t46 - t43;
t2 = pkin(4) * t9;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, -t59, t60, 0, 0, 0, 0, 0, t9 * qJD(4) - t42 * t59, t10 * qJD(4) + t40 * t59, -t33 * t60, t65 + (t45 * t27 + t28 * t36) * qJD(3) + t2 * qJD(4) + t5 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9 * qJD(3) - t49, t10 * qJD(3) + t28 * t56, 0, -pkin(4) * t49 + t2 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t56, -t55, 0, -t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6 * qJD(5) - t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, t40 * t55, t34 * qJD(4), 0, 0, 0, -pkin(3) * t56, -pkin(3) * t55, t33 * qJD(5), t4 * qJD(4) - qJD(5) * t45; 0, 0, 0, 0, 0, t50, t57, t55, -t56, 0, -pkin(6) * t55 - t54, pkin(6) * t56 - t53, -pkin(4) * t55, t30 * pkin(4) * qJD(4) + t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -t50, -t57, 0, 0, 0, t54, t53, 0, -qJD(5) * t68 - t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t58, -t44 + t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t1;
