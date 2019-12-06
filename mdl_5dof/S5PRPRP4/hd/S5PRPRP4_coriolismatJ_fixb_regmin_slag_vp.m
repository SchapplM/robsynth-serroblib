% Calculate minimal parameter regressor of coriolis matrix for
% S5PRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x16]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5PRPRP4_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP4_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP4_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP4_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:36:09
% EndTime: 2019-12-05 15:36:10
% DurationCPUTime: 0.37s
% Computational Cost: add. (317->56), mult. (726->86), div. (0->0), fcn. (755->6), ass. (0->53)
t42 = sin(qJ(4));
t43 = cos(qJ(4));
t50 = t43 * pkin(4) + t42 * qJ(5);
t59 = t43 * qJD(5);
t75 = t50 * qJD(4) - t59;
t74 = t42 * pkin(4);
t73 = cos(qJ(2));
t72 = sin(qJ(2));
t40 = t42 ^ 2;
t41 = t43 ^ 2;
t71 = -t40 - t41;
t67 = sin(pkin(8));
t68 = cos(pkin(8));
t28 = t67 * t72 - t68 * t73;
t29 = -t67 * t73 - t68 * t72;
t1 = (-0.1e1 - t71) * t29 * t28;
t70 = t1 * qJD(1);
t69 = t43 * qJ(5);
t66 = qJD(2) * t42;
t65 = qJD(2) * t43;
t64 = qJD(4) * t42;
t38 = qJD(4) * t43;
t37 = -t68 * pkin(2) - pkin(3);
t27 = -t50 + t37;
t30 = t69 - t74;
t17 = t27 * t43 - t30 * t42;
t63 = t17 * qJD(2);
t18 = -t27 * t42 - t30 * t43;
t62 = t18 * qJD(2);
t31 = t41 - t40;
t61 = t31 * qJD(2);
t60 = t40 * qJD(2);
t58 = qJD(4) * qJ(5);
t36 = t67 * pkin(2) + pkin(6);
t57 = t36 * t64;
t56 = t36 * t38;
t55 = t27 * t66;
t54 = t37 * t66;
t53 = t37 * t65;
t52 = t71 * t28;
t49 = t30 * qJD(4) + t42 * qJD(5);
t48 = -t69 / 0.2e1 + t74 / 0.2e1;
t2 = (t30 / 0.2e1 + t48) * t28;
t47 = t27 * t30 * qJD(2) + t2 * qJD(1);
t14 = t28 * t42;
t5 = t14 * qJD(2) + t29 * t38;
t16 = t28 * t43;
t46 = -t16 * qJD(2) + t29 * t64;
t45 = -t16 * qJD(4) + t29 * t66;
t32 = t42 * t65;
t4 = t14 * qJD(4) + t29 * t65;
t3 = (-t30 / 0.2e1 + t48) * t28;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1 * qJD(2); 0, 0, -qJD(2) * t72, -qJD(2) * t73, (-t67 * t28 + t68 * t29) * qJD(2) * pkin(2), 0, 0, 0, 0, 0, t4, -t45, t4, qJD(2) * t52, t45, t70 + (-t29 * t27 + t36 * t52) * qJD(2) + t3 * qJD(4) - t14 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, -t46, t5, 0, t46, t3 * qJD(2) + t75 * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2 * qJD(4) - t70; 0, 0, 0, 0, 0, t42 * t38, t31 * qJD(4), 0, 0, 0, t37 * t64, t37 * t38, -t18 * qJD(4) + t42 * t59, 0, -t17 * qJD(4) + t40 * qJD(5), -t49 * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, t32, t61, t38, -t64, 0, -t56 + t54, t53 + t57, -t56 - t62, -t75, -t57 - t63, -t36 * t75 - t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, t38, t60, -t55 + t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t64, -t38, -t64, 0, t38, t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2 * qJD(2); 0, 0, 0, 0, 0, -t32, -t61, 0, 0, 0, -t54, -t53, t62, 0, t63, t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(5), qJ(5) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4), t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32, 0, -t60, t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(4), -t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t6;
