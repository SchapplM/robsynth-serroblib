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
% cmat_reg [(5*%NQJ)%x16]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 14:48
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
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
% StartTime: 2021-01-15 14:48:19
% EndTime: 2021-01-15 14:48:20
% DurationCPUTime: 0.31s
% Computational Cost: add. (332->54), mult. (829->83), div. (0->0), fcn. (803->6), ass. (0->57)
t58 = cos(qJ(4));
t86 = t58 * pkin(4);
t85 = cos(qJ(3));
t55 = sin(pkin(8));
t57 = sin(qJ(3));
t80 = cos(pkin(8));
t41 = t57 * t55 - t85 * t80;
t56 = sin(qJ(4));
t21 = t41 * t56;
t51 = -pkin(3) - t86;
t84 = t51 * t56;
t83 = pkin(6) + qJ(5);
t53 = t56 ^ 2;
t54 = t58 ^ 2;
t48 = t53 + t54;
t62 = t57 * t80;
t63 = t85 * t55;
t42 = t63 + t62;
t3 = (0.1e1 - t48) * t42 * t41;
t82 = t3 * qJD(1);
t8 = pkin(4) * t84;
t81 = t8 * qJD(3);
t37 = t56 * t86 - t84;
t79 = t37 * qJD(3);
t40 = t53 * pkin(4) + t51 * t58;
t78 = t40 * qJD(3);
t77 = t41 * qJD(3);
t44 = t83 * t58;
t76 = t44 * qJD(4);
t75 = t48 * qJD(3);
t49 = t54 - t53;
t74 = t49 * qJD(3);
t73 = t56 * qJD(3);
t72 = t56 * qJD(4);
t71 = t56 * qJD(5);
t70 = t58 * qJD(3);
t52 = t58 * qJD(4);
t69 = pkin(3) * t73;
t68 = pkin(3) * t70;
t67 = pkin(4) * t72;
t66 = pkin(4) * t73;
t65 = t56 * t70;
t64 = t42 * t52;
t61 = (t54 / 0.2e1 + t53 / 0.2e1) * t42;
t43 = t83 * t56;
t28 = t43 * t56 + t44 * t58;
t59 = -t63 / 0.2e1 - t62 / 0.2e1;
t10 = t61 + t59;
t60 = t10 * qJD(1) + t28 * qJD(3);
t22 = t41 * t58;
t9 = t61 - t59;
t7 = t22 * qJD(3) + t42 * t72;
t6 = t21 * qJD(3) - t64;
t5 = t21 * qJD(4) - t42 * t70;
t4 = t22 * qJD(4) + t42 * t73;
t2 = pkin(4) * t21;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, -t42 * qJD(3), t77, 0, 0, 0, 0, 0, t5, t4, t5, t4, -t48 * t77, t82 + (-t28 * t41 + t42 * t51) * qJD(3) + t2 * qJD(4) + t9 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, t7, t6, t7, 0, -pkin(4) * t64 + t2 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t72, -t52, -t72, -t52, 0, -t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10 * qJD(5) - t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, t56 * t52, t49 * qJD(4), 0, 0, 0, -pkin(3) * t72, -pkin(3) * t52, -t37 * qJD(4), t40 * qJD(4), t48 * qJD(5), t8 * qJD(4) + t28 * qJD(5); 0, 0, 0, 0, 0, t65, t74, t52, -t72, 0, -pkin(6) * t52 - t69, pkin(6) * t72 - t68, -t76 - t79, t43 * qJD(4) + t78, -pkin(4) * t52, -pkin(4) * t76 + t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -t65, -t74, 0, 0, 0, t69, t68, -t71 + t79, -t58 * qJD(5) - t78, 0, -pkin(4) * t71 - t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t73, -t70, 0, -t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72, t52, -t75, -t60 + t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73, t70, 0, t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t1;
