% Calculate minimal parameter regressor of coriolis matrix for
% S5RPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x17]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:19
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RPRPR2_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR2_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR2_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR2_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:19:09
% EndTime: 2022-01-23 09:19:10
% DurationCPUTime: 0.47s
% Computational Cost: add. (619->95), mult. (1270->118), div. (0->0), fcn. (1168->8), ass. (0->82)
t62 = sin(pkin(9));
t60 = t62 ^ 2;
t64 = cos(pkin(9));
t61 = t64 ^ 2;
t54 = t60 + t61;
t48 = t54 * qJ(4);
t66 = sin(qJ(3));
t75 = cos(pkin(8)) * pkin(1) + pkin(2);
t72 = t66 * t75;
t101 = cos(qJ(3));
t104 = pkin(1) * sin(pkin(8));
t76 = t101 * t104;
t38 = t76 + t72;
t36 = qJ(4) + t38;
t11 = t54 * t36;
t100 = cos(qJ(5));
t79 = t100 * t62;
t65 = sin(qJ(5));
t98 = t65 * t64;
t46 = t79 + t98;
t103 = t64 * pkin(4);
t37 = -t101 * t75 + t66 * t104;
t74 = -pkin(3) + t37;
t29 = t74 - t103;
t58 = -pkin(3) - t103;
t81 = t58 / 0.2e1 + t29 / 0.2e1;
t112 = t81 * t46;
t78 = t100 * t64;
t99 = t65 * t62;
t44 = -t78 + t99;
t15 = t44 ^ 2 - t46 ^ 2;
t86 = qJD(1) + qJD(3);
t111 = t86 * t15;
t110 = t86 * t44;
t77 = t86 * t46;
t109 = t86 * t54;
t14 = t54 * t37;
t51 = t54 * qJD(4);
t97 = -t14 * qJD(3) + t51;
t1 = -t36 * t14 + t74 * t38;
t96 = t1 * qJD(1);
t95 = qJD(1) * t11;
t94 = t14 * qJD(1);
t93 = t37 * qJD(1);
t92 = t38 * qJD(1);
t35 = t38 * qJD(3);
t91 = t44 * qJD(1);
t90 = t44 * qJD(3);
t40 = t44 * qJD(5);
t89 = t46 * qJD(1);
t88 = t46 * qJD(3);
t87 = t46 * qJD(5);
t85 = t29 * t91;
t84 = t29 * t89;
t83 = t38 * t91;
t82 = t38 * t89;
t67 = t76 / 0.2e1 + t72 / 0.2e1;
t8 = t67 + (qJ(4) + t36) * (-t61 / 0.2e1 - t60 / 0.2e1);
t73 = t8 * qJD(1) - t48 * qJD(3);
t69 = (t98 / 0.2e1 + t79 / 0.2e1) * t37;
t2 = t69 - t112;
t71 = -t2 * qJD(1) + t58 * t88;
t68 = (t78 / 0.2e1 - t99 / 0.2e1) * t37;
t3 = t81 * t44 + t68;
t70 = -t3 * qJD(1) - t58 * t90;
t59 = t64 * pkin(7);
t50 = t64 * qJ(4) + t59;
t49 = (-pkin(7) - qJ(4)) * t62;
t43 = t46 * qJD(4);
t39 = t44 * qJD(4);
t34 = t37 * qJD(3);
t22 = t64 * t36 + t59;
t21 = (-pkin(7) - t36) * t62;
t20 = t44 * t87;
t19 = t38 * t88;
t18 = t38 * t90;
t13 = t15 * qJD(5);
t10 = t44 * t77;
t9 = t67 + t11 / 0.2e1 + t48 / 0.2e1;
t5 = t69 + t112;
t4 = t68 - (t29 + t58) * t44 / 0.2e1;
t6 = [0, 0, 0, 0, 0, -t35, t34, -t64 * t35, t97, qJD(3) * t1 + qJD(4) * t11, -t20, t13, 0, 0, 0, t29 * t87 + t18, -t29 * t40 + t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -t35 - t92, t34 + t93, -t86 * t64 * t38, -t94 + t97, t96 + (-t38 * pkin(3) - qJ(4) * t14) * qJD(3) + t9 * qJD(4), -t20, t13, 0, 0, 0, t5 * qJD(5) + t18 + t83, t4 * qJD(5) + t19 + t82; 0, 0, 0, 0, 0, 0, 0, 0, t109, qJD(3) * t9 + t95, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, t111, -t40, -t87, 0, t84 + t5 * qJD(3) + (-t100 * t22 - t65 * t21) * qJD(5), -t85 + t4 * qJD(3) + (-t100 * t21 + t65 * t22) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t87, t40; 0, 0, 0, 0, 0, t92, -t93, t64 * t92, t51 + t94, -qJD(4) * t8 - t96, -t20, t13, 0, 0, 0, -t2 * qJD(5) - t83, -t3 * qJD(5) - t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, t51, t48 * qJD(4), -t20, t13, 0, 0, 0, t58 * t87, -t58 * t40; 0, 0, 0, 0, 0, 0, 0, 0, t109, -t73, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, t111, -t40, -t87, 0, (-t100 * t50 - t65 * t49) * qJD(5) + t71, (-t100 * t49 + t65 * t50) * qJD(5) + t70; 0, 0, 0, 0, 0, 0, 0, 0, -t109, qJD(3) * t8 - t95, 0, 0, 0, 0, 0, t87, -t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, -t109, t73, 0, 0, 0, 0, 0, t87, -t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77, -t110; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t111, 0, 0, 0, t2 * qJD(3) - t43 - t84, t3 * qJD(3) + t39 + t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t111, 0, 0, 0, -t43 - t71, t39 - t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t77, t110; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t6;
