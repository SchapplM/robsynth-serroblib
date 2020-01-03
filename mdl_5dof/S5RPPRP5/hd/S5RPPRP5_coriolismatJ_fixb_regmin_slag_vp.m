% Calculate minimal parameter regressor of coriolis matrix for
% S5RPPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x22]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RPPRP5_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP5_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP5_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP5_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:53:45
% EndTime: 2019-12-31 17:53:47
% DurationCPUTime: 0.48s
% Computational Cost: add. (598->90), mult. (1205->111), div. (0->0), fcn. (1278->4), ass. (0->67)
t63 = sin(pkin(7));
t64 = cos(pkin(7));
t65 = cos(qJ(4));
t98 = sin(qJ(4));
t42 = t63 * t98 + t64 * t65;
t100 = t42 ^ 2;
t44 = t63 * t65 - t64 * t98;
t40 = t44 ^ 2;
t99 = -t63 / 0.2e1;
t97 = -pkin(6) + qJ(2);
t61 = t63 ^ 2;
t55 = t64 ^ 2 + t61;
t49 = -t64 * pkin(2) - t63 * qJ(3) - pkin(1);
t31 = t64 * pkin(3) - t49;
t67 = t42 * pkin(4) - t44 * qJ(5);
t10 = t67 + t31;
t68 = t44 * pkin(4) + t42 * qJ(5);
t1 = t10 * t68;
t96 = t1 * qJD(1);
t3 = t10 * t44 + t42 * t68;
t95 = t3 * qJD(1);
t4 = t10 * t42 - t44 * t68;
t94 = t4 * qJD(1);
t51 = t97 * t64;
t70 = t97 * t63;
t21 = t98 * t51 - t65 * t70;
t22 = t65 * t51 + t98 * t70;
t5 = -t21 * t44 + t22 * t42;
t93 = t5 * qJD(1);
t8 = -t40 - t100;
t92 = t8 * qJD(1);
t91 = t68 * qJD(1);
t11 = -t40 + t100;
t90 = t11 * qJD(1);
t66 = -t42 * t98 / 0.2e1 - t44 * t65 / 0.2e1;
t13 = t99 + t66;
t89 = t13 * qJD(1);
t88 = t21 * qJD(4);
t20 = t22 * qJD(4);
t87 = t40 * qJD(1);
t86 = t42 * qJD(1);
t33 = t42 * qJD(4);
t85 = t44 * qJD(1);
t37 = t44 * qJD(4);
t84 = t44 * qJD(5);
t50 = t55 * qJ(2);
t83 = t50 * qJD(1);
t82 = t55 * qJD(1);
t81 = t61 * qJD(1);
t80 = t63 * qJD(1);
t79 = t63 * qJD(3);
t78 = t65 * qJD(4);
t77 = qJD(4) * qJ(5);
t76 = t10 * t80;
t75 = t31 * t85;
t74 = t42 * t85;
t73 = t44 * t79;
t72 = t44 * t80;
t71 = t64 * t80;
t69 = t98 * qJD(4);
t52 = t55 * qJD(2);
t46 = t50 * qJD(2);
t35 = t44 * qJD(2);
t24 = t42 * t79;
t23 = t42 * t80;
t12 = t99 - t66;
t2 = [0, 0, 0, 0, 0, t52, t46, t64 * t79, t52, t61 * qJD(3), -t49 * t79 + t46, -t42 * t37, t11 * qJD(4), 0, 0, 0, t31 * t37 + t24, -t31 * t33 + t73, t3 * qJD(4) - t42 * t84 + t24, t8 * qJD(2), t4 * qJD(4) + t40 * qJD(5) - t73, t5 * qJD(2) + t1 * qJD(4) + (t79 - t84) * t10; 0, 0, 0, 0, 0, t82, t83, 0, t82, 0, t83, 0, 0, 0, 0, 0, 0, 0, 0, t92, 0, t12 * qJD(3) + t93; 0, 0, 0, 0, 0, 0, 0, t71, 0, t81, -t49 * t80, 0, 0, 0, 0, 0, t23, t72, t23, 0, -t72, t12 * qJD(2) + t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t74, t90, -t33, -t37, 0, -t20 + t75, -t31 * t86 + t88, -t20 + t95, t67 * qJD(4) - t42 * qJD(5), -t88 + t94, t96 + (-t22 * pkin(4) - t21 * qJ(5)) * qJD(4) + t22 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t74, -t33, t87, -t10 * t85 + t20; 0, 0, 0, 0, 0, -t82, -t83, 0, -t82, 0, -t83 - t79, 0, 0, 0, 0, 0, -t37, t33, -t37, -t92, -t33, t13 * qJD(3) - qJD(4) * t68 + t84 - t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t80, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t85, t86, -t85, 0, -t86, -t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t85; 0, 0, 0, 0, 0, 0, 0, -t71, 0, -t81, (qJD(1) * t49 + qJD(2)) * t63, 0, 0, 0, 0, 0, -t23, -t72, -t23, 0, t72, -t13 * qJD(2) - t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t69, -t78, -t69, 0, t78, (-t98 * pkin(4) + t65 * qJ(5)) * qJD(4) + t98 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74, -t90, 0, 0, 0, t35 - t75, (qJD(1) * t31 - qJD(2)) * t42, t35 - t95, 0, t42 * qJD(2) - t94, qJD(2) * t68 - t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t85, -t86, t85, 0, t86, t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(5), qJ(5) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4), t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74, 0, -t87, (qJD(1) * t10 - qJD(2)) * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(4), -t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t2;
