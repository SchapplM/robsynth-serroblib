% Calculate minimal parameter regressor of coriolis matrix for
% S4RRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,theta3]';
% 
% Output:
% cmat_reg [(4*%NQJ)%x16]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4RRPP3_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP3_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP3_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP3_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:57:52
% EndTime: 2019-12-31 16:57:53
% DurationCPUTime: 0.37s
% Computational Cost: add. (596->79), mult. (1194->125), div. (0->0), fcn. (1206->4), ass. (0->65)
t50 = sin(pkin(6));
t51 = sin(qJ(2));
t52 = cos(qJ(2));
t77 = cos(pkin(6));
t37 = t50 * t52 + t77 * t51;
t34 = t37 ^ 2;
t35 = t50 * t51 - t77 * t52;
t97 = t35 ^ 2 + t34;
t100 = qJD(3) * t97;
t99 = t97 * qJD(1);
t98 = qJ(3) + pkin(5);
t41 = t98 * t51;
t42 = t98 * t52;
t20 = t77 * t41 + t50 * t42;
t40 = t77 * t42;
t89 = t50 * t41;
t59 = t40 - t89;
t56 = t20 * t37 - t35 * t59;
t95 = qJD(3) * t56;
t93 = t56 * qJD(1);
t92 = -t35 / 0.2e1;
t58 = t40 / 0.2e1;
t91 = t37 * pkin(3);
t90 = t51 * pkin(2);
t88 = qJD(2) * pkin(2);
t60 = -t52 * pkin(2) - pkin(1);
t14 = t35 * pkin(3) - t37 * qJ(4) + t60;
t84 = t35 * qJ(4);
t15 = t84 + t90 + t91;
t1 = t14 * t15;
t87 = t1 * qJD(1);
t4 = t60 * t90;
t83 = t4 * qJD(1);
t5 = t14 * t37 + t15 * t35;
t82 = t5 * qJD(1);
t6 = t14 * t35 - t15 * t37;
t81 = t6 * qJD(1);
t44 = t50 * pkin(2) + qJ(4);
t46 = -t77 * pkin(2) - pkin(3);
t47 = t90 / 0.2e1;
t9 = t47 + (pkin(3) / 0.2e1 - t46 / 0.2e1) * t37 + (qJ(4) / 0.2e1 + t44 / 0.2e1) * t35;
t78 = t9 * qJD(1);
t76 = qJD(1) * t52;
t54 = t50 * t92 - t77 * t37 / 0.2e1;
t12 = (-t51 / 0.2e1 + t54) * pkin(2);
t75 = t12 * qJD(1);
t72 = t34 * qJD(1);
t71 = t35 * qJD(1);
t70 = t35 * qJD(2);
t69 = t37 * qJD(1);
t68 = t37 * qJD(4);
t43 = -t51 ^ 2 + t52 ^ 2;
t67 = t43 * qJD(1);
t66 = t51 * qJD(2);
t65 = t52 * qJD(2);
t64 = pkin(1) * t51 * qJD(1);
t63 = pkin(1) * t76;
t62 = t35 * t69;
t61 = t51 * t76;
t18 = t58 - t40 / 0.2e1;
t55 = t18 * qJD(1) + t44 * qJD(2);
t13 = 0.2e1 * t58 - t89;
t11 = t54 * pkin(2) + t47;
t10 = t44 * t92 + t46 * t37 / 0.2e1 + t47 + t84 / 0.2e1 + t91 / 0.2e1;
t2 = [0, 0, 0, t51 * t65, t43 * qJD(2), 0, 0, 0, -pkin(1) * t66, -pkin(1) * t65, t100, t4 * qJD(2) + t95, t5 * qJD(2) - t35 * t68, t100, t6 * qJD(2) + t34 * qJD(4), t1 * qJD(2) - t14 * t68 + t95; 0, 0, 0, t61, t67, t65, -t66, 0, -pkin(5) * t65 - t64, pkin(5) * t66 - t63, (t77 * t35 - t37 * t50) * t88, t83 + (-t20 * t50 - t59 * t77) * t88 + t11 * qJD(3), -qJD(2) * t59 + t82, (-t46 * t35 - t44 * t37) * qJD(2) - qJD(4) * t35, -qJD(2) * t20 + t81, t87 + (-t20 * t44 + t46 * t59) * qJD(2) + t10 * qJD(3) + t13 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t99, t11 * qJD(2) + t93, 0, t99, 0, t10 * qJD(2) + t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t62, -t70, t72, t13 * qJD(2) - t14 * t69; 0, 0, 0, -t61, -t67, 0, 0, 0, t64, t63, 0, t12 * qJD(3) - t83, -t37 * qJD(3) - t82, 0, -t35 * qJD(3) - t81, -t9 * qJD(3) + t18 * qJD(4) - t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4), t44 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, -t69, 0, -t71, -t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2), t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t99, -t12 * qJD(2) - t93, t37 * qJD(2), -t99, t70, t9 * qJD(2) - t68 - t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t75, t69, 0, t71, t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, 0, -t72, -t18 * qJD(2) + (qJD(1) * t14 + qJD(3)) * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(2), -t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t2;
