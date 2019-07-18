% Calculate inertial parameters regressor of coriolis matrix for
% S4RRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
% 
% Output:
% cmat_reg [(4*4)x(4*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 18:16
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4RRPR2_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR2_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR2_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPR2_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 18:16:37
% EndTime: 2019-07-18 18:16:38
% DurationCPUTime: 0.43s
% Computational Cost: add. (797->101), mult. (1169->118), div. (0->0), fcn. (853->4), ass. (0->75)
t72 = sin(qJ(2));
t108 = t72 * pkin(1);
t73 = cos(qJ(4));
t57 = t73 * t108;
t110 = t57 / 0.2e1;
t109 = -pkin(2) - pkin(3);
t74 = cos(qJ(2));
t107 = t74 * pkin(1);
t58 = qJ(3) + t108;
t71 = sin(qJ(4));
t106 = t71 * t58;
t105 = t71 * t72;
t104 = t73 * t58;
t103 = t73 * t74;
t86 = t71 * t107;
t37 = -t57 + t86;
t65 = t71 * qJD(3);
t102 = t37 * qJD(2) + t65;
t38 = (t103 + t105) * pkin(1);
t66 = t73 * qJD(3);
t101 = t38 * qJD(2) + t66;
t100 = pkin(1) * qJD(1);
t99 = pkin(1) * qJD(2);
t79 = -pkin(2) - t107;
t77 = -pkin(3) + t79;
t48 = t73 * t77;
t26 = -t48 + t106;
t76 = t71 * t77;
t27 = t76 + t104;
t6 = t26 * t37 + t27 * t38;
t98 = t6 * qJD(1);
t97 = t71 * qJ(3);
t96 = t73 * qJ(3);
t78 = -qJ(3) / 0.2e1 - t58 / 0.2e1;
t84 = t71 * t109;
t9 = t78 * t73 + t110 - t84;
t95 = t9 * qJD(1);
t61 = t73 * t109;
t80 = -t48 / 0.2e1 - t61 / 0.2e1;
t83 = -t103 / 0.2e1;
t10 = pkin(1) * t83 + (-t108 / 0.2e1 - t78) * t71 + t80;
t94 = t10 * qJD(1);
t13 = t26 * t71 + t27 * t73;
t93 = t13 * qJD(1);
t25 = (t58 * t74 + t79 * t72) * pkin(1);
t92 = t25 * qJD(1);
t91 = t37 * qJD(1);
t90 = t38 * qJD(1);
t64 = t74 * t99;
t89 = t64 + qJD(3);
t69 = qJD(1) + qJD(2);
t88 = qJD(1) - qJD(4);
t87 = qJD(2) - qJD(4);
t85 = t72 * t99;
t42 = -t61 + t97;
t82 = t42 / 0.2e1 + t26 / 0.2e1;
t43 = t84 + t96;
t81 = t43 / 0.2e1 + t27 / 0.2e1;
t46 = t69 * t71;
t47 = t69 * t73;
t21 = t71 * t42 + t73 * t43;
t4 = (t37 / 0.2e1 + t81) * t73 + (-t38 / 0.2e1 + t82) * t71;
t75 = -t4 * qJD(1) - t21 * qJD(2);
t70 = qJ(3) * qJD(3);
t63 = t74 * t100;
t62 = t72 * t100;
t54 = t69 * qJ(3);
t53 = t58 * qJD(3);
t39 = -t62 - t85;
t36 = -qJD(4) * t73 + t47;
t35 = -qJD(4) * t71 + t46;
t12 = t96 / 0.2e1 + t84 / 0.2e1 + t104 / 0.2e1 + t76 / 0.2e1 - t86 / 0.2e1 + t110;
t11 = -t97 / 0.2e1 - t106 / 0.2e1 + (t83 - t105 / 0.2e1) * pkin(1) - t80;
t3 = (-t37 / 0.2e1 + t81) * t73 + (t38 / 0.2e1 + t82) * t71;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t85, -t64, 0, 0, 0, 0, 0, 0, 0, 0, -t85, 0, t89, t25 * qJD(2) + t53, 0, 0, 0, 0, 0, 0, t27 * qJD(4) + t102, -t26 * qJD(4) + t101, 0, t6 * qJD(2) + t13 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, -t63 - t64, 0, 0, 0, 0, 0, 0, 0, 0, t39, 0, t63 + t89, t92 + t53 + (-pkin(2) * t72 + qJ(3) * t74) * t99, 0, 0, 0, 0, 0, 0, t12 * qJD(4) + t102 + t91, t11 * qJD(4) + t101 + t90, 0, t98 + (t37 * t42 + t38 * t43) * qJD(2) + t3 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, t69 * t58, 0, 0, 0, 0, 0, 0, t46, t47, 0, t3 * qJD(2) + t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12 * qJD(2) + t88 * t27, t11 * qJD(2) - t88 * t26, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, t63, 0, 0, 0, 0, 0, 0, 0, 0, t62, 0, -t63 + qJD(3), t70 - t92, 0, 0, 0, 0, 0, 0, -t9 * qJD(4) + t65 - t91, -t10 * qJD(4) + t66 - t90, 0, t4 * qJD(3) - t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3), t70, 0, 0, 0, 0, 0, 0, t43 * qJD(4) + t65, -t42 * qJD(4) + t66, 0, t21 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, t54, 0, 0, 0, 0, 0, 0, t46, t47, 0, -t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t87 * t43 - t95, -t87 * t42 - t94, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t69, -qJ(3) * qJD(2) - t58 * qJD(1), 0, 0, 0, 0, 0, 0, -t35, -t36, 0, -t4 * qJD(2) - t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t69, -t54, 0, 0, 0, 0, 0, 0, -t35, -t36, 0, t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, t36, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t27 * qJD(1) + t9 * qJD(2) - t65, t26 * qJD(1) + t10 * qJD(2) - t66, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43 * qJD(2) - t65 + t95, t42 * qJD(2) - t66 + t94, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46, -t47, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t1;
