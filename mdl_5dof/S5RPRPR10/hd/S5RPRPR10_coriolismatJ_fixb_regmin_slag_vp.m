% Calculate minimal parameter regressor of coriolis matrix for
% S5RPRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x17]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RPRPR10_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR10_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR10_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR10_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:26:06
% EndTime: 2019-12-31 18:26:08
% DurationCPUTime: 0.54s
% Computational Cost: add. (836->96), mult. (1303->133), div. (0->0), fcn. (1256->6), ass. (0->78)
t103 = sin(qJ(3));
t104 = cos(qJ(3));
t105 = -pkin(1) - pkin(2);
t48 = t103 * qJ(2) - t104 * t105;
t49 = t104 * qJ(2) + t103 * t105;
t95 = sin(pkin(8));
t74 = t95 * t49;
t96 = cos(pkin(8));
t30 = -t96 * t48 - t74;
t108 = -t30 / 0.2e1;
t47 = -pkin(3) - t48;
t27 = t96 * t47 - t74;
t25 = pkin(4) - t27;
t82 = t96 * pkin(3);
t63 = -t82 - pkin(4);
t114 = -t25 / 0.2e1 + t63 / 0.2e1;
t113 = t108 - t114;
t43 = t96 * t49;
t29 = -t95 * t48 + t43;
t85 = qJD(1) - qJD(3);
t112 = t29 * t85;
t65 = sin(qJ(5));
t66 = cos(qJ(5));
t56 = -t65 ^ 2 + t66 ^ 2;
t111 = t85 * t56;
t46 = -t96 * t103 - t95 * t104;
t102 = t27 * t46;
t28 = t95 * t47 + t43;
t45 = t95 * t103 - t96 * t104;
t101 = t28 * t45;
t100 = pkin(3) * qJD(3);
t81 = t95 * pkin(3);
t1 = (t82 / 0.2e1 + t30 / 0.2e1 - t27 / 0.2e1) * t46 + (-t81 / 0.2e1 + t28 / 0.2e1 - t29 / 0.2e1) * t45;
t99 = t1 * qJD(1);
t3 = -t27 * t29 + t28 * t30;
t98 = t3 * qJD(1);
t4 = -t101 + t102;
t97 = t4 * qJD(1);
t94 = qJD(1) * t65;
t93 = qJD(1) * t66;
t92 = qJD(2) * t46;
t91 = qJD(3) * t65;
t90 = qJD(3) * t66;
t89 = qJD(5) * t65;
t88 = qJD(5) * t66;
t87 = t56 * qJD(5);
t86 = qJ(2) * qJD(1);
t84 = t25 * t94;
t83 = t25 * t93;
t80 = t104 * qJD(1);
t79 = t104 * qJD(2);
t78 = t96 * t46;
t77 = t103 * qJD(1);
t76 = t103 * qJD(2);
t75 = t95 * t45;
t73 = t108 + t114;
t72 = qJD(1) * t29 - t92;
t71 = -qJD(3) * t29 + t92;
t5 = t73 * t65;
t70 = t5 * qJD(1) - t63 * t91;
t6 = t73 * t66;
t69 = t6 * qJD(1) - t63 * t90;
t20 = t45 * t65;
t40 = t46 * t93;
t68 = t20 * qJD(5) + t46 * t90 - t40;
t22 = t45 * t66;
t39 = t46 * t94;
t67 = -t22 * qJD(5) + t46 * t91 - t39;
t62 = t81 + pkin(7);
t58 = t65 * t88;
t53 = -qJD(3) * t104 + t80;
t52 = -qJD(3) * t103 + t77;
t42 = (-t90 + t93) * t65;
t26 = -pkin(7) + t28;
t8 = t113 * t66;
t7 = t113 * t65;
t2 = t46 * t108 - t101 / 0.2e1 + t29 * t45 / 0.2e1 + t102 / 0.2e1 + (-t75 / 0.2e1 + t78 / 0.2e1) * pkin(3);
t9 = [0, 0, 0, 0, qJD(2), qJ(2) * qJD(2), 0, t49 * qJD(3) + t76, -t48 * qJD(3) + t79, t4 * qJD(2) + t3 * qJD(3), t58, t87, 0, 0, 0, -t25 * t89 - t71 * t66, -t25 * t88 + t71 * t65; 0, 0, 0, 0, qJD(1), t86, 0, t77, t80, t2 * qJD(3) + t97, 0, 0, 0, 0, 0, -t40, t39; 0, 0, 0, 0, 0, 0, 0, t85 * t49, -t85 * t48, t98 + t2 * qJD(2) + (-t96 * t29 + t95 * t30) * t100, -t58, -t87, 0, 0, 0, t7 * qJD(5) + t66 * t112, t8 * qJD(5) - t65 * t112; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, t111, -t88, t89, 0, t7 * qJD(3) - t26 * t88 - t84, t8 * qJD(3) + t26 * t89 - t83; 0, 0, 0, 0, -qJD(1), -t86, 0, -t52, -t53, -t1 * qJD(3) - t97, 0, 0, 0, 0, 0, -t68, t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, t52, t53, -t99 + (-t75 + t78) * t100, 0, 0, 0, 0, 0, t68, -t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t85 * t20 + t46 * t88, -t85 * t22 - t46 * t89; 0, 0, 0, 0, 0, 0, 0, -t49 * qJD(1) - t76, t48 * qJD(1) - t79, t1 * qJD(2) - t98, -t58, -t87, 0, 0, 0, -t5 * qJD(5) - t72 * t66, -t6 * qJD(5) + t72 * t65; 0, 0, 0, 0, 0, 0, 0, -t77, -t80, t99, 0, 0, 0, 0, 0, t40, -t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, t87, 0, 0, 0, t63 * t89, t63 * t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t42, -t111, t88, -t89, 0, -t62 * t88 - t70, t62 * t89 - t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t89, -t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t42, -t111, 0, 0, 0, t20 * qJD(2) + t5 * qJD(3) + t84, t22 * qJD(2) + t6 * qJD(3) + t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20 * qJD(1), t22 * qJD(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, t111, 0, 0, 0, t70, t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t9;
