% Calculate minimal parameter regressor of coriolis matrix for
% S4RRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% cmat_reg [(4*%NQJ)%x15]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4RRRP2_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP2_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP2_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP2_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:13:10
% EndTime: 2019-12-31 17:13:11
% DurationCPUTime: 0.38s
% Computational Cost: add. (527->86), mult. (1097->124), div. (0->0), fcn. (759->4), ass. (0->81)
t69 = cos(qJ(2));
t107 = t69 * pkin(1);
t57 = -pkin(2) - t107;
t112 = t57 / 0.2e1 - pkin(2) / 0.2e1;
t66 = sin(qJ(3));
t64 = t66 ^ 2;
t68 = cos(qJ(3));
t65 = t68 ^ 2;
t53 = t65 + t64;
t89 = qJD(1) + qJD(2);
t114 = t89 * t53;
t54 = t65 - t64;
t113 = t89 * t54;
t109 = pkin(3) * t66;
t67 = sin(qJ(2));
t108 = t67 * pkin(1);
t56 = pkin(6) + t108;
t94 = qJ(4) + t56;
t34 = t94 * t66;
t106 = t34 * t66;
t35 = t94 * t68;
t105 = t35 * t68;
t102 = -qJ(4) - pkin(6);
t45 = t102 * t66;
t104 = t45 * t66;
t46 = t102 * t68;
t103 = t46 * t68;
t33 = t53 * t107;
t47 = t53 * qJD(4);
t101 = t33 * qJD(2) + t47;
t100 = pkin(1) * qJD(1);
t99 = pkin(1) * qJD(2);
t98 = pkin(2) * qJD(2);
t97 = qJD(3) * pkin(3);
t58 = -t68 * pkin(3) - pkin(2);
t44 = t58 - t107;
t3 = t44 * t109;
t96 = t3 * qJD(1);
t11 = t105 + t106;
t8 = (t11 * t69 + t44 * t67) * pkin(1);
t95 = t8 * qJD(1);
t93 = qJD(1) * t57;
t92 = qJD(3) * t66;
t63 = qJD(3) * t68;
t91 = t11 * qJD(1);
t90 = t33 * qJD(1);
t88 = t67 * t99;
t87 = pkin(3) * t63;
t86 = qJD(4) * t109;
t85 = t67 * t100;
t83 = t108 / 0.2e1;
t82 = -t107 / 0.2e1;
t81 = t66 * t93;
t80 = t68 * t93;
t79 = pkin(1) * t89;
t78 = t66 * t85;
t77 = t66 * t82;
t76 = t89 * t66;
t75 = t67 * t79;
t13 = -t103 - t104;
t1 = (t82 - t58 / 0.2e1 - t44 / 0.2e1) * t109;
t9 = t58 * t109;
t74 = -t1 * qJD(1) + t9 * qJD(2);
t4 = t83 + (t46 / 0.2e1 - t35 / 0.2e1) * t68 + (t45 / 0.2e1 - t34 / 0.2e1) * t66;
t73 = -t4 * qJD(1) + t13 * qJD(2);
t72 = t82 - t112;
t18 = t72 * t66;
t71 = t18 * qJD(1) + t66 * t98;
t19 = t72 * t68;
t70 = t19 * qJD(1) + t68 * t98;
t62 = pkin(3) * t92;
t55 = t66 * t63;
t52 = t66 * t88;
t48 = t54 * qJD(3);
t40 = pkin(3) * t76;
t32 = t68 * t76;
t21 = (t82 + t112) * t68;
t20 = t112 * t66 + t77;
t5 = -t103 / 0.2e1 + t105 / 0.2e1 - t104 / 0.2e1 + t106 / 0.2e1 + t83;
t2 = pkin(3) * t77 + (t44 + t58) * t109 / 0.2e1;
t6 = [0, 0, 0, 0, -t88, -t69 * t99, t55, t48, 0, 0, 0, t57 * t92 - t68 * t88, t57 * t63 + t52, t101, t8 * qJD(2) + t3 * qJD(3) + t11 * qJD(4); 0, 0, 0, 0, -t75, -t69 * t79, t55, t48, 0, 0, 0, t20 * qJD(3) - t68 * t75, t21 * qJD(3) + t52 + t78, t90 + t101, t95 + t2 * qJD(3) + t5 * qJD(4) + (t13 * t69 + t58 * t67) * t99; 0, 0, 0, 0, 0, 0, t32, t113, t63, -t92, 0, t20 * qJD(2) - t56 * t63 + t81, t21 * qJD(2) + t56 * t92 + t80, -t87, t2 * qJD(2) - t35 * t97 + t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t114, t5 * qJD(2) + t91; 0, 0, 0, 0, t85, t69 * t100, t55, t48, 0, 0, 0, -t18 * qJD(3) + t68 * t85, -t19 * qJD(3) - t78, t47 - t90, -t1 * qJD(3) - t4 * qJD(4) - t95; 0, 0, 0, 0, 0, 0, t55, t48, 0, 0, 0, -pkin(2) * t92, -pkin(2) * t63, t47, t9 * qJD(3) + t13 * qJD(4); 0, 0, 0, 0, 0, 0, t32, t113, t63, -t92, 0, -pkin(6) * t63 - t71, pkin(6) * t92 - t70, -t87, t46 * t97 + t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t114, t73; 0, 0, 0, 0, 0, 0, -t32, -t113, 0, 0, 0, t18 * qJD(2) - t81, t19 * qJD(2) - t80, 0, t1 * qJD(2) - t86 - t96; 0, 0, 0, 0, 0, 0, -t32, -t113, 0, 0, 0, t71, t70, 0, -t74 - t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t114, t4 * qJD(2) + t62 - t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t114, t62 - t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t6;
