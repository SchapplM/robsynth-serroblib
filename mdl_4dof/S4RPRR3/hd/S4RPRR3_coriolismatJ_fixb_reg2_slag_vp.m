% Calculate inertial parameters regressor of coriolis matrix for
% S4RPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% 
% Output:
% cmat_reg [(4*4)x(4*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4RPRR3_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR3_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR3_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR3_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:49:21
% EndTime: 2019-12-31 16:49:22
% DurationCPUTime: 0.51s
% Computational Cost: add. (1176->55), mult. (2226->76), div. (0->0), fcn. (2239->6), ass. (0->47)
t87 = qJD(3) + qJD(4);
t100 = cos(qJ(4));
t70 = sin(pkin(7)) * pkin(1) + pkin(5);
t101 = pkin(6) + t70;
t75 = sin(qJ(3));
t62 = t101 * t75;
t76 = cos(qJ(3));
t63 = t101 * t76;
t74 = sin(qJ(4));
t80 = t100 * t62 + t74 * t63;
t107 = t87 * t80;
t41 = t100 * t63 - t74 * t62;
t106 = t87 * t41;
t64 = -t100 * t76 + t74 * t75;
t105 = t87 * t64;
t103 = pkin(3) * t74;
t102 = pkin(3) * t75;
t99 = qJD(3) * pkin(3);
t71 = -cos(pkin(7)) * pkin(1) - pkin(2);
t67 = -t76 * pkin(3) + t71;
t96 = qJD(1) * t67;
t95 = qJD(1) * t76;
t94 = qJD(4) * t67;
t66 = t100 * t75 + t74 * t76;
t31 = t64 * t102 + t67 * t66;
t93 = t31 * qJD(1);
t32 = t66 * t102 - t67 * t64;
t92 = t32 * qJD(1);
t33 = t64 ^ 2 - t66 ^ 2;
t91 = t33 * qJD(1);
t68 = -t75 ^ 2 + t76 ^ 2;
t90 = t68 * qJD(1);
t89 = t75 * qJD(3);
t88 = t76 * qJD(3);
t86 = t66 * t64 * qJD(1);
t85 = t64 * t96;
t84 = t66 * t96;
t83 = t71 * t75 * qJD(1);
t82 = t71 * t95;
t81 = t75 * t88;
t79 = t100 * qJD(3);
t78 = t100 * qJD(4);
t46 = t87 * t66;
t7 = t67 * t102;
t77 = t7 * qJD(1);
t69 = t75 * t95;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81, t68 * qJD(3), 0, -t81, 0, 0, t71 * t89, t71 * t88, 0, 0, -t64 * t46, t87 * t33, 0, t66 * t105, 0, 0, t31 * qJD(3) + t66 * t94, t32 * qJD(3) - t64 * t94, 0, t7 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, t90, t88, -t69, -t89, 0, -t70 * t88 + t83, t70 * t89 + t82, 0, 0, -t86, t91, -t105, t86, -t46, 0, -t106 + t93, t107 + t92, (t100 * t64 - t66 * t74) * t99, (-t100 * t41 - t74 * t80) * t99 + t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t86, t91, -t105, t86, -t46, 0, -t106 + t84, t107 - t85, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t89, -t88, 0, 0, 0, 0, 0, 0, 0, 0, -t46, t105, 0, (-t100 * t66 - t64 * t74) * t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46, t105, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t69, -t90, 0, t69, 0, 0, -t83, -t82, 0, 0, t86, -t91, 0, -t86, 0, 0, -t93, -t92, 0, -t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(4) * t103, -pkin(3) * t78, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t87 * t103, (-t79 - t78) * pkin(3), 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86, -t91, 0, -t86, 0, 0, -t84, t85, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74 * t99, pkin(3) * t79, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t1;
