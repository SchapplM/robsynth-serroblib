% Calculate inertial parameters regressor of coriolis matrix for
% S4RPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
% 
% Output:
% cmat_reg [(4*4)x(4*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4RPRR8_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR8_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR8_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR8_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:55:19
% EndTime: 2019-12-31 16:55:19
% DurationCPUTime: 0.56s
% Computational Cost: add. (1200->69), mult. (2031->80), div. (0->0), fcn. (2089->4), ass. (0->53)
t87 = qJD(3) + qJD(4);
t76 = -pkin(1) - pkin(5);
t110 = -pkin(6) + t76;
t73 = sin(qJ(3));
t63 = t110 * t73;
t75 = cos(qJ(3));
t64 = t110 * t75;
t72 = sin(qJ(4));
t74 = cos(qJ(4));
t80 = t72 * t63 - t74 * t64;
t112 = t87 * t80;
t45 = t74 * t63 + t72 * t64;
t111 = t87 * t45;
t62 = -t72 * t73 + t74 * t75;
t42 = t87 * t62;
t60 = t72 * t75 + t74 * t73;
t109 = t60 * t42;
t103 = qJD(3) * pkin(3);
t108 = (t60 * t74 - t62 * t72) * t103;
t105 = pkin(3) * t75;
t104 = pkin(3) * qJD(4);
t67 = t73 * pkin(3) + qJ(2);
t100 = qJD(4) * t67;
t30 = t60 ^ 2 - t62 ^ 2;
t99 = t30 * qJD(1);
t34 = t60 * t105 + t67 * t62;
t98 = t34 * qJD(1);
t35 = t62 * t105 - t67 * t60;
t97 = t35 * qJD(1);
t96 = t60 * qJD(1);
t95 = t62 * qJD(1);
t65 = t73 ^ 2 - t75 ^ 2;
t94 = t65 * qJD(1);
t93 = t67 * qJD(1);
t92 = t73 * qJD(1);
t91 = t73 * qJD(3);
t90 = t75 * qJD(1);
t89 = t75 * qJD(3);
t88 = qJ(2) * qJD(3);
t70 = qJD(1) * qJ(2);
t86 = t60 * t95;
t85 = t60 * t93;
t84 = t62 * t93;
t83 = t73 * t89;
t82 = t73 * t70;
t81 = t75 * t70;
t79 = pkin(3) * t87;
t8 = t67 * t105;
t77 = t8 * qJD(1);
t71 = qJ(2) * qJD(2);
t66 = t73 * t90;
t41 = t87 * t60;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2), t71, -t83, t65 * qJD(3), 0, t83, 0, 0, qJD(2) * t73 + t75 * t88, qJD(2) * t75 - t73 * t88, 0, t71, -t109, t87 * t30, 0, t109, 0, 0, t60 * qJD(2) + t34 * qJD(3) + t62 * t100, t62 * qJD(2) + t35 * qJD(3) - t60 * t100, 0, t67 * qJD(2) + t8 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(1), t70, 0, 0, 0, 0, 0, 0, t92, t90, 0, t70, 0, 0, 0, 0, 0, 0, t96, t95, 0, t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t66, t94, -t91, t66, -t89, 0, -t76 * t91 + t81, -t76 * t89 - t82, 0, 0, -t86, t99, -t41, t86, -t42, 0, -t111 + t98, t112 + t97, t108, (-t45 * t74 - t72 * t80) * t103 + t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t86, t99, -t41, t86, -t42, 0, -t111 + t84, t112 - t85, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(1), -t70, 0, 0, 0, 0, 0, 0, -t92, -t90, 0, -t70, 0, 0, 0, 0, 0, 0, -t96, -t95, 0, -t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t91, -t89, 0, 0, 0, 0, 0, 0, 0, 0, -t41, -t42, 0, -t108; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41, -t42, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66, -t94, 0, -t66, 0, 0, -t81, t82, 0, 0, t86, -t99, 0, -t86, 0, 0, -t98, -t97, 0, -t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t72 * t104, -t74 * t104, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t72 * t79, -t74 * t79, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86, -t99, 0, -t86, 0, 0, -t84, t85, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72 * t103, t74 * t103, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t1;
