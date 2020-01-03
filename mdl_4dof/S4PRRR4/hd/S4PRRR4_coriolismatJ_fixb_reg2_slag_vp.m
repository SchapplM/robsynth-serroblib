% Calculate inertial parameters regressor of coriolis matrix for
% S4PRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% 
% Output:
% cmat_reg [(4*4)x(4*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4PRRR4_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR4_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR4_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR4_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:32:42
% EndTime: 2019-12-31 16:32:43
% DurationCPUTime: 0.49s
% Computational Cost: add. (897->53), mult. (1947->74), div. (0->0), fcn. (1960->4), ass. (0->45)
t84 = qJD(3) + qJD(4);
t100 = -pkin(6) - pkin(5);
t72 = sin(qJ(3));
t65 = t100 * t72;
t73 = cos(qJ(3));
t66 = t100 * t73;
t71 = sin(qJ(4));
t97 = cos(qJ(4));
t77 = -t65 * t97 - t66 * t71;
t104 = t84 * t77;
t50 = t65 * t71 - t66 * t97;
t103 = t84 * t50;
t61 = t71 * t72 - t73 * t97;
t102 = t84 * t61;
t99 = pkin(3) * t71;
t98 = pkin(3) * t72;
t96 = qJD(3) * pkin(3);
t70 = -pkin(3) * t73 - pkin(2);
t93 = qJD(2) * t70;
t92 = qJD(2) * t73;
t91 = qJD(4) * t70;
t63 = t71 * t73 + t72 * t97;
t27 = t61 ^ 2 - t63 ^ 2;
t90 = t27 * qJD(2);
t32 = t61 * t98 + t63 * t70;
t89 = t32 * qJD(2);
t33 = -t61 * t70 + t63 * t98;
t88 = t33 * qJD(2);
t67 = -t72 ^ 2 + t73 ^ 2;
t87 = t67 * qJD(2);
t86 = t72 * qJD(3);
t85 = t73 * qJD(3);
t83 = pkin(2) * t72 * qJD(2);
t82 = pkin(2) * t92;
t81 = t63 * t61 * qJD(2);
t80 = t61 * t93;
t79 = t63 * t93;
t78 = t72 * t85;
t76 = t97 * qJD(3);
t75 = t97 * qJD(4);
t40 = t84 * t63;
t7 = t70 * t98;
t74 = t7 * qJD(2);
t68 = t72 * t92;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t86, -t85, 0, 0, 0, 0, 0, 0, 0, 0, -t40, t102, 0, (-t61 * t71 - t63 * t97) * t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t40, t102, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t78, t67 * qJD(3), 0, -t78, 0, 0, -pkin(2) * t86, -pkin(2) * t85, 0, 0, -t61 * t40, t84 * t27, 0, t63 * t102, 0, 0, qJD(3) * t32 + t63 * t91, qJD(3) * t33 - t61 * t91, 0, t7 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68, t87, t85, -t68, -t86, 0, -pkin(5) * t85 - t83, pkin(5) * t86 - t82, 0, 0, -t81, t90, -t102, t81, -t40, 0, -t103 + t89, t104 + t88, (t61 * t97 - t63 * t71) * t96, (-t50 * t97 - t71 * t77) * t96 + t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t81, t90, -t102, t81, -t40, 0, -t103 + t79, t104 - t80, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t68, -t87, 0, t68, 0, 0, t83, t82, 0, 0, t81, -t90, 0, -t81, 0, 0, -t89, -t88, 0, -t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(4) * t99, -pkin(3) * t75, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t84 * t99, (-t76 - t75) * pkin(3), 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81, -t90, 0, -t81, 0, 0, -t79, t80, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71 * t96, pkin(3) * t76, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t1;
