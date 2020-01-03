% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S2RR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% qJD [2x1]
%   Generalized joint velocities
% qJDD [2x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d2]';
% 
% Output:
% tauB_reg [6x(3*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S2RR1_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(2,1),zeros(3,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR1_invdynB_fixb_reg2_snew_vp: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'S2RR1_invdynB_fixb_reg2_snew_vp: qJD has to be [2x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [2 1]), ...
  'S2RR1_invdynB_fixb_reg2_snew_vp: qJDD has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S2RR1_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S2RR1_invdynB_fixb_reg2_snew_vp: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:19:11
% EndTime: 2020-01-03 11:19:11
% DurationCPUTime: 0.24s
% Computational Cost: add. (243->84), mult. (598->123), div. (0->0), fcn. (400->4), ass. (0->68)
t74 = sin(qJ(1));
t76 = cos(qJ(1));
t61 = g(1) * t74 + g(3) * t76;
t48 = -qJDD(1) * pkin(1) - t61;
t73 = sin(qJ(2));
t75 = cos(qJ(2));
t44 = -t75 * g(2) + t48 * t73;
t45 = g(2) * t73 + t48 * t75;
t80 = t44 * t75 - t45 * t73;
t102 = pkin(1) * t80;
t101 = t74 * g(2);
t100 = t76 * g(2);
t71 = t73 ^ 2;
t78 = qJD(1) ^ 2;
t99 = t71 * t78;
t72 = t75 ^ 2;
t98 = t72 * t78;
t62 = g(1) * t76 - t74 * g(3);
t49 = pkin(1) * t78 + t62;
t97 = t73 * t49;
t67 = t75 * t78 * t73;
t59 = qJDD(2) + t67;
t96 = t73 * t59;
t60 = qJDD(2) - t67;
t95 = t73 * t60;
t94 = t75 * t49;
t93 = t75 * t59;
t92 = t75 * t60;
t91 = -t71 - t72;
t90 = t73 * qJDD(1);
t89 = t74 * qJDD(1);
t88 = t75 * qJDD(1);
t87 = t76 * qJDD(1);
t86 = qJD(1) * qJD(2);
t85 = t73 * t86;
t84 = t75 * t86;
t83 = -t76 * t61 + t62 * t74;
t82 = t74 * t67;
t81 = t76 * t67;
t32 = t44 * t73 + t45 * t75;
t79 = t61 * t74 + t62 * t76;
t77 = qJD(2) ^ 2;
t66 = -t77 - t98;
t65 = -t77 + t98;
t64 = -t77 - t99;
t63 = t77 - t99;
t58 = t91 * t78;
t57 = (t71 - t72) * t78;
t56 = -t74 * t78 + t87;
t55 = t76 * t78 + t89;
t54 = t91 * qJDD(1);
t53 = t85 - t88;
t52 = -0.2e1 * t85 + t88;
t51 = -0.2e1 * t84 - t90;
t50 = -t84 - t90;
t46 = t91 * t86;
t43 = t50 * t75 + t71 * t86;
t42 = -t53 * t73 + t72 * t86;
t41 = -t64 * t73 - t92;
t40 = -t63 * t73 + t93;
t39 = t66 * t75 - t96;
t38 = t65 * t75 - t95;
t37 = -t64 * t75 + t95;
t36 = -t66 * t73 - t93;
t35 = -t51 * t73 - t52 * t75;
t34 = pkin(1) * t37 - t94;
t33 = pkin(1) * t36 - t97;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t56, -t55, 0, -t79, 0, 0, 0, 0, 0, 0, t39 * t74 + t52 * t76, t41 * t74 + t51 * t76, t54 * t74 + t58 * t76, t32 * t74 - t49 * t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t36, t37, 0, t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, -t55, -t56, 0, t83, 0, 0, 0, 0, 0, 0, t39 * t76 - t52 * t74, t41 * t76 - t51 * t74, t54 * t76 - t58 * t74, t32 * t76 + t49 * t74; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t55, 0, t56, 0, t100, -t101, t83, 0, t43 * t74 - t81, t35 * t74 + t57 * t76, t40 * t74 - t73 * t87, t42 * t74 + t81, t38 * t74 - t75 * t87, qJDD(2) * t76 + t46 * t74, t33 * t74 - t44 * t76, t34 * t74 - t45 * t76, t74 * t80, t74 * t102; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, 0, 0, 0, qJDD(1), -t62, t61, 0, 0, (-t50 + t84) * t73, -t51 * t75 + t52 * t73, -t63 * t75 - t96, (-t53 - t85) * t75, -t65 * t73 - t92, 0, -pkin(1) * t39 - t94, -pkin(1) * t41 + t97, -pkin(1) * t54 - t32, -pkin(1) * t32; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, t56, 0, -t55, 0, -t101, -t100, t79, 0, t43 * t76 + t82, t35 * t76 - t57 * t74, t40 * t76 + t73 * t89, t42 * t76 - t82, t38 * t76 + t74 * t88, -qJDD(2) * t74 + t46 * t76, t33 * t76 + t44 * t74, t34 * t76 + t45 * t74, t76 * t80, t76 * t102;];
tauB_reg = t1;
