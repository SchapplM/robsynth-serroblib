% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5PPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1]';
% MDP [14x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PPRRP4_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PPRRP4_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP4_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP4_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRRP4_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [14 1]), ...
  'S5PPRRP4_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [14x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:34:42
% EndTime: 2019-12-31 17:34:43
% DurationCPUTime: 0.31s
% Computational Cost: add. (236->76), mult. (567->122), div. (0->0), fcn. (294->4), ass. (0->49)
t76 = sin(qJ(4));
t74 = t76 ^ 2;
t78 = cos(qJ(4));
t75 = t78 ^ 2;
t111 = (t74 - t75) * MDP(7);
t102 = qJD(4) * pkin(4);
t77 = sin(qJ(3));
t97 = t77 * qJD(2);
t68 = qJD(3) * pkin(6) + t97;
t87 = qJ(5) * qJD(3) + t68;
t95 = t78 * qJD(1);
t60 = -t87 * t76 - t95;
t59 = t60 + t102;
t73 = t76 * qJD(1);
t61 = t87 * t78 - t73;
t85 = t59 * t76 - t61 * t78;
t110 = t85 * MDP(14);
t79 = cos(qJ(3));
t94 = t79 * qJD(2);
t86 = qJD(5) + t94;
t93 = qJ(5) * qJD(4);
t57 = (-t76 * t68 - t95) * qJD(4) + (-t76 * t93 + t86 * t78) * qJD(3);
t109 = (qJD(4) * t59 - t57) * MDP(14);
t90 = (t74 + t75) * MDP(13);
t108 = -qJ(5) - pkin(6);
t107 = t59 - t60;
t98 = t76 * qJD(4);
t65 = (pkin(4) * t98 + t97) * qJD(3);
t80 = qJD(4) ^ 2;
t81 = qJD(3) ^ 2;
t104 = t80 + t81;
t103 = qJD(3) * pkin(3);
t92 = -t78 * pkin(4) - pkin(3);
t64 = t92 * qJD(3) + qJD(5) - t94;
t101 = MDP(14) * t64;
t100 = qJD(3) * t78;
t99 = t76 * qJD(3);
t96 = t78 * MDP(12);
t91 = qJD(4) * t108;
t84 = t103 * qJD(3);
t72 = qJD(1) * t98;
t58 = -qJD(4) * t78 * t68 + t72 + (-t86 * t76 - t78 * t93) * qJD(3);
t83 = (-qJD(4) * t61 - t58) * MDP(14);
t82 = -0.2e1 * t103;
t67 = t108 * t78;
t66 = t108 * t76;
t63 = -t76 * qJD(5) + t78 * t91;
t62 = t78 * qJD(5) + t76 * t91;
t1 = [(t80 * MDP(12) + t83) * t78 + (t80 * MDP(11) + t109) * t76; (-t65 * MDP(14) + (-MDP(5) + t90) * t81 + (-t110 + 0.2e1 * (-t76 * MDP(11) - t96) * qJD(4)) * qJD(3)) * t79 + (qJD(3) * t101 - t81 * MDP(4) + (-t104 * MDP(11) - t109) * t78 + (t104 * MDP(12) + t83) * t76) * t77; (t62 * t100 + t57 * t78 - t58 * t76 - t63 * t99) * MDP(13) + (-t57 * t67 + t58 * t66 + t59 * t63 + t61 * t62 + t65 * t92) * MDP(14) + (-t77 * t101 + (-qJD(3) * t90 + t110) * t79) * qJD(2) + (t78 * MDP(8) - t76 * MDP(9) + (-t78 * MDP(11) + t76 * MDP(12)) * pkin(6)) * t80 + (-0.2e1 * qJD(3) * t111 + (t82 * MDP(12) + (-qJD(3) * t66 - t59) * MDP(13)) * t78 + (0.2e1 * MDP(6) * t100 + t82 * MDP(11) + (qJD(3) * t67 - t61) * MDP(13) + pkin(4) * t101) * t76) * qJD(4); (-t73 * qJD(4) + t76 * t84 + t72) * MDP(11) + t84 * t96 + (-t102 + t107) * MDP(13) * t100 + (t107 * t61 + (-t64 * t99 + t58) * pkin(4)) * MDP(14) + (-t76 * t78 * MDP(6) + t111) * t81; (qJD(3) * t85 + t65) * MDP(14) - t81 * t90;];
tauc = t1;
