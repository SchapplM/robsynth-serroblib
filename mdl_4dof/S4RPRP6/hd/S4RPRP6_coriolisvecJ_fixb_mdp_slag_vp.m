% Calculate Coriolis joint torque vector for
% S4RPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPRP6_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 10:27
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RPRP6_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP6_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP6_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP6_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [5x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S4RPRP6_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 10:27:42
% EndTime: 2021-01-15 10:27:44
% DurationCPUTime: 0.29s
% Computational Cost: add. (225->82), mult. (461->123), div. (0->0), fcn. (169->2), ass. (0->39)
t84 = cos(qJ(3));
t101 = qJD(3) * t84;
t75 = pkin(3) * t101 + qJD(2);
t85 = -pkin(1) - pkin(5);
t74 = t85 * qJD(1) + qJD(2);
t95 = qJ(4) * qJD(1);
t66 = (t74 - t95) * t84;
t62 = qJD(3) * pkin(3) + t66;
t83 = sin(qJ(3));
t81 = t83 ^ 2;
t82 = t84 ^ 2;
t108 = (t81 - t82) * MDP(8);
t107 = qJD(1) * t75;
t103 = qJ(4) - t85;
t102 = qJD(3) * t83;
t100 = t83 * qJD(4);
t99 = t84 * qJD(4);
t78 = t83 * pkin(3) + qJ(2);
t69 = qJD(1) * t78 + qJD(4);
t98 = qJD(4) + t69;
t97 = qJ(2) * MDP(12);
t96 = qJ(2) * MDP(13);
t94 = qJD(3) * MDP(15);
t92 = t83 * t95;
t71 = t103 * t84;
t91 = -qJ(2) * MDP(6) - MDP(5);
t76 = qJD(3) * t92;
t60 = -qJD(1) * t99 - t74 * t102 + t76;
t61 = t74 * t101 + (-qJ(4) * t101 - t100) * qJD(1);
t90 = -t60 * t84 - t61 * t83;
t65 = t83 * t74 - t92;
t89 = t62 * t83 - t65 * t84;
t88 = t83 * MDP(14) + t84 * MDP(15);
t87 = qJD(1) ^ 2;
t86 = qJD(3) ^ 2;
t70 = t103 * t83;
t64 = -qJD(3) * t71 - t100;
t63 = t103 * t102 - t99;
t1 = [t90 * MDP(16) + (-t60 * t71 - t61 * t70 + t62 * t63 + t65 * t64 + t69 * t75) * MDP(17) + (t78 * MDP(17) + t88) * t107 + ((t69 * t84 + t63) * MDP(14) + (-t69 * t83 - t64) * MDP(15) + t89 * MDP(16)) * qJD(3) + ((-MDP(13) * t85 - MDP(10)) * t84 + (-MDP(12) * t85 - MDP(9)) * t83) * t86 + ((-t63 * t84 - t64 * t83) * MDP(16) + t88 * t75 + 0.2e1 * (MDP(12) * t83 + MDP(13) * t84 - t91) * qJD(2) + (0.2e1 * t108 + (t78 * MDP(14) + t70 * MDP(16) + 0.2e1 * t97) * t84 + (-t78 * MDP(15) - t71 * MDP(16) - 0.2e1 * t84 * MDP(7) - 0.2e1 * t96) * t83) * qJD(3)) * qJD(1); (-t69 * qJD(1) - t89 * qJD(3) - t90) * MDP(17) + t91 * t87 + ((MDP(13) + MDP(15)) * t84 + (MDP(12) + MDP(14)) * t83) * (-t86 - t87); (t65 * qJD(3) + t76) * MDP(14) + t66 * t94 + (t60 * pkin(3) + (t62 - t66) * t65) * MDP(17) + (-pkin(3) * t82 * MDP(15) - t108) * t87 + (-t87 * t97 - t74 * t94 + (-pkin(3) * t69 * MDP(17) - t98 * MDP(14) + qJ(4) * t94) * qJD(1)) * t84 + (-qJD(3) * t74 * MDP(14) + (t96 + (-MDP(14) * pkin(3) + MDP(7)) * t84) * t87 + t98 * MDP(15) * qJD(1)) * t83; t107 * MDP(17) + (-t81 - t82) * MDP(16) * t87 + ((t62 * t84 + t65 * t83) * MDP(17) + 0.2e1 * (t84 * MDP(14) - t83 * MDP(15)) * qJD(3)) * qJD(1);];
tauc = t1;
