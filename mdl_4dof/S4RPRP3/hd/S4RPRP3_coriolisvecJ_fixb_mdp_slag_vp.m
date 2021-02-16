% Calculate Coriolis joint torque vector for
% S4RPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
% MDP [15x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPRP3_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 10:20
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RPRP3_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(15,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP3_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP3_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP3_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'S4RPRP3_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 10:20:37
% EndTime: 2021-01-15 10:20:38
% DurationCPUTime: 0.32s
% Computational Cost: add. (232->73), mult. (551->114), div. (0->0), fcn. (252->4), ass. (0->43)
t65 = sin(pkin(6)) * pkin(1) + pkin(5);
t63 = t65 * qJD(1);
t72 = sin(qJ(3));
t86 = t72 * qJD(3);
t61 = t63 * t86;
t68 = t72 ^ 2;
t73 = cos(qJ(3));
t69 = t73 ^ 2;
t98 = (t68 - t69) * MDP(6);
t78 = qJ(4) * qJD(1) + t63;
t55 = t72 * qJD(2) + t78 * t73;
t97 = pkin(3) * t68;
t75 = qJD(1) ^ 2;
t96 = t73 * t75;
t67 = t73 * qJD(2);
t54 = -t78 * t72 + t67;
t93 = qJD(3) * pkin(3);
t53 = t54 + t93;
t95 = t53 - t54;
t92 = qJ(4) + t65;
t66 = -cos(pkin(6)) * pkin(1) - pkin(2);
t91 = qJD(1) * (-t73 * pkin(3) + t66);
t64 = qJD(1) * t66;
t90 = qJD(1) * t72;
t89 = qJD(1) * t73;
t88 = t64 * qJD(1);
t85 = t72 * qJD(4);
t84 = t73 * qJD(4);
t58 = qJD(4) + t91;
t83 = -qJD(4) - t58;
t82 = qJ(4) * t86;
t81 = t58 + t91;
t80 = 0.2e1 * t64;
t79 = qJD(3) * t92;
t77 = t53 * t72 - t55 * t73;
t74 = qJD(3) ^ 2;
t60 = t92 * t73;
t59 = t92 * t72;
t57 = -t73 * t79 - t85;
t56 = -t72 * t79 + t84;
t52 = -qJD(1) * t85 - t55 * qJD(3);
t51 = t67 * qJD(3) - t61 + (-t82 + t84) * qJD(1);
t1 = [(t51 * t73 - t52 * t72 + t56 * t89 - t57 * t90) * MDP(14) + (t51 * t60 - t52 * t59 + t53 * t57 + t55 * t56) * MDP(15) + (t73 * MDP(7) - t72 * MDP(8) + (-MDP(10) * t73 + MDP(11) * t72) * t65) * t74 + (t57 * MDP(12) - t56 * MDP(13) + 0.2e1 * (MDP(13) * t97 - t98) * qJD(1) + (t80 * MDP(11) + t81 * MDP(13) + (qJD(1) * t59 - t53) * MDP(14)) * t73 + (0.2e1 * MDP(5) * t89 + t80 * MDP(10) + t81 * MDP(12) + (-qJD(1) * t60 - t55) * MDP(14) + (-0.2e1 * MDP(12) * t89 + t81 * MDP(15)) * pkin(3)) * t72) * qJD(3); (-t77 * qJD(3) + t51 * t72 + t52 * t73) * MDP(15) + ((-MDP(11) - MDP(13)) * t73 + (-MDP(10) - MDP(12)) * t72) * t74; t75 * t98 - t73 * t88 * MDP(11) + (-t75 * t97 + t61 + (t54 - t67) * qJD(3) + (t83 * t73 + t82) * qJD(1)) * MDP(13) + (-t93 + t95) * MDP(14) * t89 + (t95 * t55 + (-t58 * t90 + t52) * pkin(3)) * MDP(15) + (-MDP(5) * t96 - MDP(10) * t88 + (pkin(3) * t96 + t83 * qJD(1)) * MDP(12)) * t72; (-t68 - t69) * MDP(14) * t75 + (t77 * MDP(15) + (0.2e1 * t73 * MDP(13) + (MDP(15) * pkin(3) + 0.2e1 * MDP(12)) * t72) * qJD(3)) * qJD(1);];
tauc = t1;
