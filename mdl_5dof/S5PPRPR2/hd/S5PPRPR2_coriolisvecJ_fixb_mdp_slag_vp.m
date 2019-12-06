% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5PPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2]';
% MDP [15x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PPRPR2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:03
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PPRPR2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(15,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'S5PPRPR2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:03:22
% EndTime: 2019-12-05 15:03:24
% DurationCPUTime: 0.27s
% Computational Cost: add. (151->56), mult. (404->85), div. (0->0), fcn. (262->6), ass. (0->40)
t91 = MDP(4) - MDP(6);
t75 = sin(qJ(5));
t77 = cos(qJ(5));
t108 = t75 * MDP(14) + t77 * MDP(15);
t73 = sin(pkin(8));
t74 = cos(pkin(8));
t76 = sin(qJ(3));
t78 = cos(qJ(3));
t64 = t73 * t78 + t74 * t76;
t60 = t64 * qJD(1);
t57 = qJD(3) * qJ(4) + t60;
t89 = qJD(1) * qJD(3);
t86 = t78 * t89;
t87 = t76 * t89;
t58 = t73 * t86 + t74 * t87;
t107 = -t57 * qJD(3) + t58;
t106 = (t75 ^ 2 - t77 ^ 2) * MDP(10);
t105 = MDP(5) - MDP(7);
t101 = t78 * t74;
t102 = t76 * t73;
t63 = -t101 + t102;
t104 = qJD(5) * (t57 - t60);
t80 = qJD(5) ^ 2;
t103 = t63 * t80;
t97 = t77 * MDP(9);
t62 = t64 * qJD(3);
t96 = qJD(5) * t62;
t93 = t77 * MDP(14);
t88 = qJD(1) * t102;
t59 = -qJD(1) * t101 + t88;
t90 = qJD(4) + t59;
t67 = t74 * t86;
t70 = qJD(3) * qJD(4);
t55 = t73 * t87 - t67 - t70;
t85 = -(-pkin(3) - pkin(6)) * t80 - t55;
t82 = -MDP(15) * t75 + t93;
t81 = qJD(3) ^ 2;
t61 = t63 * qJD(3);
t56 = -qJD(3) * pkin(3) + t90;
t1 = [(-t55 * t64 + t56 * t62 - t57 * t61 + t58 * t63) * MDP(8) + (-t103 * t75 + t77 * t96) * MDP(14) + (-t103 * t77 - t75 * t96) * MDP(15) + (-t91 * t62 + t82 * t64 * qJD(5) + (t105 - t108) * t61) * qJD(3); -t82 * t80; -t67 * MDP(5) + (t67 + 0.2e1 * t70) * MDP(7) + (-t55 * qJ(4) - t56 * t60 + t57 * t90) * MDP(8) + (-t80 * MDP(12) + MDP(14) * t104 + MDP(15) * t85) * t77 + (-t80 * MDP(11) + MDP(14) * t85 - MDP(15) * t104) * t75 + (t91 * t60 + (qJ(4) * t82 - 0.2e1 * t75 * t97 + 0.2e1 * t106) * qJD(5) + t108 * t90 + t105 * (-t59 + t88)) * qJD(3) + (-MDP(8) * pkin(3) - t91) * t58; -t81 * MDP(7) + t107 * MDP(8) + t108 * (-t80 - t81); -t81 * t106 + t107 * t93 + (-MDP(15) * t107 + t81 * t97) * t75;];
tauc = t1;
