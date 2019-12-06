% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5PPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
% MDP [16x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PPRPR1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PPRPR1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [16 1]), ...
  'S5PPRPR1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [16x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:01:22
% EndTime: 2019-12-05 15:01:25
% DurationCPUTime: 0.63s
% Computational Cost: add. (296->77), mult. (831->127), div. (0->0), fcn. (626->8), ass. (0->42)
t109 = sin(pkin(8));
t111 = cos(pkin(8));
t113 = sin(qJ(3));
t115 = cos(qJ(3));
t141 = -t109 * t113 + t111 * t115;
t142 = t141 * qJD(1);
t95 = t109 * t115 + t111 * t113;
t91 = t95 * qJD(3);
t108 = sin(pkin(9));
t110 = cos(pkin(9));
t140 = MDP(6) * t110 - MDP(7) * t108 + MDP(4);
t128 = t108 ^ 2 + t110 ^ 2;
t112 = sin(qJ(5));
t114 = cos(qJ(5));
t92 = t108 * t112 - t110 * t114;
t88 = t92 * qJD(5);
t94 = t108 * t114 + t110 * t112;
t86 = t94 * qJD(3);
t117 = qJD(4) - t142;
t138 = t128 * qJD(3);
t137 = pkin(6) + qJ(4);
t82 = qJD(1) * t91;
t126 = qJD(3) * t112;
t125 = qJD(3) * t114;
t103 = -pkin(4) * t110 - pkin(3);
t122 = t108 * t126;
t121 = t110 * t125;
t77 = (qJD(4) + t142) * qJD(3);
t120 = t128 * t77;
t87 = t95 * qJD(1);
t118 = t128 * (qJD(3) * qJ(4) + t87);
t89 = t94 * qJD(5);
t98 = qJD(5) * t121;
t97 = t137 * t110;
t96 = t137 * t108;
t90 = t141 * qJD(3);
t83 = -t121 + t122;
t81 = qJD(3) * t89;
t80 = -qJD(5) * t122 + t98;
t78 = -qJD(3) * pkin(3) + t117;
t75 = t103 * qJD(3) + t117;
t1 = [t90 * MDP(8) * t138 + (t118 * t90 + t95 * t120 - t141 * t82 + t78 * t91) * MDP(9) + (-t141 * t81 + t91 * t83 + (t95 * t88 - t94 * t90) * qJD(5)) * MDP(15) + (-t141 * t80 + t91 * t86 + (t95 * t89 + t92 * t90) * qJD(5)) * MDP(16) + (-t90 * MDP(5) - t140 * t91) * qJD(3); (-MDP(15) * t89 + MDP(16) * t88) * qJD(5); (t117 * t138 + t120) * MDP(8) + (-pkin(3) * t82 + qJ(4) * t120 + t117 * t118 - t78 * t87) * MDP(9) + (t80 * t94 - t86 * t88) * MDP(10) + (-t80 * t92 - t81 * t94 + t83 * t88 - t86 * t89) * MDP(11) - t88 * qJD(5) * MDP(12) - t89 * qJD(5) * MDP(13) + (t103 * t81 + t75 * t89 + t82 * t92 - t87 * t83 + ((t112 * t96 - t114 * t97) * qJD(5) - t117 * t94) * qJD(5)) * MDP(15) + (t103 * t80 - t75 * t88 + t82 * t94 - t87 * t86 + ((t112 * t97 + t114 * t96) * qJD(5) + t117 * t92) * qJD(5)) * MDP(16) + t140 * (qJD(3) * t87 - t82); (-t118 * qJD(3) + t82) * MDP(9) + t98 * MDP(16) - t128 * MDP(8) * qJD(3) ^ 2 + ((t108 * t125 + t110 * t126 + t86) * MDP(15) + (-t83 - t122) * MDP(16)) * qJD(5); t86 * t83 * MDP(10) + (-t83 ^ 2 + t86 ^ 2) * MDP(11) + (t98 + (t83 - t122) * qJD(5)) * MDP(12) + (-t75 * t86 - t94 * t77) * MDP(15) + (t75 * t83 + t92 * t77) * MDP(16);];
tauc = t1;
