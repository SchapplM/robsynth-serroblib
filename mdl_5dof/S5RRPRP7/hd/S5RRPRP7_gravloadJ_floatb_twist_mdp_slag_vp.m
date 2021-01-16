% Calculate Gravitation load on the joints for
% S5RRPRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRP7_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 20:42
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRP7_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP7_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP7_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RRPRP7_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 20:41:25
% EndTime: 2021-01-15 20:41:28
% DurationCPUTime: 0.40s
% Computational Cost: add. (215->70), mult. (311->92), div. (0->0), fcn. (291->10), ass. (0->36)
t102 = MDP(12) - MDP(23);
t100 = MDP(20) + MDP(22);
t99 = MDP(21) - MDP(24);
t80 = sin(qJ(1));
t83 = cos(qJ(1));
t101 = -g(1) * t83 - g(2) * t80;
t74 = qJ(2) + pkin(8);
t70 = sin(t74);
t96 = g(3) * t70;
t77 = -qJ(3) - pkin(6);
t95 = t80 * t77;
t78 = sin(qJ(4));
t94 = t80 * t78;
t81 = cos(qJ(4));
t93 = t80 * t81;
t92 = t83 * t78;
t91 = t83 * t81;
t64 = g(1) * t80 - g(2) * t83;
t88 = -pkin(4) * t81 - qJ(5) * t78;
t71 = cos(t74);
t58 = t71 * t94 + t91;
t60 = t71 * t92 - t93;
t50 = g(1) * t60 + g(2) * t58 + t78 * t96;
t82 = cos(qJ(2));
t79 = sin(qJ(2));
t76 = cos(pkin(8));
t75 = sin(pkin(8));
t72 = t82 * pkin(2);
t69 = t72 + pkin(1);
t66 = t77 * t83;
t63 = -t75 * pkin(3) + t76 * pkin(7);
t62 = t76 * pkin(3) + t75 * pkin(7) + pkin(2);
t61 = t71 * t91 + t94;
t59 = t71 * t93 - t92;
t54 = t62 * t82 + t63 * t79 + pkin(1);
t1 = [(-g(1) * (-t80 * t69 - t66) - g(2) * (t83 * t69 - t95)) * MDP(14) + (-g(1) * (-t59 * pkin(4) - t58 * qJ(5) - t54 * t80 - t66) - g(2) * (t61 * pkin(4) + t60 * qJ(5) + t54 * t83 - t95)) * MDP(25) - (MDP(3) - MDP(13)) * t101 + t100 * (g(1) * t59 - g(2) * t61) - t99 * (g(1) * t58 - g(2) * t60) + (-t79 * MDP(10) + t71 * MDP(11) + t82 * MDP(9) - t102 * t70 + MDP(2)) * t64; (g(3) * t79 - t101 * t82) * MDP(10) + (-g(3) * (t70 * pkin(7) + t72 + (pkin(3) - t88) * t71) + t101 * (-t62 * t79 + t63 * t82 + t70 * t88)) * MDP(25) + (MDP(14) * pkin(2) + MDP(9)) * (-g(3) * t82 - t101 * t79) + t102 * (-t101 * t71 + t96) + (t100 * t81 - t99 * t78 + MDP(11)) * (-g(3) * t71 - t101 * t70); (-MDP(14) - MDP(25)) * t64; (-g(1) * (-t60 * pkin(4) + t61 * qJ(5)) - g(2) * (-t58 * pkin(4) + t59 * qJ(5)) - (-pkin(4) * t78 + qJ(5) * t81) * t96) * MDP(25) + t99 * (g(1) * t61 + g(2) * t59 + t81 * t96) + t100 * t50; -t50 * MDP(25);];
taug = t1;
