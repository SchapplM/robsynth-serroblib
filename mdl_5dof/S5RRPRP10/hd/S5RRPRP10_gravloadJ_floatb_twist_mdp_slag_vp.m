% Calculate Gravitation load on the joints for
% S5RRPRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRP10_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 21:03
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRP10_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP10_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP10_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RRPRP10_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 21:02:46
% EndTime: 2021-01-15 21:02:47
% DurationCPUTime: 0.20s
% Computational Cost: add. (135->53), mult. (267->71), div. (0->0), fcn. (237->6), ass. (0->30)
t96 = MDP(10) - MDP(13);
t95 = MDP(20) + MDP(22);
t94 = MDP(21) + MDP(23);
t73 = sin(qJ(2));
t76 = cos(qJ(2));
t74 = sin(qJ(1));
t77 = cos(qJ(1));
t80 = g(1) * t77 + g(2) * t74;
t53 = g(3) * t73 + t80 * t76;
t92 = MDP(9) - MDP(12) + MDP(24);
t88 = g(3) * t76;
t72 = sin(qJ(4));
t87 = t72 * t77;
t86 = t74 * t72;
t75 = cos(qJ(4));
t85 = t74 * t75;
t84 = t75 * t77;
t83 = t76 * pkin(2) + t73 * qJ(3);
t67 = pkin(4) * t72 + qJ(3);
t71 = pkin(2) + pkin(7) + qJ(5);
t81 = t67 * t73 + t71 * t76;
t56 = t73 * t84 - t86;
t58 = t73 * t85 + t87;
t66 = pkin(4) * t75 + pkin(3) + pkin(6);
t62 = pkin(1) + t83;
t59 = -t73 * t86 + t84;
t57 = t73 * t87 + t85;
t54 = pkin(1) + t81;
t52 = t80 * t73 - t88;
t1 = [(-g(1) * (pkin(6) * t77 - t62 * t74) - g(2) * (t74 * pkin(6) + t62 * t77)) * MDP(14) + (-g(1) * (-t54 * t74 + t66 * t77) - g(2) * (t54 * t77 + t66 * t74)) * MDP(25) + (MDP(3) - MDP(11)) * t80 + t95 * (-g(1) * t59 - g(2) * t57) + t94 * (g(1) * t58 - g(2) * t56) + (-t96 * t73 + t92 * t76 + MDP(2)) * (g(1) * t74 - g(2) * t77); (-g(3) * t83 - t80 * (-pkin(2) * t73 + qJ(3) * t76)) * MDP(14) + (-g(3) * t81 - t80 * (t67 * t76 - t71 * t73)) * MDP(25) + t92 * t52 + (-t95 * t72 - t94 * t75 + t96) * t53; (-MDP(14) - MDP(25)) * t52; t94 * (g(1) * t57 - g(2) * t59 - t72 * t88) + (pkin(4) * MDP(25) + t95) * (-g(1) * t56 - g(2) * t58 + t75 * t88); -t53 * MDP(25);];
taug = t1;
