% Calculate Gravitation load on the joints for
% S5RRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPP1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 22:15
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRPP1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S5RRRPP1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 22:14:49
% EndTime: 2021-01-15 22:14:50
% DurationCPUTime: 0.24s
% Computational Cost: add. (255->54), mult. (222->62), div. (0->0), fcn. (171->8), ass. (0->29)
t71 = qJ(3) + pkin(8);
t65 = sin(t71);
t66 = cos(t71);
t94 = t66 * pkin(4) + t65 * qJ(5);
t93 = MDP(14) + MDP(18);
t92 = MDP(15) - MDP(20);
t72 = qJ(1) + qJ(2);
t67 = sin(t72);
t68 = cos(t72);
t59 = g(1) * t68 + g(2) * t67;
t75 = sin(qJ(1));
t88 = t75 * pkin(1);
t73 = -qJ(4) - pkin(7);
t87 = t68 * t73;
t76 = cos(qJ(3));
t69 = t76 * pkin(3);
t64 = t69 + pkin(2);
t62 = t68 * t64;
t85 = t94 * t68 + t62;
t84 = -t67 * t73 + t62;
t58 = g(1) * t67 - g(2) * t68;
t82 = -t67 * t64 - t87;
t74 = sin(qJ(3));
t79 = (-MDP(16) - MDP(19) + MDP(6)) * t59 + (t76 * MDP(12) - t74 * MDP(13) - t92 * t65 + t93 * t66 + MDP(5)) * t58;
t78 = (-g(1) * (-t64 - t94) + g(2) * t73) * t67;
t77 = cos(qJ(1));
t70 = t77 * pkin(1);
t44 = -g(3) * t66 + t59 * t65;
t1 = [t79 + (-g(1) * (-t87 - t88) - g(2) * (t70 + t85) + t78) * MDP(21) + (g(1) * t77 + g(2) * t75) * MDP(3) + (-g(1) * (t82 - t88) - g(2) * (t70 + t84)) * MDP(17) + (g(1) * t75 - g(2) * t77) * MDP(2); (-g(1) * t82 - g(2) * t84) * MDP(17) + (g(1) * t87 - g(2) * t85 + t78) * MDP(21) + t79; (g(3) * t74 + t59 * t76) * MDP(13) + (-g(3) * (t69 + t94) + t59 * (pkin(3) * t74 + pkin(4) * t65 - qJ(5) * t66)) * MDP(21) + (pkin(3) * MDP(17) + MDP(12)) * (-g(3) * t76 + t59 * t74) + t92 * (g(3) * t65 + t59 * t66) + t93 * t44; (-MDP(17) - MDP(21)) * t58; -t44 * MDP(21);];
taug = t1;
