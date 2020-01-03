% Calculate Gravitation load on the joints for
% S5RRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPP3_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRPP3_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP3_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP3_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S5RRRPP3_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:53:42
% EndTime: 2019-12-31 20:53:43
% DurationCPUTime: 0.25s
% Computational Cost: add. (233->55), mult. (253->63), div. (0->0), fcn. (202->6), ass. (0->34)
t80 = sin(qJ(3));
t75 = t80 * qJ(4);
t82 = cos(qJ(3));
t76 = t82 * pkin(3);
t94 = t76 + t75;
t102 = MDP(12) - MDP(15) + MDP(20);
t101 = MDP(13) - MDP(16) - MDP(19);
t79 = qJ(1) + qJ(2);
t73 = sin(t79);
t74 = cos(t79);
t62 = g(1) * t74 + g(2) * t73;
t100 = t62 * t80;
t99 = pkin(3) * t80;
t98 = g(1) * t73;
t95 = pkin(3) + qJ(5);
t93 = qJ(4) * t82;
t92 = t82 * qJ(5);
t70 = t74 * pkin(7);
t81 = sin(qJ(1));
t91 = -t81 * pkin(1) + t70;
t90 = -pkin(2) - t75;
t89 = t73 * pkin(7) + (pkin(2) + t94) * t74;
t87 = t73 * pkin(4) + t74 * t92 + t89;
t86 = (t90 - t76) * t98;
t85 = (-MDP(14) - MDP(18) + MDP(6)) * t62 + (-t101 * t80 + t102 * t82 + MDP(5)) * (-g(2) * t74 + t98);
t84 = (-t95 * t82 + t90) * t98;
t83 = cos(qJ(1));
t77 = t83 * pkin(1);
t71 = t74 * pkin(4);
t66 = t74 * t93;
t63 = t73 * t93;
t49 = g(3) * t80 + t62 * t82;
t48 = -g(3) * t82 + t100;
t1 = [(-g(1) * (t71 + t91) - g(2) * (t77 + t87) - t84) * MDP(21) + t85 + (-g(1) * t91 - g(2) * (t77 + t89) - t86) * MDP(17) + (g(1) * t81 - g(2) * t83) * MDP(2) + (g(1) * t83 + g(2) * t81) * MDP(3); (-g(1) * t70 - g(2) * t89 - t86) * MDP(17) + (-g(1) * (t70 + t71) - g(2) * t87 - t84) * MDP(21) + t85; (-g(1) * (-t74 * t99 + t66) - g(2) * (-t73 * t99 + t63) - g(3) * t94) * MDP(17) + (-g(1) * t66 - g(2) * t63 - g(3) * (t92 + t94) + t95 * t100) * MDP(21) + t101 * t49 + t102 * t48; (-MDP(17) - MDP(21)) * t48; -t49 * MDP(21);];
taug = t1;
