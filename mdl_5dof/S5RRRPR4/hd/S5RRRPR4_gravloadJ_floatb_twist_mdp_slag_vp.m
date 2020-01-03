% Calculate Gravitation load on the joints for
% S5RRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPR4_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRPR4_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR4_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR4_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S5RRRPR4_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:11:40
% EndTime: 2019-12-31 21:11:40
% DurationCPUTime: 0.13s
% Computational Cost: add. (196->43), mult. (240->62), div. (0->0), fcn. (221->8), ass. (0->28)
t77 = sin(qJ(3));
t80 = cos(qJ(3));
t88 = t80 * pkin(3) + t77 * qJ(4);
t101 = -pkin(2) - t88;
t100 = MDP(12) + MDP(14);
t99 = MDP(13) - MDP(16);
t76 = sin(qJ(5));
t79 = cos(qJ(5));
t86 = t80 * t76 - t77 * t79;
t75 = qJ(1) + qJ(2);
t73 = sin(t75);
t74 = cos(t75);
t64 = g(1) * t74 + g(2) * t73;
t51 = t86 * t73;
t65 = t77 * t76 + t80 * t79;
t52 = t65 * t73;
t53 = t86 * t74;
t54 = t65 * t74;
t98 = (g(1) * t53 + g(2) * t51 + g(3) * t65) * MDP(23) + (g(1) * t54 + g(2) * t52 - g(3) * t86) * MDP(24);
t97 = g(1) * t73;
t90 = t73 * pkin(7) - t101 * t74;
t83 = t101 * t97;
t82 = (-g(1) * t51 + g(2) * t53) * MDP(24) + (g(1) * t52 - g(2) * t54) * MDP(23) + (MDP(6) - MDP(15)) * t64 + (t100 * t80 - t99 * t77 + MDP(5)) * (-g(2) * t74 + t97);
t81 = cos(qJ(1));
t78 = sin(qJ(1));
t71 = t74 * pkin(7);
t49 = -g(3) * t80 + t64 * t77;
t1 = [(g(1) * t78 - g(2) * t81) * MDP(2) + (g(1) * t81 + g(2) * t78) * MDP(3) + (-g(1) * (-t78 * pkin(1) + t71) - g(2) * (t81 * pkin(1) + t90) - t83) * MDP(17) + t82; (-g(1) * t71 - g(2) * t90 - t83) * MDP(17) + t82; (-g(3) * t88 + t64 * (pkin(3) * t77 - qJ(4) * t80)) * MDP(17) + t99 * (g(3) * t77 + t64 * t80) + t100 * t49 - t98; -t49 * MDP(17); t98;];
taug = t1;
