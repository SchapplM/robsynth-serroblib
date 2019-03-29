% Calculate Gravitation load on the joints for
% S5RRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRR2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-29 15:26
% Revision: 932832b1be1be80f59b7f1a581a1a8f328bdb39d (2019-03-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRRR2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(2,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5RRRRR2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [2x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S5RRRRR2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-29 15:26:31
% EndTime: 2019-03-29 15:26:31
% DurationCPUTime: 0.08s
% Computational Cost: add. (193->35), mult. (188->52), div. (0->0), fcn. (170->10), ass. (0->25)
t55 = qJ(3) + qJ(4);
t51 = sin(t55);
t73 = g(3) * t51;
t56 = qJ(1) + qJ(2);
t52 = sin(t56);
t57 = sin(qJ(5));
t71 = t52 * t57;
t60 = cos(qJ(5));
t70 = t52 * t60;
t54 = cos(t56);
t69 = t54 * t57;
t68 = t54 * t60;
t53 = cos(t55);
t66 = g(1) * t54 + g(2) * t52;
t67 = (t66 * t53 + t73) * MDP(20) + (t60 * MDP(26) - t57 * MDP(27) + MDP(19)) * (-g(3) * t53 + t66 * t51);
t43 = t53 * t71 + t68;
t44 = -t53 * t70 + t69;
t45 = -t53 * t69 + t70;
t46 = t53 * t68 + t71;
t58 = sin(qJ(3));
t61 = cos(qJ(3));
t64 = t66 * MDP(6) + (-g(1) * t44 - g(2) * t46) * MDP(26) + (-g(1) * t43 - g(2) * t45) * MDP(27) + (t61 * MDP(12) - t58 * MDP(13) + t53 * MDP(19) - t51 * MDP(20) + MDP(5)) * (g(1) * t52 - g(2) * t54);
t62 = cos(qJ(1));
t59 = sin(qJ(1));
t1 = [(g(1) * t59 - g(2) * t62) * MDP(2) + (g(1) * t62 + g(2) * t59) * MDP(3) + t64; t64; (-g(3) * t61 + t66 * t58) * MDP(12) + (g(3) * t58 + t66 * t61) * MDP(13) + t67; t67; (-g(1) * t45 + g(2) * t43 + t57 * t73) * MDP(26) + (g(1) * t46 - g(2) * t44 + t60 * t73) * MDP(27);];
taug  = t1;
