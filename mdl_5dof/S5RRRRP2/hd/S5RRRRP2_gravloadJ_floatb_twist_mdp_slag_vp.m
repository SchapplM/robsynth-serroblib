% Calculate Gravitation load on the joints for
% S5RRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRP2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:49
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRRP2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S5RRRRP2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:49:28
% EndTime: 2022-01-20 11:49:28
% DurationCPUTime: 0.12s
% Computational Cost: add. (222->39), mult. (192->50), div. (0->0), fcn. (147->8), ass. (0->23)
t71 = MDP(19) + MDP(21);
t70 = MDP(20) + MDP(22);
t59 = qJ(3) + qJ(4);
t55 = cos(t59);
t63 = cos(qJ(3));
t69 = t63 * pkin(3) + pkin(4) * t55;
t49 = pkin(2) + t69;
t60 = qJ(1) + qJ(2);
t54 = sin(t60);
t56 = cos(t60);
t58 = qJ(5) + pkin(7) + pkin(8);
t68 = t56 * t49 + t54 * t58;
t67 = -t49 * t54 + t58 * t56;
t48 = g(1) * t56 + g(2) * t54;
t53 = sin(t59);
t33 = -g(3) * t55 + t48 * t53;
t66 = t70 * (g(3) * t53 + t48 * t55) + t71 * t33;
t47 = g(1) * t54 - g(2) * t56;
t61 = sin(qJ(3));
t65 = (-MDP(23) + MDP(6)) * t48 + (t63 * MDP(12) - t61 * MDP(13) - t70 * t53 + t71 * t55 + MDP(5)) * t47;
t64 = cos(qJ(1));
t62 = sin(qJ(1));
t1 = [(g(1) * t62 - g(2) * t64) * MDP(2) + (g(1) * t64 + g(2) * t62) * MDP(3) + (-g(1) * (-pkin(1) * t62 + t67) - g(2) * (pkin(1) * t64 + t68)) * MDP(24) + t65; (-g(1) * t67 - g(2) * t68) * MDP(24) + t65; (-g(3) * t63 + t48 * t61) * MDP(12) + (g(3) * t61 + t48 * t63) * MDP(13) + (-g(3) * t69 - t48 * (-pkin(3) * t61 - pkin(4) * t53)) * MDP(24) + t66; t33 * MDP(24) * pkin(4) + t66; -t47 * MDP(24);];
taug = t1;
