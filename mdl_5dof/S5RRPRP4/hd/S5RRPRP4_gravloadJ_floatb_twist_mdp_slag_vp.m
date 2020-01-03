% Calculate Gravitation load on the joints for
% S5RRPRP4
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
% MDP [20x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRP4_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:53
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRP4_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP4_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP4_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [20 1]), ...
  'S5RRPRP4_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [20x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:52:56
% EndTime: 2019-12-31 19:52:57
% DurationCPUTime: 0.14s
% Computational Cost: add. (174->43), mult. (178->51), div. (0->0), fcn. (135->6), ass. (0->21)
t64 = sin(qJ(4));
t66 = cos(qJ(4));
t82 = -t64 * pkin(4) + t66 * qJ(5);
t81 = -MDP(15) - MDP(17);
t80 = MDP(16) - MDP(19);
t63 = qJ(1) + qJ(2);
t59 = sin(t63);
t60 = cos(t63);
t79 = -g(1) * t59 + g(2) * t60;
t65 = sin(qJ(1));
t75 = t65 * pkin(1);
t74 = t60 * pkin(2) + t59 * qJ(3);
t55 = t60 * qJ(3);
t72 = -pkin(2) * t59 + t55;
t70 = t60 * pkin(7) - t59 * t82 + t74;
t69 = -(MDP(18) + MDP(5) - MDP(7)) * t79 + (t64 * t81 - t66 * t80 + MDP(6) - MDP(8)) * (g(1) * t60 + g(2) * t59);
t68 = t55 + (-pkin(2) - pkin(7)) * t59 - t82 * t60;
t67 = cos(qJ(1));
t61 = t67 * pkin(1);
t37 = -g(3) * t64 - t66 * t79;
t1 = [(g(1) * t65 - g(2) * t67) * MDP(2) + (g(1) * t67 + g(2) * t65) * MDP(3) + (-g(1) * (t72 - t75) - g(2) * (t61 + t74)) * MDP(9) + (-g(1) * (t68 - t75) - g(2) * (t61 + t70)) * MDP(20) + t69; (-g(1) * t72 - g(2) * t74) * MDP(9) + (-g(1) * t68 - g(2) * t70) * MDP(20) + t69; -(-MDP(20) - MDP(9)) * t79; (-g(3) * t82 + t79 * (pkin(4) * t66 + qJ(5) * t64)) * MDP(20) + t81 * t37 + t80 * (g(3) * t66 - t64 * t79); t37 * MDP(20);];
taug = t1;
