% Calculate Gravitation load on the joints for
% S5RRPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR12_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRR12_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR12_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR12_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRPRR12_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:30:37
% EndTime: 2019-12-31 20:30:38
% DurationCPUTime: 0.25s
% Computational Cost: add. (122->48), mult. (302->75), div. (0->0), fcn. (315->8), ass. (0->27)
t58 = sin(qJ(1));
t56 = sin(qJ(4));
t57 = sin(qJ(2));
t60 = cos(qJ(2));
t74 = cos(qJ(4));
t82 = -t60 * t56 + t57 * t74;
t44 = t82 * t58;
t50 = t57 * t56 + t60 * t74;
t45 = t50 * t58;
t61 = cos(qJ(1));
t46 = t82 * t61;
t47 = t50 * t61;
t55 = sin(qJ(5));
t59 = cos(qJ(5));
t75 = g(3) * t82;
t83 = (MDP(27) * t59 - MDP(28) * t55 + MDP(20)) * (-g(1) * t46 - g(2) * t44 + g(3) * t50) + (g(1) * t47 + g(2) * t45 + t75) * MDP(21);
t52 = g(1) * t61 + g(2) * t58;
t80 = MDP(9) + MDP(11);
t79 = MDP(10) - MDP(13);
t68 = pkin(2) * t60 + qJ(3) * t57;
t66 = t45 * t59 + t55 * t61;
t65 = t45 * t55 - t59 * t61;
t64 = pkin(1) + t68;
t42 = -g(3) * t60 + t52 * t57;
t41 = t47 * t59 - t55 * t58;
t40 = -t47 * t55 - t58 * t59;
t1 = [((-g(1) * pkin(6) - g(2) * t64) * t61 + (-g(2) * pkin(6) + g(1) * t64) * t58) * MDP(14) + (g(1) * t45 - g(2) * t47) * MDP(20) + (g(1) * t44 - g(2) * t46) * MDP(21) + (g(1) * t66 - g(2) * t41) * MDP(27) + (-g(1) * t65 - g(2) * t40) * MDP(28) + (MDP(3) - MDP(12)) * t52 + (-t79 * t57 + t80 * t60 + MDP(2)) * (g(1) * t58 - g(2) * t61); (-g(3) * t68 + t52 * (pkin(2) * t57 - qJ(3) * t60)) * MDP(14) + t79 * (g(3) * t57 + t52 * t60) + t80 * t42 - t83; -t42 * MDP(14); t83; (-g(1) * t40 + g(2) * t65 + t55 * t75) * MDP(27) + (g(1) * t41 + g(2) * t66 + t59 * t75) * MDP(28);];
taug = t1;
