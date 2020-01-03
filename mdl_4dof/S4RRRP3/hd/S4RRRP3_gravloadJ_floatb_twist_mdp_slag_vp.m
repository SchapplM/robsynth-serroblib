% Calculate Gravitation load on the joints for
% S4RRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRRP3_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:14
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RRRP3_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP3_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP3_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S4RRRP3_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:14:15
% EndTime: 2019-12-31 17:14:16
% DurationCPUTime: 0.08s
% Computational Cost: add. (125->31), mult. (140->40), div. (0->0), fcn. (109->6), ass. (0->19)
t51 = sin(qJ(3));
t53 = cos(qJ(3));
t58 = t53 * pkin(3) + t51 * qJ(4);
t68 = -pkin(2) - t58;
t67 = MDP(12) + MDP(14);
t66 = MDP(13) - MDP(16);
t50 = qJ(1) + qJ(2);
t48 = sin(t50);
t49 = cos(t50);
t42 = g(1) * t49 + g(2) * t48;
t65 = g(1) * t48;
t60 = t48 * pkin(6) - t68 * t49;
t56 = (-MDP(15) + MDP(6)) * t42 + (-t66 * t51 + t67 * t53 + MDP(5)) * (-g(2) * t49 + t65);
t55 = t68 * t65;
t54 = cos(qJ(1));
t52 = sin(qJ(1));
t46 = t49 * pkin(6);
t31 = -g(3) * t53 + t42 * t51;
t1 = [(g(1) * t52 - g(2) * t54) * MDP(2) + (g(1) * t54 + g(2) * t52) * MDP(3) + (-g(1) * (-t52 * pkin(1) + t46) - g(2) * (t54 * pkin(1) + t60) - t55) * MDP(17) + t56; (-g(1) * t46 - g(2) * t60 - t55) * MDP(17) + t56; (-g(3) * t58 + t42 * (pkin(3) * t51 - qJ(4) * t53)) * MDP(17) + t66 * (g(3) * t51 + t42 * t53) + t67 * t31; -t31 * MDP(17);];
taug = t1;
