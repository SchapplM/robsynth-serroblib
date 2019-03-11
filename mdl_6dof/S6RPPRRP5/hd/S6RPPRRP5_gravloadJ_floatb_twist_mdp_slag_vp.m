% Calculate Gravitation load on the joints for
% S6RPPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRRP5_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPPRRP5_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP5_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP5_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S6RPPRRP5_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:08:47
% EndTime: 2019-03-09 02:08:48
% DurationCPUTime: 0.21s
% Computational Cost: add. (103->54), mult. (204->71), div. (0->0), fcn. (171->6), ass. (0->29)
t59 = sin(qJ(1));
t62 = cos(qJ(1));
t49 = g(1) * t62 + g(2) * t59;
t82 = MDP(16) - MDP(24);
t58 = sin(qJ(4));
t61 = cos(qJ(4));
t42 = -g(3) * t58 + t49 * t61;
t75 = g(3) * t61;
t57 = sin(qJ(5));
t74 = t59 * t57;
t60 = cos(qJ(5));
t73 = t59 * t60;
t72 = t62 * t57;
t71 = t62 * t60;
t70 = -pkin(1) - qJ(3);
t69 = t62 * pkin(1) + t59 * qJ(2);
t68 = -MDP(25) - MDP(9);
t66 = pkin(5) * t57 + pkin(7);
t65 = g(2) * (t62 * qJ(3) + t69);
t48 = g(1) * t59 - g(2) * t62;
t50 = t60 * pkin(5) + pkin(4);
t56 = -qJ(6) - pkin(8);
t63 = t58 * t50 + t61 * t56;
t45 = -t58 * t72 - t73;
t43 = t58 * t74 - t71;
t53 = t62 * qJ(2);
t46 = t58 * t71 - t74;
t44 = -t58 * t73 - t72;
t1 = [(-g(1) * (-t59 * pkin(1) + t53) - g(2) * t69) * MDP(6) + (-g(1) * (t70 * t59 + t53) - t65) * MDP(9) + (-g(1) * t44 - g(2) * t46) * MDP(22) + (-g(1) * t43 - g(2) * t45) * MDP(23) + (-g(1) * t53 - t65 + (g(1) * t66 - g(2) * t63) * t62 + (-g(1) * (-t63 + t70) + g(2) * t66) * t59) * MDP(25) + (MDP(3) - MDP(5) - MDP(7)) * t49 + (t58 * MDP(15) + t82 * t61 + MDP(2) - MDP(4) + MDP(8)) * t48; (-MDP(6) + t68) * t48; t68 * t49; (g(3) * t63 - t49 * (t50 * t61 - t56 * t58)) * MDP(25) + t82 * (t49 * t58 + t75) + (-MDP(22) * t60 + MDP(23) * t57 - MDP(15)) * t42; (g(1) * t46 - g(2) * t44 + t60 * t75) * MDP(23) + (pkin(5) * MDP(25) + MDP(22)) * (-g(1) * t45 + g(2) * t43 + t57 * t75); t42 * MDP(25);];
taug  = t1;
