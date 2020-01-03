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
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPP1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRPP1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5RRRPP1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:49:58
% EndTime: 2019-12-31 20:49:59
% DurationCPUTime: 0.17s
% Computational Cost: add. (215->52), mult. (190->62), div. (0->0), fcn. (145->8), ass. (0->27)
t62 = qJ(3) + pkin(8);
t56 = sin(t62);
t57 = cos(t62);
t83 = t57 * pkin(4) + t56 * qJ(5);
t63 = qJ(1) + qJ(2);
t58 = sin(t63);
t59 = cos(t63);
t50 = g(1) * t59 + g(2) * t58;
t66 = sin(qJ(1));
t79 = t66 * pkin(1);
t64 = -qJ(4) - pkin(7);
t78 = t59 * t64;
t67 = cos(qJ(3));
t60 = t67 * pkin(3);
t55 = t60 + pkin(2);
t53 = t59 * t55;
t76 = t83 * t59 + t53;
t75 = -t58 * t64 + t53;
t49 = g(1) * t58 - g(2) * t59;
t73 = -t58 * t55 - t78;
t65 = sin(qJ(3));
t71 = (-MDP(14) - MDP(17) + MDP(6)) * t50 + (t67 * MDP(12) - t65 * MDP(13) + MDP(16) * t57 + t56 * MDP(18) + MDP(5)) * t49;
t69 = (-g(1) * (-t55 - t83) + g(2) * t64) * t58;
t68 = cos(qJ(1));
t61 = t68 * pkin(1);
t40 = -g(3) * t57 + t50 * t56;
t1 = [(g(1) * t66 - g(2) * t68) * MDP(2) + (g(1) * t68 + g(2) * t66) * MDP(3) + (-g(1) * (t73 - t79) - g(2) * (t61 + t75)) * MDP(15) + (-g(1) * (-t78 - t79) - g(2) * (t61 + t76) + t69) * MDP(19) + t71; (-g(1) * t73 - g(2) * t75) * MDP(15) + (g(1) * t78 - g(2) * t76 + t69) * MDP(19) + t71; (g(3) * t65 + t50 * t67) * MDP(13) + t40 * MDP(16) + (-g(3) * t56 - t50 * t57) * MDP(18) + (-g(3) * (t60 + t83) + t50 * (pkin(3) * t65 + pkin(4) * t56 - qJ(5) * t57)) * MDP(19) + (pkin(3) * MDP(15) + MDP(12)) * (-g(3) * t67 + t50 * t65); (-MDP(15) - MDP(19)) * t49; -t40 * MDP(19);];
taug = t1;
