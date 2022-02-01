% Calculate Gravitation load on the joints for
% S5RPPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRP1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:13
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPRP1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5RPPRP1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:12:46
% EndTime: 2022-01-23 09:12:46
% DurationCPUTime: 0.17s
% Computational Cost: add. (142->46), mult. (167->67), div. (0->0), fcn. (151->8), ass. (0->26)
t57 = sin(pkin(8));
t58 = cos(pkin(8));
t62 = cos(qJ(4));
t82 = -t57 * (-qJ(5) - pkin(6)) + (t62 * pkin(4) + pkin(3)) * t58;
t81 = MDP(13) + MDP(15);
t80 = MDP(14) + MDP(16);
t76 = g(3) * t57;
t56 = qJ(1) + pkin(7);
t54 = cos(t56);
t60 = sin(qJ(4));
t74 = t54 * t60;
t72 = t58 * t60;
t71 = t58 * t62;
t70 = MDP(18) + MDP(7);
t53 = sin(t56);
t63 = cos(qJ(1));
t68 = t63 * pkin(1) + t54 * pkin(2) + t53 * qJ(3);
t61 = sin(qJ(1));
t67 = -t61 * pkin(1) + t54 * qJ(3);
t66 = -g(1) * t54 - g(2) * t53;
t65 = g(1) * t53 - g(2) * t54;
t46 = t53 * t62 - t54 * t72;
t44 = t53 * t72 + t54 * t62;
t47 = t53 * t60 + t54 * t71;
t45 = -t53 * t71 + t74;
t1 = [(g(1) * t63 + g(2) * t61) * MDP(3) + t66 * MDP(6) + (-g(1) * (-t53 * pkin(2) + t67) - g(2) * t68) * MDP(7) + (-g(1) * (pkin(4) * t74 + t67) - g(2) * (t82 * t54 + t68) + (-g(1) * (-pkin(2) - t82) - g(2) * pkin(4) * t60) * t53) * MDP(18) + (t57 * MDP(17) + t58 * MDP(5)) * t65 + (pkin(1) * MDP(4) + MDP(2)) * (g(1) * t61 - g(2) * t63) + t81 * (-g(1) * t45 - g(2) * t47) + t80 * (-g(1) * t44 - g(2) * t46); (-MDP(4) - t70) * g(3); -t70 * t65; t80 * (g(1) * t47 - g(2) * t45 + t62 * t76) + (pkin(4) * MDP(18) + t81) * (-g(1) * t46 + g(2) * t44 + t60 * t76); (g(3) * t58 + t66 * t57) * MDP(18);];
taug = t1;
