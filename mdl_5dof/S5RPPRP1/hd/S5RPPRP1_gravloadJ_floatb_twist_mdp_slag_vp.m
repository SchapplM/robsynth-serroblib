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
% Datum: 2021-01-15 16:56
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
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
% StartTime: 2021-01-15 16:55:48
% EndTime: 2021-01-15 16:55:50
% DurationCPUTime: 0.22s
% Computational Cost: add. (142->45), mult. (167->67), div. (0->0), fcn. (151->8), ass. (0->26)
t58 = sin(pkin(8));
t59 = cos(pkin(8));
t63 = cos(qJ(4));
t83 = -t58 * (-qJ(5) - pkin(6)) + (t63 * pkin(4) + pkin(3)) * t59;
t82 = MDP(13) + MDP(15);
t81 = MDP(14) + MDP(16);
t79 = g(1) * t58;
t57 = qJ(1) + pkin(7);
t53 = sin(t57);
t61 = sin(qJ(4));
t75 = t53 * t61;
t73 = t59 * t61;
t72 = t59 * t63;
t62 = sin(qJ(1));
t71 = t62 * pkin(1) + t53 * pkin(2);
t70 = MDP(18) + MDP(7);
t54 = cos(t57);
t64 = cos(qJ(1));
t68 = t64 * pkin(1) + t54 * pkin(2) + t53 * qJ(3);
t67 = g(2) * t54 + g(3) * t53;
t66 = -g(2) * t53 + g(3) * t54;
t46 = -t53 * t63 + t54 * t73;
t44 = -t53 * t73 - t54 * t63;
t47 = t54 * t72 + t75;
t45 = t53 * t72 - t54 * t61;
t1 = [(g(2) * t62 - g(3) * t64) * MDP(3) + t66 * MDP(6) + (-g(2) * t68 - g(3) * (-t54 * qJ(3) + t71)) * MDP(7) + (-g(2) * (pkin(4) * t75 + t68) - g(3) * (t83 * t53 + t71) + (-g(2) * t83 - g(3) * (-pkin(4) * t61 - qJ(3))) * t54) * MDP(18) + (-t58 * MDP(17) - t59 * MDP(5)) * t67 + (MDP(4) * pkin(1) + MDP(2)) * (-g(2) * t64 - g(3) * t62) + t82 * (-g(2) * t47 - g(3) * t45) + t81 * (g(2) * t46 - g(3) * t44); (-MDP(4) - t70) * g(1); t70 * t67; t81 * (g(2) * t45 - g(3) * t47 + t63 * t79) + (pkin(4) * MDP(18) + t82) * (-g(2) * t44 - g(3) * t46 + t61 * t79); (g(1) * t59 + t58 * t66) * MDP(18);];
taug = t1;
