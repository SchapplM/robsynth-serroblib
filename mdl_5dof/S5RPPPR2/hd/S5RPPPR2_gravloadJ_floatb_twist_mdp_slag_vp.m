% Calculate Gravitation load on the joints for
% S5RPPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPPR2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPPR2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RPPPR2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:22:56
% EndTime: 2020-01-03 11:22:59
% DurationCPUTime: 0.36s
% Computational Cost: add. (122->60), mult. (280->97), div. (0->0), fcn. (306->10), ass. (0->37)
t76 = sin(pkin(7));
t78 = cos(pkin(7));
t104 = pkin(2) * t78 + qJ(3) * t76;
t82 = cos(qJ(1));
t91 = cos(pkin(8));
t88 = t82 * t91;
t75 = sin(pkin(8));
t80 = sin(qJ(1));
t95 = t80 * t75;
t62 = t78 * t88 + t95;
t74 = sin(pkin(9));
t77 = cos(pkin(9));
t96 = t76 * t82;
t55 = t62 * t77 + t74 * t96;
t89 = t80 * t91;
t94 = t82 * t75;
t61 = t78 * t94 - t89;
t79 = sin(qJ(5));
t81 = cos(qJ(5));
t103 = t55 * t79 - t61 * t81;
t102 = t55 * t81 + t61 * t79;
t98 = t75 * t76;
t97 = t76 * t80;
t93 = t82 * pkin(1) + t80 * qJ(2);
t90 = MDP(11) + MDP(15);
t87 = t80 * pkin(1) - t82 * qJ(2);
t86 = t104 * t82 + t93;
t65 = g(2) * t82 + g(3) * t80;
t84 = g(2) * t80 - g(3) * t82;
t83 = t104 * t80 + t87;
t60 = t78 * t89 - t94;
t59 = t78 * t95 + t88;
t58 = t76 * t77 * t91 - t78 * t74;
t54 = t60 * t77 + t74 * t97;
t53 = t54 * t81 + t59 * t79;
t52 = -t54 * t79 + t59 * t81;
t1 = [(-g(2) * t93 - g(3) * t87) * MDP(7) + (-g(2) * t62 - g(3) * t60) * MDP(8) + (-g(2) * t86 - g(3) * t83) * MDP(11) + (-g(2) * t55 - g(3) * t54) * MDP(12) + (-g(2) * (-t62 * t74 + t77 * t96) - g(3) * (-t60 * t74 + t77 * t97)) * MDP(13) + (-g(2) * (t62 * pkin(3) + t61 * qJ(4) + t86) - g(3) * (t60 * pkin(3) + t59 * qJ(4) + t83)) * MDP(15) + (-g(2) * t102 - g(3) * t53) * MDP(21) + (g(2) * t103 - g(3) * t52) * MDP(22) + (MDP(3) - MDP(6)) * t84 + (MDP(9) - MDP(14)) * (g(2) * t61 + g(3) * t59) + (-MDP(2) - MDP(4) * t78 + (MDP(5) - MDP(10)) * t76) * t65; (MDP(7) + t90) * t65; t90 * (g(1) * t78 - t76 * t84); (-g(1) * t98 - g(2) * t59 + g(3) * t61) * MDP(15); (-g(1) * (-t58 * t79 + t81 * t98) - g(2) * t52 - g(3) * t103) * MDP(21) + (-g(1) * (-t58 * t81 - t79 * t98) + g(2) * t53 - g(3) * t102) * MDP(22);];
taug = t1;
