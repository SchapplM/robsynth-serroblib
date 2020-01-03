% Calculate Gravitation load on the joints for
% S5RPRPR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR15_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRPR15_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR15_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR15_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S5RPRPR15_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:37:25
% EndTime: 2019-12-31 18:37:26
% DurationCPUTime: 0.32s
% Computational Cost: add. (107->54), mult. (189->78), div. (0->0), fcn. (168->8), ass. (0->29)
t60 = sin(qJ(3));
t62 = cos(qJ(3));
t65 = t60 * pkin(3) - t62 * qJ(4);
t61 = sin(qJ(1));
t63 = cos(qJ(1));
t83 = -g(1) * t61 + g(2) * t63;
t82 = MDP(13) - MDP(16);
t47 = -g(3) * t60 - t62 * t83;
t78 = g(3) * t62;
t57 = pkin(8) + qJ(5);
t51 = sin(t57);
t76 = t61 * t51;
t52 = cos(t57);
t75 = t61 * t52;
t58 = sin(pkin(8));
t74 = t61 * t58;
t59 = cos(pkin(8));
t73 = t61 * t59;
t72 = t63 * t51;
t71 = t63 * t52;
t70 = t63 * t58;
t69 = t63 * t59;
t68 = t63 * pkin(1) + t61 * qJ(2);
t54 = t63 * qJ(2);
t45 = t60 * t71 - t76;
t44 = t60 * t72 + t75;
t43 = t60 * t75 + t72;
t42 = -t60 * t76 + t71;
t1 = [(-g(1) * (-t61 * pkin(1) + t54) - g(2) * t68) * MDP(6) + (-g(1) * (t60 * t69 - t74) - g(2) * (t60 * t73 + t70)) * MDP(14) + (-g(1) * (-t60 * t70 - t73) - g(2) * (-t60 * t74 + t69)) * MDP(15) + (-g(1) * (t65 * t63 + t54) - g(2) * (t63 * pkin(6) + t68) + (-g(1) * (-pkin(1) - pkin(6)) - g(2) * t65) * t61) * MDP(17) + (-g(1) * t45 - g(2) * t43) * MDP(23) + (g(1) * t44 - g(2) * t42) * MDP(24) - (MDP(2) - MDP(4)) * t83 + (-t60 * MDP(12) - t82 * t62 + MDP(3) - MDP(5)) * (g(1) * t63 + g(2) * t61); -(-MDP(17) - MDP(6)) * t83; (g(3) * t65 + t83 * (pkin(3) * t62 + qJ(4) * t60)) * MDP(17) + t82 * (-t60 * t83 + t78) + (-MDP(14) * t59 + MDP(15) * t58 - MDP(23) * t52 + MDP(24) * t51 - MDP(12)) * t47; t47 * MDP(17); (-g(1) * t42 - g(2) * t44 + t51 * t78) * MDP(23) + (g(1) * t43 - g(2) * t45 + t52 * t78) * MDP(24);];
taug = t1;
