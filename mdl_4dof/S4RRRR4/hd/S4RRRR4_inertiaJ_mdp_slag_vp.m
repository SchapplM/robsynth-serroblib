% Calculate joint inertia matrix for
% S4RRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRRR4_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRRR4_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR4_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR4_inertiaJ_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S4RRRR4_inertiaJ_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:26:15
% EndTime: 2019-12-31 17:26:16
% DurationCPUTime: 0.20s
% Computational Cost: add. (197->68), mult. (391->100), div. (0->0), fcn. (375->6), ass. (0->38)
t73 = sin(qJ(4));
t76 = cos(qJ(4));
t91 = t73 * MDP(20) + t76 * MDP(21);
t82 = MDP(23) * t76 - MDP(24) * t73;
t78 = cos(qJ(2));
t67 = -pkin(2) * t78 - pkin(1);
t97 = 0.2e1 * t67;
t96 = 0.2e1 * t78;
t95 = pkin(5) + pkin(6);
t75 = sin(qJ(2));
t61 = t95 * t75;
t62 = t95 * t78;
t74 = sin(qJ(3));
t77 = cos(qJ(3));
t51 = t61 * t77 + t62 * t74;
t48 = t51 * t73;
t93 = t51 * t76;
t92 = t73 * t76;
t59 = t74 * t75 - t77 * t78;
t90 = MDP(22) * t59;
t87 = MDP(19) * t92;
t71 = t73 ^ 2;
t86 = t71 * MDP(18) + MDP(15) + 0.2e1 * t87;
t60 = t74 * t78 + t77 * t75;
t85 = -pkin(3) * t60 - pkin(7) * t59;
t65 = pkin(2) * t74 + pkin(7);
t66 = -pkin(2) * t77 - pkin(3);
t84 = -t59 * t65 + t60 * t66;
t83 = MDP(20) * t76 - MDP(21) * t73;
t81 = -MDP(23) * t73 - MDP(24) * t76;
t80 = (MDP(16) * t77 - MDP(17) * t74) * pkin(2);
t52 = -t61 * t74 + t62 * t77;
t72 = t76 ^ 2;
t79 = -t51 * MDP(16) - t52 * MDP(17) + ((-t71 + t72) * MDP(19) + MDP(18) * t92 + MDP(13)) * t60 + (-MDP(14) + t91) * t59;
t47 = t59 * pkin(3) - t60 * pkin(7) + t67;
t45 = t47 * t73 + t52 * t76;
t44 = t47 * t76 - t52 * t73;
t1 = [t60 * MDP(17) * t97 + pkin(1) * MDP(9) * t96 + MDP(1) + (MDP(16) * t97 + t90 + 0.2e1 * (-MDP(12) + t83) * t60) * t59 + 0.2e1 * (t44 * t59 + t60 * t48) * MDP(23) + 0.2e1 * (-t45 * t59 + t60 * t93) * MDP(24) + (-0.2e1 * pkin(1) * MDP(10) + MDP(4) * t75 + MDP(5) * t96) * t75 + (t72 * MDP(18) + MDP(11) - 0.2e1 * t87) * t60 ^ 2; t75 * MDP(6) + t78 * MDP(7) + (t84 * t73 - t93) * MDP(23) + (t84 * t76 + t48) * MDP(24) + (-MDP(10) * t78 - MDP(9) * t75) * pkin(5) + t79; -0.2e1 * t66 * t82 + MDP(8) + 0.2e1 * t80 + t86; (t85 * t73 - t93) * MDP(23) + (t85 * t76 + t48) * MDP(24) + t79; t80 + t86 + t82 * (pkin(3) - t66); 0.2e1 * pkin(3) * t82 + t86; t44 * MDP(23) - t45 * MDP(24) + t83 * t60 + t90; t81 * t65 + t91; t81 * pkin(7) + t91; MDP(22);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
