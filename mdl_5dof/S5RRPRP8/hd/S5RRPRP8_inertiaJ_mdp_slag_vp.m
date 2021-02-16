% Calculate joint inertia matrix for
% S5RRPRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRP8_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 20:52
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRP8_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP8_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP8_inertiaJ_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RRPRP8_inertiaJ_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 20:52:21
% EndTime: 2021-01-15 20:52:22
% DurationCPUTime: 0.22s
% Computational Cost: add. (256->91), mult. (412->123), div. (0->0), fcn. (348->4), ass. (0->33)
t67 = sin(qJ(4));
t68 = sin(qJ(2));
t69 = cos(qJ(4));
t70 = cos(qJ(2));
t49 = -t68 * t67 - t70 * t69;
t85 = 0.2e1 * t49;
t82 = -pkin(2) - pkin(3);
t55 = t69 * qJ(3) + t67 * t82;
t75 = MDP(21) + MDP(23);
t84 = t75 * t55;
t83 = 2 * MDP(22);
t81 = pkin(6) - pkin(7);
t65 = t68 ^ 2;
t80 = t70 ^ 2 + t65;
t79 = MDP(25) * pkin(4);
t54 = t67 * qJ(3) - t69 * t82;
t52 = pkin(4) + t54;
t78 = t52 * MDP(25);
t77 = t54 * MDP(20);
t76 = MDP(20) + MDP(22);
t56 = -t70 * pkin(2) - t68 * qJ(3) - pkin(1);
t74 = t81 * t68;
t57 = t81 * t70;
t46 = t67 * t57 - t69 * t74;
t48 = t70 * pkin(3) - t56;
t73 = -pkin(2) * MDP(14) - MDP(11);
t47 = t69 * t57 + t67 * t74;
t50 = -t70 * t67 + t68 * t69;
t42 = t50 * qJ(5) + t46;
t44 = t49 * qJ(5) + t47;
t72 = t50 * MDP(17) + t49 * MDP(18) - t46 * MDP(20) - t47 * MDP(21) - t42 * MDP(22) - t44 * MDP(23);
t45 = -t49 * pkin(4) + t48;
t1 = [MDP(1) + t65 * MDP(4) + (t80 * pkin(6) ^ 2 + t56 ^ 2) * MDP(14) + (t42 ^ 2 + t44 ^ 2 + t45 ^ 2) * MDP(25) + 0.2e1 * t80 * MDP(12) * pkin(6) + (MDP(15) * t50 + MDP(16) * t85 + 0.2e1 * t48 * MDP(21) + 0.2e1 * t45 * MDP(23) + 0.2e1 * MDP(24) * t42) * t50 + 0.2e1 * (-t56 * MDP(11) + pkin(1) * MDP(9)) * t70 + (-t48 * MDP(20) - t45 * MDP(22) + MDP(24) * t44) * t85 + 0.2e1 * (-pkin(1) * MDP(10) - t56 * MDP(13) + t70 * MDP(5)) * t68; t68 * MDP(6) + t70 * MDP(7) + (-t68 * pkin(2) + t70 * qJ(3)) * MDP(12) + (t55 * t49 + t52 * t50) * MDP(24) + (t42 * t52 + t44 * t55) * MDP(25) + ((MDP(14) * qJ(3) - MDP(10) + MDP(13)) * t70 + (-MDP(9) + t73) * t68) * pkin(6) - t72; MDP(8) + 0.2e1 * pkin(2) * MDP(11) + 0.2e1 * qJ(3) * MDP(13) + (pkin(2) ^ 2 + qJ(3) ^ 2) * MDP(14) + MDP(19) + (t52 ^ 2 + t55 ^ 2) * MDP(25) + 0.2e1 * t77 + t52 * t83 + 0.2e1 * t84; (t67 * t49 - t69 * t50) * MDP(24) + (-t42 * t69 + t44 * t67) * MDP(25) + (pkin(6) * MDP(14) + MDP(12)) * t68; (-t76 - t78) * t69 + (t55 * MDP(25) + t75) * t67 + t73; MDP(14) + (t67 ^ 2 + t69 ^ 2) * MDP(25); (-t50 * MDP(24) - t42 * MDP(25)) * pkin(4) + t72; -MDP(19) - t77 + (-0.2e1 * pkin(4) - t54) * MDP(22) - pkin(4) * t78 - t84; -t75 * t67 + (t76 + t79) * t69; MDP(19) + (t83 + t79) * pkin(4); -t49 * MDP(22) + t50 * MDP(23) + t45 * MDP(25); 0; 0; 0; MDP(25);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
