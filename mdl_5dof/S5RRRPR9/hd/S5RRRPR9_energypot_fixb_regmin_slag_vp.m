% Calculate minimal parameter regressor of potential energy for
% S5RRRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% U_reg [1x28]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 23:25
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRPR9_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR9_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR9_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR9_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 23:23:54
% EndTime: 2021-01-15 23:23:54
% DurationCPUTime: 0.08s
% Computational Cost: add. (80->43), mult. (97->64), div. (0->0), fcn. (106->10), ass. (0->27)
t86 = sin(qJ(2));
t99 = g(3) * t86;
t87 = sin(qJ(1));
t89 = cos(qJ(2));
t98 = t87 * t89;
t83 = qJ(3) + pkin(9);
t82 = qJ(5) + t83;
t76 = sin(t82);
t90 = cos(qJ(1));
t97 = t90 * t76;
t77 = cos(t82);
t96 = t90 * t77;
t80 = sin(t83);
t95 = t90 * t80;
t81 = cos(t83);
t94 = t90 * t81;
t85 = sin(qJ(3));
t93 = t90 * t85;
t88 = cos(qJ(3));
t92 = t90 * t88;
t91 = g(1) * t90 + g(2) * t87;
t84 = qJ(4) + pkin(7);
t79 = t88 * pkin(3) + pkin(2);
t78 = t85 * pkin(3) + pkin(6);
t75 = t79 * t89 + t84 * t86 + pkin(1);
t74 = -g(3) * t89 + t91 * t86;
t1 = [0, -t91, g(1) * t87 - g(2) * t90, 0, 0, 0, 0, 0, -t91 * t89 - t99, t74, 0, 0, 0, 0, 0, -g(1) * (t87 * t85 + t89 * t92) - g(2) * (t88 * t98 - t93) - t88 * t99, -g(1) * (t87 * t88 - t89 * t93) - g(2) * (-t85 * t98 - t92) + t85 * t99, -g(1) * (t87 * t80 + t89 * t94) - g(2) * (t81 * t98 - t95) - t81 * t99, -g(1) * (t87 * t81 - t89 * t95) - g(2) * (-t80 * t98 - t94) + t80 * t99, -t74, -g(1) * (t75 * t90 + t78 * t87) - g(2) * (t75 * t87 - t78 * t90) - g(3) * (t86 * t79 - t89 * t84 + pkin(5)), 0, 0, 0, 0, 0, -g(1) * (t87 * t76 + t89 * t96) - g(2) * (t77 * t98 - t97) - t77 * t99, -g(1) * (t87 * t77 - t89 * t97) - g(2) * (-t76 * t98 - t96) + t76 * t99;];
U_reg = t1;
