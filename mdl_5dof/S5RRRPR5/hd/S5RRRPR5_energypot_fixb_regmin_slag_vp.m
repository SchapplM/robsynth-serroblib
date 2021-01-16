% Calculate minimal parameter regressor of potential energy for
% S5RRRPR5
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
% Datum: 2021-01-15 23:12
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRPR5_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR5_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR5_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR5_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 23:10:42
% EndTime: 2021-01-15 23:10:42
% DurationCPUTime: 0.06s
% Computational Cost: add. (70->29), mult. (69->41), div. (0->0), fcn. (70->10), ass. (0->22)
t83 = qJ(2) + qJ(3);
t79 = pkin(9) + t83;
t77 = sin(t79);
t95 = g(3) * t77;
t84 = sin(qJ(5));
t86 = sin(qJ(1));
t94 = t86 * t84;
t87 = cos(qJ(5));
t93 = t86 * t87;
t89 = cos(qJ(1));
t92 = t89 * t84;
t91 = t89 * t87;
t90 = g(1) * t89 + g(2) * t86;
t88 = cos(qJ(2));
t85 = sin(qJ(2));
t82 = -qJ(4) - pkin(7) - pkin(6);
t81 = cos(t83);
t80 = sin(t83);
t78 = cos(t79);
t76 = g(1) * t86 - g(2) * t89;
t75 = t88 * pkin(2) + pkin(3) * t81 + pkin(1);
t1 = [0, -t90, t76, 0, 0, 0, 0, 0, -g(3) * t85 - t90 * t88, -g(3) * t88 + t90 * t85, 0, 0, 0, 0, 0, -g(3) * t80 - t90 * t81, -g(3) * t81 + t90 * t80, -t90 * t78 - t95, -g(3) * t78 + t90 * t77, -t76, -g(1) * (t89 * t75 - t86 * t82) - g(2) * (t86 * t75 + t89 * t82) - g(3) * (t85 * pkin(2) + pkin(3) * t80 + pkin(5)), 0, 0, 0, 0, 0, -g(1) * (t78 * t91 + t94) - g(2) * (t78 * t93 - t92) - t87 * t95, -g(1) * (-t78 * t92 + t93) - g(2) * (-t78 * t94 - t91) + t84 * t95;];
U_reg = t1;
