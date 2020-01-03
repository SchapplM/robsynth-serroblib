% Calculate minimal parameter regressor of potential energy for
% S5PRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
% 
% Output:
% U_reg [1x17]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRPRR9_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR9_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR9_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR9_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:39:44
% EndTime: 2019-12-31 17:39:44
% DurationCPUTime: 0.03s
% Computational Cost: add. (54->18), mult. (48->26), div. (0->0), fcn. (52->8), ass. (0->13)
t93 = cos(qJ(4));
t92 = sin(qJ(4));
t88 = pkin(8) + qJ(2);
t86 = sin(t88);
t87 = cos(t88);
t80 = -t86 * t92 - t87 * t93;
t81 = t86 * t93 - t87 * t92;
t91 = g(1) * t80 - g(2) * t81;
t90 = cos(qJ(5));
t89 = sin(qJ(5));
t83 = -g(1) * t87 - g(2) * t86;
t82 = g(1) * t86 - g(2) * t87;
t1 = [-g(3) * qJ(1), 0, t83, t82, t83, -t82, -g(1) * (t87 * pkin(2) + t86 * qJ(3) + cos(pkin(8)) * pkin(1)) - g(2) * (t86 * pkin(2) - t87 * qJ(3) + sin(pkin(8)) * pkin(1)) - g(3) * (pkin(5) + qJ(1)), 0, t91, -g(1) * t81 - g(2) * t80, 0, 0, 0, 0, 0, g(3) * t89 + t91 * t90, g(3) * t90 - t91 * t89;];
U_reg = t1;
