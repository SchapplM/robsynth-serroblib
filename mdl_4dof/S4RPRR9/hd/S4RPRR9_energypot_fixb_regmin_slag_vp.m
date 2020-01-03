% Calculate minimal parameter regressor of potential energy for
% S4RPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
% 
% Output:
% U_reg [1x20]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RPRR9_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR9_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR9_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR9_energypot_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:56:28
% EndTime: 2019-12-31 16:56:28
% DurationCPUTime: 0.04s
% Computational Cost: add. (20->18), mult. (45->29), div. (0->0), fcn. (46->6), ass. (0->14)
t85 = cos(qJ(3));
t91 = g(3) * t85;
t81 = sin(qJ(4));
t83 = sin(qJ(1));
t90 = t83 * t81;
t84 = cos(qJ(4));
t89 = t83 * t84;
t86 = cos(qJ(1));
t88 = t86 * t81;
t87 = t86 * t84;
t79 = g(1) * t83 - g(2) * t86;
t82 = sin(qJ(3));
t80 = g(1) * t86 + g(2) * t83;
t1 = [0, -t80, t79, t80, -t79, -g(1) * (t86 * pkin(1) + t83 * qJ(2)) - g(2) * (t83 * pkin(1) - t86 * qJ(2)) - g(3) * pkin(4), 0, 0, 0, 0, 0, -t79 * t82 - t91, g(3) * t82 - t79 * t85, 0, 0, 0, 0, 0, -g(1) * (t82 * t89 + t88) - g(2) * (-t82 * t87 + t90) - t84 * t91, -g(1) * (-t82 * t90 + t87) - g(2) * (t82 * t88 + t89) + t81 * t91;];
U_reg = t1;
