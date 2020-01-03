% Calculate minimal parameter regressor of potential energy for
% S5RRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% 
% Output:
% U_reg [1x17]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPPR3_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR3_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR3_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR3_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:26:35
% EndTime: 2019-12-31 19:26:35
% DurationCPUTime: 0.04s
% Computational Cost: add. (58->22), mult. (40->29), div. (0->0), fcn. (34->8), ass. (0->15)
t103 = g(3) * (qJ(3) + pkin(6) + pkin(5));
t95 = qJ(1) + qJ(2);
t90 = sin(t95);
t97 = sin(qJ(1));
t102 = t97 * pkin(1) + pkin(2) * t90;
t91 = cos(t95);
t99 = cos(qJ(1));
t101 = t99 * pkin(1) + pkin(2) * t91;
t89 = pkin(8) + t95;
t85 = sin(t89);
t86 = cos(t89);
t100 = -g(1) * t85 + g(2) * t86;
t98 = cos(qJ(5));
t96 = sin(qJ(5));
t1 = [0, -g(1) * t99 - g(2) * t97, g(1) * t97 - g(2) * t99, 0, -g(1) * t91 - g(2) * t90, g(1) * t90 - g(2) * t91, -g(1) * t101 - g(2) * t102 - t103, g(1) * t86 + g(2) * t85, t100, -g(1) * (t86 * pkin(3) + t85 * qJ(4) + t101) - g(2) * (t85 * pkin(3) - t86 * qJ(4) + t102) - t103, 0, 0, 0, 0, 0, -g(3) * t98 + t100 * t96, g(3) * t96 + t100 * t98;];
U_reg = t1;
