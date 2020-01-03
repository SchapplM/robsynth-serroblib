% Calculate minimal parameter regressor of potential energy for
% S4RRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% U_reg [1x17]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RRPR4_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR4_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR4_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR4_energypot_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:02:33
% EndTime: 2019-12-31 17:02:33
% DurationCPUTime: 0.03s
% Computational Cost: add. (44->19), mult. (39->25), div. (0->0), fcn. (36->8), ass. (0->13)
t90 = qJ(1) + qJ(2);
t87 = sin(t90);
t88 = cos(t90);
t95 = g(1) * t88 + g(2) * t87;
t94 = cos(qJ(1));
t93 = sin(qJ(1));
t92 = cos(pkin(7));
t91 = sin(pkin(7));
t89 = pkin(7) + qJ(4);
t86 = cos(t89);
t85 = sin(t89);
t84 = g(1) * t87 - g(2) * t88;
t1 = [0, -g(1) * t94 - g(2) * t93, g(1) * t93 - g(2) * t94, 0, -t95, t84, -g(3) * t91 - t95 * t92, -g(3) * t92 + t95 * t91, -t84, -g(1) * (t94 * pkin(1) + t88 * pkin(2) + t87 * qJ(3)) - g(2) * (t93 * pkin(1) + t87 * pkin(2) - t88 * qJ(3)) - g(3) * (pkin(5) + pkin(4)), 0, 0, 0, 0, 0, -g(3) * t85 - t95 * t86, -g(3) * t86 + t95 * t85;];
U_reg = t1;
