% Calculate minimal parameter regressor of potential energy for
% S4PRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
% 
% Output:
% U_reg [1x15]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:24
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4PRPR6_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR6_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR6_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR6_energypot_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:24:41
% EndTime: 2019-12-31 16:24:41
% DurationCPUTime: 0.06s
% Computational Cost: add. (43->30), mult. (73->48), div. (0->0), fcn. (77->8), ass. (0->16)
t91 = sin(qJ(2));
t97 = g(3) * t91;
t88 = sin(pkin(6));
t92 = cos(qJ(2));
t96 = t88 * t92;
t90 = cos(pkin(6));
t95 = t90 * t92;
t94 = g(1) * t90 + g(2) * t88;
t93 = pkin(2) * t92 + qJ(3) * t91 + pkin(1);
t89 = cos(pkin(7));
t87 = sin(pkin(7));
t86 = pkin(7) + qJ(4);
t85 = cos(t86);
t84 = sin(t86);
t83 = -g(3) * t92 + t94 * t91;
t1 = [-g(3) * qJ(1), 0, -t94 * t92 - t97, t83, -g(1) * (t88 * t87 + t89 * t95) - g(2) * (-t90 * t87 + t89 * t96) - t89 * t97, -g(1) * (-t87 * t95 + t88 * t89) - g(2) * (-t87 * t96 - t90 * t89) + t87 * t97, -t83, -g(3) * (t91 * pkin(2) - t92 * qJ(3) + qJ(1)) + (g(2) * pkin(4) - g(1) * t93) * t90 + (-g(1) * pkin(4) - g(2) * t93) * t88, 0, 0, 0, 0, 0, -g(1) * (t88 * t84 + t85 * t95) - g(2) * (-t90 * t84 + t85 * t96) - t85 * t97, -g(1) * (-t84 * t95 + t88 * t85) - g(2) * (-t84 * t96 - t90 * t85) + t84 * t97;];
U_reg = t1;
