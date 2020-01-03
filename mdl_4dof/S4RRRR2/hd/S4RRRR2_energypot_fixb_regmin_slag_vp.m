% Calculate minimal parameter regressor of potential energy for
% S4RRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% 
% Output:
% U_reg [1x20]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RRRR2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR2_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR2_energypot_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:23:18
% EndTime: 2019-12-31 17:23:18
% DurationCPUTime: 0.03s
% Computational Cost: add. (30->11), mult. (28->16), div. (0->0), fcn. (28->8), ass. (0->12)
t94 = qJ(1) + qJ(2);
t90 = sin(t94);
t92 = cos(t94);
t99 = g(1) * t92 + g(2) * t90;
t98 = cos(qJ(1));
t97 = cos(qJ(3));
t96 = sin(qJ(1));
t95 = sin(qJ(3));
t93 = qJ(3) + qJ(4);
t91 = cos(t93);
t89 = sin(t93);
t1 = [0, -g(1) * t98 - g(2) * t96, g(1) * t96 - g(2) * t98, 0, -t99, g(1) * t90 - g(2) * t92, 0, 0, 0, 0, 0, -g(3) * t95 - t99 * t97, -g(3) * t97 + t99 * t95, 0, 0, 0, 0, 0, -g(3) * t89 - t99 * t91, -g(3) * t91 + t99 * t89;];
U_reg = t1;
