% Calculate minimal parameter regressor of potential energy for
% S5RRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% U_reg [1x17]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPRR3_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR3_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR3_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR3_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:00:34
% EndTime: 2020-01-03 12:00:34
% DurationCPUTime: 0.03s
% Computational Cost: add. (46->17), mult. (29->23), div. (0->0), fcn. (26->8), ass. (0->12)
t100 = qJ(1) + qJ(2);
t97 = pkin(9) + qJ(4) + t100;
t95 = sin(t97);
t96 = cos(t97);
t105 = g(2) * t95 - g(3) * t96;
t104 = cos(qJ(1));
t103 = cos(qJ(5));
t102 = sin(qJ(1));
t101 = sin(qJ(5));
t99 = cos(t100);
t98 = sin(t100);
t1 = [0, -g(2) * t102 + g(3) * t104, -g(2) * t104 - g(3) * t102, 0, -g(2) * t98 + g(3) * t99, -g(2) * t99 - g(3) * t98, -g(1) * (qJ(3) + pkin(6) + pkin(5)) - g(2) * (t102 * pkin(1) + pkin(2) * t98) - g(3) * (-t104 * pkin(1) - pkin(2) * t99), 0, -t105, -g(2) * t96 - g(3) * t95, 0, 0, 0, 0, 0, -g(1) * t101 - t105 * t103, -g(1) * t103 + t105 * t101;];
U_reg = t1;
