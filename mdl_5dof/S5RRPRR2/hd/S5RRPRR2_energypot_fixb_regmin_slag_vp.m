% Calculate minimal parameter regressor of potential energy for
% S5RRPRR2
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
% U_reg [1x26]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPRR2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR2_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR2_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:28:26
% EndTime: 2019-12-05 18:28:26
% DurationCPUTime: 0.04s
% Computational Cost: add. (54->19), mult. (46->25), div. (0->0), fcn. (43->8), ass. (0->15)
t167 = qJ(2) + pkin(9) + qJ(4);
t170 = sin(qJ(1));
t172 = cos(qJ(1));
t173 = g(1) * t172 + g(2) * t170;
t171 = cos(qJ(2));
t169 = sin(qJ(2));
t168 = -pkin(6) - qJ(3);
t166 = t171 * pkin(2) + pkin(1);
t165 = qJ(5) + t167;
t164 = cos(t167);
t163 = sin(t167);
t162 = cos(t165);
t161 = sin(t165);
t160 = g(1) * t170 - g(2) * t172;
t1 = [0, -t173, t160, 0, 0, 0, 0, 0, -g(3) * t169 - t173 * t171, -g(3) * t171 + t173 * t169, -t160, -g(1) * (t172 * t166 - t170 * t168) - g(2) * (t170 * t166 + t172 * t168) - g(3) * (t169 * pkin(2) + pkin(5)), 0, 0, 0, 0, 0, -g(3) * t163 - t173 * t164, -g(3) * t164 + t173 * t163, 0, 0, 0, 0, 0, -g(3) * t161 - t173 * t162, -g(3) * t162 + t173 * t161;];
U_reg = t1;
