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
% U_reg [1x26]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2019-12-31 21:14:16
% EndTime: 2019-12-31 21:14:16
% DurationCPUTime: 0.05s
% Computational Cost: add. (54->27), mult. (59->38), div. (0->0), fcn. (60->10), ass. (0->21)
t172 = qJ(2) + qJ(3);
t168 = pkin(9) + t172;
t184 = g(3) * sin(t168);
t173 = sin(qJ(5));
t175 = sin(qJ(1));
t183 = t175 * t173;
t176 = cos(qJ(5));
t182 = t175 * t176;
t178 = cos(qJ(1));
t181 = t178 * t173;
t180 = t178 * t176;
t179 = g(1) * t178 + g(2) * t175;
t177 = cos(qJ(2));
t174 = sin(qJ(2));
t171 = -qJ(4) - pkin(7) - pkin(6);
t170 = cos(t172);
t169 = sin(t172);
t167 = cos(t168);
t165 = g(1) * t175 - g(2) * t178;
t164 = t177 * pkin(2) + pkin(3) * t170 + pkin(1);
t1 = [0, -t179, t165, 0, 0, 0, 0, 0, -g(3) * t174 - t179 * t177, -g(3) * t177 + t179 * t174, 0, 0, 0, 0, 0, -g(3) * t169 - t179 * t170, -g(3) * t170 + t179 * t169, -t165, -g(1) * (t178 * t164 - t175 * t171) - g(2) * (t175 * t164 + t178 * t171) - g(3) * (t174 * pkin(2) + pkin(3) * t169 + pkin(5)), 0, 0, 0, 0, 0, -g(1) * (t167 * t180 + t183) - g(2) * (t167 * t182 - t181) - t176 * t184, -g(1) * (-t167 * t181 + t182) - g(2) * (-t167 * t183 - t180) + t173 * t184;];
U_reg = t1;
