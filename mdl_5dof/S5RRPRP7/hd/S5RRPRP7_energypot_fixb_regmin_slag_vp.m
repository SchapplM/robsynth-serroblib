% Calculate minimal parameter regressor of potential energy for
% S5RRPRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% U_reg [1x23]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPRP7_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP7_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP7_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP7_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:01:31
% EndTime: 2019-12-31 20:01:31
% DurationCPUTime: 0.06s
% Computational Cost: add. (86->37), mult. (111->50), div. (0->0), fcn. (117->8), ass. (0->29)
t174 = qJ(2) + pkin(8);
t172 = cos(t174);
t192 = pkin(3) * t172;
t171 = sin(t174);
t191 = g(3) * t171;
t177 = sin(qJ(2));
t190 = t177 * pkin(2) + pkin(5);
t176 = sin(qJ(4));
t178 = sin(qJ(1));
t189 = t178 * t176;
t179 = cos(qJ(4));
t188 = t178 * t179;
t181 = cos(qJ(1));
t187 = t181 * t176;
t186 = t181 * t179;
t180 = cos(qJ(2));
t170 = t180 * pkin(2) + pkin(1);
t175 = -pkin(6) - qJ(3);
t185 = t178 * t170 + t181 * t175;
t184 = t181 * t170 - t178 * t175;
t183 = g(1) * t181 + g(2) * t178;
t161 = t172 * t189 + t186;
t163 = t172 * t187 - t188;
t182 = g(1) * t163 + g(2) * t161 + t176 * t191;
t165 = g(1) * t178 - g(2) * t181;
t164 = t172 * t186 + t189;
t162 = t172 * t188 - t187;
t160 = -g(1) * t164 - g(2) * t162 - t179 * t191;
t1 = [0, -t183, t165, 0, 0, 0, 0, 0, -g(3) * t177 - t183 * t180, -g(3) * t180 + t183 * t177, -t165, -g(1) * t184 - g(2) * t185 - g(3) * t190, 0, 0, 0, 0, 0, t160, t182, t160, g(3) * t172 - t183 * t171, -t182, -g(1) * (t164 * pkin(4) + t163 * qJ(5) + t181 * t192 + t184) - g(2) * (t162 * pkin(4) + t161 * qJ(5) + t178 * t192 + t185) - g(3) * (-t172 * pkin(7) + t190) + (-g(3) * (pkin(4) * t179 + qJ(5) * t176 + pkin(3)) - t183 * pkin(7)) * t171;];
U_reg = t1;
