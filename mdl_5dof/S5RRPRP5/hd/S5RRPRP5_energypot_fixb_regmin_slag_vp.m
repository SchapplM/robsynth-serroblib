% Calculate minimal parameter regressor of potential energy for
% S5RRPRP5
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
% Datum: 2019-12-31 19:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPRP5_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP5_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP5_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP5_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:54:56
% EndTime: 2019-12-31 19:54:56
% DurationCPUTime: 0.05s
% Computational Cost: add. (84->30), mult. (71->34), div. (0->0), fcn. (65->8), ass. (0->18)
t154 = -pkin(6) - qJ(3);
t155 = sin(qJ(2));
t161 = t155 * pkin(2) + pkin(5);
t157 = cos(qJ(2));
t147 = t157 * pkin(2) + pkin(1);
t153 = qJ(2) + pkin(8);
t156 = sin(qJ(1));
t158 = cos(qJ(1));
t160 = g(1) * t158 + g(2) * t156;
t148 = qJ(4) + t153;
t145 = sin(t148);
t146 = cos(t148);
t159 = pkin(4) * t146 + qJ(5) * t145 + pkin(3) * cos(t153) + t147;
t152 = -pkin(7) + t154;
t144 = g(1) * t156 - g(2) * t158;
t142 = -g(3) * t145 - t160 * t146;
t141 = -g(3) * t146 + t160 * t145;
t1 = [0, -t160, t144, 0, 0, 0, 0, 0, -g(3) * t155 - t160 * t157, -g(3) * t157 + t160 * t155, -t144, -g(1) * (t158 * t147 - t156 * t154) - g(2) * (t156 * t147 + t158 * t154) - g(3) * t161, 0, 0, 0, 0, 0, t142, t141, t142, -t144, -t141, -g(3) * (t145 * pkin(4) - t146 * qJ(5) + pkin(3) * sin(t153) + t161) + (-g(1) * t159 - g(2) * t152) * t158 + (g(1) * t152 - g(2) * t159) * t156;];
U_reg = t1;
