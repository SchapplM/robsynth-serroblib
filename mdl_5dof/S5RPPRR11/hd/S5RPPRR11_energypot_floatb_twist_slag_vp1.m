% Calculate potential energy for
% S5RPPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPPRR11_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR11_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RPPRR11_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR11_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR11_energypot_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR11_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR11_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:05:29
% EndTime: 2019-12-31 18:05:29
% DurationCPUTime: 0.30s
% Computational Cost: add. (98->74), mult. (109->78), div. (0->0), fcn. (89->6), ass. (0->21)
t27 = rSges(6,3) + pkin(7);
t8 = sin(qJ(4));
t9 = sin(qJ(1));
t26 = t9 * t8;
t12 = cos(qJ(1));
t25 = t12 * t8;
t10 = cos(qJ(5));
t24 = t9 * t10;
t22 = t10 * t12;
t21 = pkin(5) + r_base(3);
t20 = t9 * pkin(1) + r_base(2);
t19 = pkin(2) + t21;
t18 = t9 * qJ(3) + t20;
t17 = t12 * pkin(1) + t9 * qJ(2) + r_base(1);
t16 = pkin(3) + t19;
t15 = t12 * pkin(6) + t18;
t14 = t12 * qJ(3) + t17;
t11 = cos(qJ(4));
t13 = rSges(5,1) * t8 + rSges(5,2) * t11;
t7 = sin(qJ(5));
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t12 - t9 * rSges(2,2) + r_base(1)) + g(2) * (t9 * rSges(2,1) + rSges(2,2) * t12 + r_base(2)) + g(3) * (rSges(2,3) + t21)) - m(3) * (g(1) * (-rSges(3,2) * t12 + t9 * rSges(3,3) + t17) + g(2) * (-t9 * rSges(3,2) + (-rSges(3,3) - qJ(2)) * t12 + t20) + g(3) * (rSges(3,1) + t21)) - m(4) * (g(1) * (t9 * rSges(4,2) + rSges(4,3) * t12 + t14) + g(2) * (t9 * rSges(4,3) + (-rSges(4,2) - qJ(2)) * t12 + t18) + g(3) * (rSges(4,1) + t19)) - m(5) * (g(1) * t14 + g(2) * t15 + g(3) * (rSges(5,1) * t11 - rSges(5,2) * t8 + t16) + (g(1) * (-rSges(5,3) - pkin(6)) + g(2) * t13) * t9 + (g(1) * t13 + g(2) * (rSges(5,3) - qJ(2))) * t12) - m(6) * (g(1) * (pkin(4) * t25 - t9 * pkin(6) + (t22 * t8 - t9 * t7) * rSges(6,1) + (-t25 * t7 - t24) * rSges(6,2) + t14) + g(2) * (pkin(4) * t26 - t12 * qJ(2) + (t12 * t7 + t24 * t8) * rSges(6,1) + (-t26 * t7 + t22) * rSges(6,2) + t15) + g(3) * (t27 * t8 + t16) + (g(3) * (rSges(6,1) * t10 - rSges(6,2) * t7 + pkin(4)) - (g(1) * t12 + g(2) * t9) * t27) * t11);
U = t1;
