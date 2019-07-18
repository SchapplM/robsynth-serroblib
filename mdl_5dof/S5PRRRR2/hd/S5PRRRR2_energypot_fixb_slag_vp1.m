% Calculate potential energy for
% S5PRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,d5]';
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
% Datum: 2019-07-18 13:30
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRRRR2_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(6,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR2_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR2_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5PRRRR2_energypot_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR2_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRR2_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:30:09
% EndTime: 2019-07-18 13:30:09
% DurationCPUTime: 0.08s
% Computational Cost: add. (86->47), mult. (60->50), div. (0->0), fcn. (36->8), ass. (0->19)
t45 = rSges(6,3) + pkin(6);
t38 = cos(qJ(2));
t44 = t38 * pkin(2) + pkin(1);
t43 = pkin(4) + qJ(1);
t34 = qJ(2) + qJ(3);
t29 = sin(t34);
t36 = sin(qJ(2));
t32 = t36 * pkin(2);
t42 = pkin(3) * t29 + t32;
t30 = cos(t34);
t41 = pkin(3) * t30 + t44;
t40 = pkin(5) + t43;
t35 = sin(qJ(5));
t37 = cos(qJ(5));
t39 = rSges(6,1) * t37 - rSges(6,2) * t35;
t31 = qJ(4) + t34;
t28 = cos(t31);
t27 = sin(t31);
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * rSges(2,1) + g(2) * rSges(2,2) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (t38 * rSges(3,1) - t36 * rSges(3,2) + pkin(1)) + g(2) * (t36 * rSges(3,1) + t38 * rSges(3,2)) + g(3) * (qJ(1) + rSges(3,3))) - m(4) * (g(1) * (t30 * rSges(4,1) - t29 * rSges(4,2) + t44) + g(2) * (t29 * rSges(4,1) + t30 * rSges(4,2) + t32) + g(3) * (rSges(4,3) + t43)) - m(5) * (g(1) * (t28 * rSges(5,1) - t27 * rSges(5,2) + t41) + g(2) * (t27 * rSges(5,1) + t28 * rSges(5,2) + t42) + g(3) * (rSges(5,3) + t40)) - m(6) * (g(1) * t41 + g(2) * t42 + g(3) * (t35 * rSges(6,1) + t37 * rSges(6,2) + t40) + (g(1) * t39 - g(2) * t45) * t28 + (g(1) * t45 + g(2) * t39) * t27);
U  = t1;
