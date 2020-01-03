% Calculate potential energy for
% S5RPPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
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
% Datum: 2019-12-31 17:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPPRP3_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP3_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP3_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP3_energypot_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP3_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRP3_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:50:41
% EndTime: 2019-12-31 17:50:41
% DurationCPUTime: 0.16s
% Computational Cost: add. (111->55), mult. (89->59), div. (0->0), fcn. (65->6), ass. (0->20)
t49 = rSges(6,1) + pkin(4);
t48 = rSges(5,3) + pkin(6);
t47 = rSges(6,3) + qJ(5) + pkin(6);
t46 = pkin(5) + qJ(2);
t34 = qJ(1) + pkin(7);
t30 = sin(t34);
t37 = sin(qJ(1));
t32 = t37 * pkin(1);
t45 = t30 * pkin(2) + t32;
t44 = pkin(3) + t46;
t31 = cos(t34);
t39 = cos(qJ(1));
t33 = t39 * pkin(1);
t43 = t31 * pkin(2) + t30 * qJ(3) + t33;
t36 = sin(qJ(4));
t38 = cos(qJ(4));
t42 = rSges(5,1) * t36 + rSges(5,2) * t38;
t41 = rSges(6,2) * t38 + t49 * t36;
t40 = g(1) * t43 + g(2) * t45;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t39 * rSges(2,1) - t37 * rSges(2,2)) + g(2) * (t37 * rSges(2,1) + t39 * rSges(2,2)) + g(3) * (pkin(5) + rSges(2,3))) - m(3) * (g(1) * (t31 * rSges(3,1) - t30 * rSges(3,2) + t33) + g(2) * (t30 * rSges(3,1) + t31 * rSges(3,2) + t32) + g(3) * (rSges(3,3) + t46)) - m(4) * (g(1) * (-t31 * rSges(4,2) + t30 * rSges(4,3) + t43) + g(2) * (-t30 * rSges(4,2) + (-rSges(4,3) - qJ(3)) * t31 + t45) + g(3) * (rSges(4,1) + t46)) - m(5) * (g(3) * (t38 * rSges(5,1) - t36 * rSges(5,2) + t44) + (g(1) * t42 + g(2) * t48) * t30 + (g(1) * t48 + g(2) * (-qJ(3) - t42)) * t31 + t40) - m(6) * (g(3) * (-t36 * rSges(6,2) + t49 * t38 + t44) + (g(1) * t41 + g(2) * t47) * t30 + (g(1) * t47 + g(2) * (-qJ(3) - t41)) * t31 + t40);
U = t1;
