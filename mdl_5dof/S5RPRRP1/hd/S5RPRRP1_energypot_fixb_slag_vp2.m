% Calculate potential energy for
% S5RPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:00
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRRP1_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP1_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP1_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP1_energypot_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP1_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP1_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:59:10
% EndTime: 2019-12-05 17:59:10
% DurationCPUTime: 0.25s
% Computational Cost: add. (96->43), mult. (110->29), div. (0->0), fcn. (79->6), ass. (0->16)
t59 = m(5) + m(6);
t58 = mrSges(5,2) + mrSges(6,2);
t57 = m(6) * pkin(4) + mrSges(5,1) + mrSges(6,1);
t55 = m(3) + m(4) + t59;
t40 = qJ(3) + qJ(4);
t33 = sin(t40);
t34 = cos(t40);
t41 = sin(qJ(3));
t43 = cos(qJ(3));
t54 = t43 * mrSges(4,2) + t57 * t33 + t58 * t34 - mrSges(2,2) + mrSges(3,3) + (pkin(3) * t59 + mrSges(4,1)) * t41;
t45 = -pkin(7) - pkin(6);
t53 = -m(4) * pkin(6) + m(5) * t45 + m(6) * (-qJ(5) + t45) - mrSges(2,1) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3) - mrSges(6,3);
t52 = pkin(2) + pkin(5);
t44 = cos(qJ(1));
t42 = sin(qJ(1));
t1 = (-m(4) * t52 - t43 * mrSges(4,1) + t41 * mrSges(4,2) - mrSges(3,1) - mrSges(1,3) - mrSges(2,3) - t59 * (t43 * pkin(3) + t52) - t57 * t34 + t58 * t33 + (-m(2) - m(3)) * pkin(5)) * g(3) + (-mrSges(1,2) + (t55 * qJ(2) + t54) * t44 + (-t55 * pkin(1) + t53) * t42) * g(2) + (-mrSges(1,1) - t55 * (t44 * pkin(1) + t42 * qJ(2)) + t53 * t44 - t54 * t42) * g(1);
U = t1;
