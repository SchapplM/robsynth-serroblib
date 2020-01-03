% Calculate potential energy for
% S5RRPRP1
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
% Datum: 2020-01-03 11:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPRP1_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP1_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_energypot_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP1_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP1_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:58:51
% EndTime: 2020-01-03 11:58:51
% DurationCPUTime: 0.25s
% Computational Cost: add. (130->45), mult. (92->34), div. (0->0), fcn. (61->8), ass. (0->20)
t57 = m(5) + m(6);
t56 = -mrSges(5,2) - mrSges(6,2);
t55 = m(6) * pkin(4) + mrSges(5,1) + mrSges(6,1);
t54 = -m(4) - t57;
t42 = sin(qJ(4));
t44 = cos(qJ(4));
t53 = t57 * pkin(3) + t56 * t42 + t55 * t44 + mrSges(4,1);
t52 = m(5) * pkin(7) - m(6) * (-qJ(5) - pkin(7)) - mrSges(4,2) + mrSges(5,3) + mrSges(6,3);
t51 = pkin(6) + pkin(5);
t43 = sin(qJ(1));
t38 = t43 * pkin(1);
t45 = cos(qJ(1));
t50 = t45 * pkin(1);
t40 = qJ(1) + qJ(2);
t37 = cos(t40);
t36 = sin(t40);
t35 = pkin(8) + t40;
t32 = cos(t35);
t31 = sin(t35);
t1 = (m(3) * t50 + t45 * mrSges(2,1) + t37 * mrSges(3,1) - t43 * mrSges(2,2) - t36 * mrSges(3,2) - mrSges(1,3) + t54 * (-pkin(2) * t37 - t50) + t53 * t32 + t52 * t31) * g(3) + (-m(3) * t38 - t43 * mrSges(2,1) - t36 * mrSges(3,1) - t45 * mrSges(2,2) - t37 * mrSges(3,2) - mrSges(1,2) + t54 * (pkin(2) * t36 + t38) + t52 * t32 - t53 * t31) * g(2) + (-m(2) * pkin(5) - m(3) * t51 - mrSges(1,1) - mrSges(2,3) - mrSges(3,3) - mrSges(4,3) + t56 * t44 - t55 * t42 + t54 * (qJ(3) + t51)) * g(1);
U = t1;
