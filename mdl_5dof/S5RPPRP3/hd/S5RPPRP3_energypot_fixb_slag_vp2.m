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
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPPRP3_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP3_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP3_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP3_energypot_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP3_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP3_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:50:41
% EndTime: 2019-12-31 17:50:41
% DurationCPUTime: 0.23s
% Computational Cost: add. (111->44), mult. (96->31), div. (0->0), fcn. (65->6), ass. (0->17)
t56 = -m(6) * pkin(4) - mrSges(5,1) - mrSges(6,1);
t55 = mrSges(5,2) + mrSges(6,2);
t54 = -m(5) - m(6);
t53 = m(4) - t54;
t51 = -m(5) * pkin(6) + m(6) * (-qJ(5) - pkin(6)) - mrSges(3,1) + mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t41 = sin(qJ(4));
t43 = cos(qJ(4));
t50 = t56 * t41 - t55 * t43 + mrSges(3,2) - mrSges(4,3);
t42 = sin(qJ(1));
t36 = t42 * pkin(1);
t44 = cos(qJ(1));
t37 = t44 * pkin(1);
t40 = qJ(2) + pkin(5);
t38 = qJ(1) + pkin(7);
t35 = cos(t38);
t34 = sin(t38);
t1 = (-m(2) * pkin(5) - mrSges(4,1) - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) + t54 * (pkin(3) + t40) + t56 * t43 + t55 * t41 + (-m(3) - m(4)) * t40) * g(3) + (-m(3) * t36 - t42 * mrSges(2,1) - t44 * mrSges(2,2) - mrSges(1,2) - t53 * (t34 * pkin(2) + t36) + (t53 * qJ(3) - t50) * t35 + t51 * t34) * g(2) + (-m(3) * t37 - t44 * mrSges(2,1) + t42 * mrSges(2,2) - mrSges(1,1) - t53 * (t35 * pkin(2) + t34 * qJ(3) + t37) + t51 * t35 + t50 * t34) * g(1);
U = t1;
