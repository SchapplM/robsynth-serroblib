% Calculate potential energy for
% S4PRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4PRRR3_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR3_energypot_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR3_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR3_energypot_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR3_energypot_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRR3_energypot_fixb_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:31:32
% EndTime: 2019-12-31 16:31:33
% DurationCPUTime: 0.14s
% Computational Cost: add. (84->38), mult. (63->32), div. (0->0), fcn. (38->8), ass. (0->17)
t48 = -m(4) - m(5);
t40 = sin(qJ(4));
t41 = cos(qJ(4));
t47 = -m(5) * pkin(3) - t41 * mrSges(5,1) + t40 * mrSges(5,2) - mrSges(4,1);
t46 = m(5) * pkin(6) - mrSges(4,2) + mrSges(5,3);
t38 = sin(pkin(7));
t34 = t38 * pkin(1);
t39 = cos(pkin(7));
t35 = t39 * pkin(1);
t45 = pkin(4) + qJ(1);
t37 = pkin(7) + qJ(2);
t33 = qJ(3) + t37;
t32 = cos(t37);
t31 = sin(t37);
t30 = cos(t33);
t29 = sin(t33);
t1 = (-m(2) * qJ(1) - m(3) * t45 - t40 * mrSges(5,1) - t41 * mrSges(5,2) - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) - mrSges(4,3) + t48 * (pkin(5) + t45)) * g(3) + (-m(3) * t34 - t38 * mrSges(2,1) - t31 * mrSges(3,1) - t39 * mrSges(2,2) - t32 * mrSges(3,2) - mrSges(1,2) + t48 * (pkin(2) * t31 + t34) + t46 * t30 + t47 * t29) * g(2) + (-m(3) * t35 - t39 * mrSges(2,1) - t32 * mrSges(3,1) + t38 * mrSges(2,2) + t31 * mrSges(3,2) - mrSges(1,1) + t48 * (pkin(2) * t32 + t35) + t47 * t30 - t46 * t29) * g(1);
U = t1;
