% Calculate potential energy for
% S5RPPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
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
% Datum: 2019-12-31 18:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPPRR11_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR11_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR11_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR11_energypot_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR11_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR11_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:05:29
% EndTime: 2019-12-31 18:05:29
% DurationCPUTime: 0.27s
% Computational Cost: add. (80->48), mult. (116->35), div. (0->0), fcn. (89->6), ass. (0->19)
t38 = sin(qJ(5));
t41 = cos(qJ(5));
t59 = -m(6) * pkin(4) - mrSges(6,1) * t41 + mrSges(6,2) * t38 - mrSges(5,1);
t58 = -m(6) * pkin(7) + mrSges(5,2) - mrSges(6,3);
t57 = m(5) + m(6);
t39 = sin(qJ(4));
t42 = cos(qJ(4));
t55 = t39 * t59 - t42 * t58 - mrSges(2,1) + mrSges(3,2) - mrSges(4,3);
t54 = t38 * mrSges(6,1) + t41 * mrSges(6,2) + mrSges(2,2) - mrSges(4,2) - mrSges(3,3) + mrSges(5,3);
t53 = pkin(2) + pkin(5);
t40 = sin(qJ(1));
t43 = cos(qJ(1));
t52 = t43 * pkin(1) + t40 * qJ(2);
t35 = t40 * pkin(1);
t49 = -t43 * qJ(2) + t35;
t32 = t40 * qJ(3);
t48 = t32 + t49;
t36 = t43 * pkin(6);
t1 = (-m(4) * t53 - mrSges(3,1) - mrSges(4,1) - mrSges(1,3) - mrSges(2,3) - t57 * (pkin(3) + t53) + t59 * t42 + t58 * t39 + (-m(2) - m(3)) * pkin(5)) * g(3) + (-mrSges(1,2) - m(3) * t49 - m(4) * t48 - m(5) * (t36 + t48) - m(6) * (t32 + t35 + t36) + (m(6) * qJ(2) - t54) * t43 + t55 * t40) * g(2) + (-m(3) * t52 - mrSges(1,1) + (-m(4) - t57) * (t43 * qJ(3) + t52) + t55 * t43 + (pkin(6) * t57 + t54) * t40) * g(1);
U = t1;
