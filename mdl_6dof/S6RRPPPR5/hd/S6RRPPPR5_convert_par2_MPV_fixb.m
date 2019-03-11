% Return the minimum parameter vector for
% S6RRPPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% MPV [29x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RRPPPR5_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR5_convert_par2_MPV_fixb: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR5_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPPR5_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPPR5_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t144 = 2 * pkin(8) * mrSges(7,3) + Ifges(7,2);
t136 = sin(pkin(9));
t137 = cos(pkin(9));
t143 = t136 * t137;
t142 = m(7) * pkin(8) + mrSges(7,3);
t132 = Ifges(4,4) + Ifges(5,6) - Ifges(6,5);
t141 = t132 * t143;
t139 = (pkin(5) ^ 2);
t140 = t139 * m(7) + Ifges(5,1) + Ifges(3,2) + Ifges(6,2) + Ifges(4,3);
t138 = pkin(8) ^ 2;
t134 = t137 ^ 2;
t133 = t136 ^ 2;
t131 = Ifges(4,5) - Ifges(5,4) - Ifges(6,6);
t130 = pkin(5) * t142 + Ifges(6,4) - Ifges(5,5) + Ifges(4,6);
t129 = m(7) * t138 + Ifges(6,1) + Ifges(4,2) + Ifges(5,3) + t144;
t128 = Ifges(4,1) + Ifges(5,2) + Ifges(6,3) + (t138 + t139) * m(7) + t144;
t1 = [Ifges(2,3) + 2 * pkin(7) * mrSges(3,3) + (pkin(1) ^ 2 + pkin(7) ^ 2) * m(3) + t140; m(3) * pkin(1) + mrSges(2,1); -pkin(7) * m(3) + mrSges(2,2) - mrSges(3,3); t128 * t134 + t129 * t133 + Ifges(3,1) - t140 - 0.2e1 * t141; t130 * t136 - t131 * t137 + Ifges(3,4); Ifges(3,5) + (t134 - t133) * t132 + (t128 - t129) * t143; -t130 * t137 - t131 * t136 + Ifges(3,6); t128 * t133 + t129 * t134 + Ifges(3,3) + 0.2e1 * t141; mrSges(3,1); mrSges(3,2); mrSges(4,1); mrSges(4,2); mrSges(4,3); m(4); mrSges(5,1); mrSges(5,2); mrSges(5,3); m(5); m(7) * pkin(5) + mrSges(6,1); mrSges(6,2) - t142; mrSges(6,3); m(6) + m(7); Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
