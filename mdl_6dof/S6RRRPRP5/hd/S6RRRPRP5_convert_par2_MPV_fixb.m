% Return the minimum parameter vector for
% S6RRRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
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
% MPV [30x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RRRPRP5_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP5_convert_par2_MPV_fixb: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP5_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP5_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP5_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t147 = (-Ifges(6,2) - Ifges(7,3));
t132 = sin(pkin(10));
t133 = cos(pkin(10));
t146 = t132 * t133;
t145 = 2 * pkin(9) * mrSges(6,3) - t147;
t144 = Ifges(5,4) * t146;
t135 = pkin(9) ^ 2;
t137 = (pkin(4) ^ 2);
t123 = Ifges(5,2) + (t135 + t137) * m(6) + t145;
t124 = m(6) * t135 + Ifges(5,1) + t145;
t128 = t132 ^ 2;
t129 = t133 ^ 2;
t143 = t123 * t129 + t124 * t128 + Ifges(4,2);
t142 = m(4) * pkin(8) + mrSges(4,3);
t141 = -m(6) * pkin(9) - mrSges(6,3);
t140 = (2 * pkin(8) * mrSges(4,3)) + t143 + 0.2e1 * t144;
t127 = m(6) * pkin(4) + mrSges(5,1);
t139 = -mrSges(5,2) * t132 + t127 * t133;
t138 = (pkin(2) ^ 2);
t136 = pkin(8) ^ 2;
t134 = (m(3) + m(4));
t126 = pkin(4) * t141 + Ifges(5,5);
t1 = [Ifges(2,3) + t138 * m(4) + Ifges(3,2) + 2 * pkin(7) * mrSges(3,3) + (pkin(1) ^ 2 + pkin(7) ^ 2) * t134; pkin(1) * t134 + mrSges(2,1); -pkin(7) * t134 + mrSges(2,2) - mrSges(3,3); Ifges(3,1) - Ifges(3,2) + ((t136 - t138) * m(4)) + t140; pkin(2) * t142 + Ifges(3,4); Ifges(3,5); Ifges(3,6); Ifges(3,3) + ((t136 + t138) * m(4)) + t140; m(4) * pkin(2) + mrSges(3,1); mrSges(3,2) - t142; t123 * t128 + t124 * t129 + Ifges(4,1) - t143 - 0.4e1 * t144; Ifges(4,4) + (t129 - t128) * Ifges(5,4) + (-t123 + t124) * t146; -Ifges(5,6) * t132 + t126 * t133 + Ifges(4,5); Ifges(5,6) * t133 + t126 * t132 + Ifges(4,6); (t137 * m(6)) + 0.2e1 * pkin(3) * t139 + Ifges(4,3) + Ifges(5,3); mrSges(4,1) + t139; mrSges(5,2) * t133 + t127 * t132 + mrSges(4,2); mrSges(5,3) - t141; m(5) + m(6); Ifges(6,1) + Ifges(7,1) + t147; Ifges(6,4) - Ifges(7,5); Ifges(6,5) + Ifges(7,4); Ifges(6,6) - Ifges(7,6); Ifges(6,3) + Ifges(7,2); mrSges(6,1); mrSges(6,2); mrSges(7,1); mrSges(7,2); mrSges(7,3); m(7);];
MPV  = t1;
