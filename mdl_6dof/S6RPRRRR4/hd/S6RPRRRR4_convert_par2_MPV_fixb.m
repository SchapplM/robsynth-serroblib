% Return the minimum parameter vector for
% S6RPRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
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
% MPV [35x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RPRRRR4_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR4_convert_par2_MPV_fixb: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR4_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR4_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR4_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t127 = (m(6) + m(7));
t144 = -pkin(9) * t127 - mrSges(6,3);
t116 = (mrSges(5,3) - t144);
t122 = (m(5) + t127);
t143 = -pkin(8) * t122 - t116;
t131 = (pkin(8) ^ 2);
t135 = (pkin(3) ^ 2);
t142 = (Ifges(4,2) + (t131 + t135) * t122);
t130 = (pkin(9) ^ 2);
t134 = (pkin(4) ^ 2);
t141 = (Ifges(5,2) + (t130 + t134) * t127);
t133 = (pkin(5) ^ 2);
t140 = (m(7) * t133 + Ifges(6,2));
t139 = 2 * pkin(10) * mrSges(7,3) + Ifges(7,2);
t118 = (m(4) + t122);
t138 = m(7) * pkin(10) + mrSges(7,3);
t137 = (mrSges(4,3) - t143);
t136 = 2 * pkin(9) * mrSges(6,3) + 2 * pkin(7) * t137 + 2 * pkin(8) * t116 + t140 + t141 + t142;
t132 = pkin(7) ^ 2;
t129 = pkin(10) ^ 2;
t126 = cos(pkin(11));
t125 = sin(pkin(11));
t1 = [Ifges(2,3) + t126 ^ 2 * (Ifges(3,2) + (pkin(2) ^ 2 + t132) * t118 + t136) + (0.2e1 * t126 * Ifges(3,4) + (t118 * t132 + Ifges(3,1) + t136) * t125) * t125; mrSges(2,1); mrSges(2,2); pkin(2) * t118 + mrSges(3,1); mrSges(3,2); pkin(7) * t118 + mrSges(3,3) + t137; m(3) + t118; t122 * t131 + Ifges(4,1) - t142; Ifges(4,4); pkin(3) * t143 + Ifges(4,5); Ifges(4,6); t122 * t135 + Ifges(4,3); pkin(3) * t122 + mrSges(4,1); mrSges(4,2); t127 * t130 + Ifges(5,1) - t141; Ifges(5,4); pkin(4) * t144 + Ifges(5,5); Ifges(5,6); t127 * t134 + Ifges(5,3); pkin(4) * t127 + mrSges(5,1); mrSges(5,2); m(7) * t129 + Ifges(6,1) + t139 - t140; pkin(5) * t138 + Ifges(6,4); Ifges(6,5); Ifges(6,6); Ifges(6,3) + (t129 + t133) * m(7) + t139; m(7) * pkin(5) + mrSges(6,1); mrSges(6,2) - t138; Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
