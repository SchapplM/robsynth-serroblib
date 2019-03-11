% Return the minimum parameter vector for
% S6RRRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,theta5]';
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
% Datum: 2019-03-09 21:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RRRRPP4_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP4_convert_par2_MPV_fixb: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP4_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP4_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP4_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t140 = (m(4) + m(5));
t141 = pkin(9) ^ 2;
t143 = pkin(3) ^ 2;
t152 = Ifges(4,2) + (t141 + t143) * m(5);
t133 = sin(pkin(10));
t134 = cos(pkin(10));
t151 = t133 * t134;
t129 = t133 ^ 2;
t130 = t134 ^ 2;
t136 = Ifges(6,2) + Ifges(7,3);
t139 = Ifges(6,1) + Ifges(7,1);
t150 = t129 * t139 + t130 * t136 + Ifges(5,2);
t149 = -pkin(9) * m(5) - mrSges(5,3);
t138 = Ifges(6,4) - Ifges(7,5);
t148 = t138 * t151;
t127 = (mrSges(4,3) - t149);
t147 = pkin(8) * t140 + t127;
t146 = t134 * mrSges(6,1) - t133 * mrSges(6,2);
t145 = (2 * pkin(9) * mrSges(5,3)) + (2 * pkin(8) * t127) + 0.2e1 * t148 + t150 + t152;
t144 = (pkin(2) ^ 2);
t142 = pkin(8) ^ 2;
t137 = Ifges(6,5) + Ifges(7,4);
t135 = Ifges(6,6) - Ifges(7,6);
t131 = (m(3) + t140);
t1 = [Ifges(2,3) + Ifges(3,2) + t144 * t140 + 2 * pkin(7) * mrSges(3,3) + (pkin(1) ^ 2 + pkin(7) ^ 2) * t131; pkin(1) * t131 + mrSges(2,1); -pkin(7) * t131 + mrSges(2,2) - mrSges(3,3); Ifges(3,1) - Ifges(3,2) + ((t142 - t144) * t140) + t145; t147 * pkin(2) + Ifges(3,4); Ifges(3,5); Ifges(3,6); Ifges(3,3) + ((t142 + t144) * t140) + t145; pkin(2) * t140 + mrSges(3,1); mrSges(3,2) - t147; t141 * m(5) + Ifges(4,1) - t152; Ifges(4,4); t149 * pkin(3) + Ifges(4,5); Ifges(4,6); t143 * m(5) + Ifges(4,3); m(5) * pkin(3) + mrSges(4,1); mrSges(4,2); t129 * t136 + t130 * t139 + Ifges(5,1) - 0.4e1 * t148 - t150; Ifges(5,4) + (t130 - t129) * t138 + (-t136 + t139) * t151; -t133 * t135 + t134 * t137 + Ifges(5,5); t133 * t137 + t134 * t135 + Ifges(5,6); 0.2e1 * pkin(4) * t146 + Ifges(7,2) + Ifges(5,3) + Ifges(6,3); mrSges(5,1) + t146; t133 * mrSges(6,1) + t134 * mrSges(6,2) + mrSges(5,2); mrSges(6,3); m(6); mrSges(7,1); mrSges(7,2); mrSges(7,3); m(7);];
MPV  = t1;
