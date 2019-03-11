% Return the minimum parameter vector for
% S6RRRPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta5]';
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
% MPV [32x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RRRPPR7_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR7_convert_par2_MPV_fixb: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR7_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR7_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR7_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t153 = 2 * pkin(9) * mrSges(7,3) + Ifges(7,2);
t140 = sin(pkin(10));
t141 = cos(pkin(10));
t152 = t140 * t141;
t151 = Ifges(6,4) * t152;
t150 = pkin(8) * m(4) + mrSges(4,3);
t149 = -pkin(9) * m(7) - mrSges(7,3);
t143 = pkin(9) ^ 2;
t145 = (pkin(5) ^ 2);
t132 = Ifges(6,2) + (t143 + t145) * m(7) + t153;
t134 = m(7) * t143 + Ifges(6,1) + t153;
t136 = t140 ^ 2;
t137 = t141 ^ 2;
t148 = t137 * t132 + t136 * t134 + Ifges(4,2) + Ifges(5,3);
t147 = (2 * pkin(8) * mrSges(4,3)) + t148 + 0.2e1 * t151;
t146 = (pkin(2) ^ 2);
t144 = pkin(8) ^ 2;
t142 = (m(3) + m(4));
t135 = t149 * pkin(5) + Ifges(6,5);
t1 = [Ifges(2,3) + t146 * m(4) + Ifges(3,2) + 2 * pkin(7) * mrSges(3,3) + (pkin(1) ^ 2 + pkin(7) ^ 2) * t142; pkin(1) * t142 + mrSges(2,1); -pkin(7) * t142 + mrSges(2,2) - mrSges(3,3); Ifges(3,1) - Ifges(3,2) + ((t144 - t146) * m(4)) + t147; t150 * pkin(2) + Ifges(3,4); Ifges(3,5); Ifges(3,6); Ifges(3,3) + ((t144 + t146) * m(4)) + t147; m(4) * pkin(2) + mrSges(3,1); mrSges(3,2) - t150; t136 * t132 + t137 * t134 + Ifges(4,1) + Ifges(5,1) - t148 - 0.4e1 * t151; Ifges(4,4) - Ifges(5,5) - (t137 - t136) * Ifges(6,4) + (t132 - t134) * t152; t140 * Ifges(6,6) - t141 * t135 + Ifges(5,4) + Ifges(4,5); t141 * Ifges(6,6) + t140 * t135 + Ifges(4,6) - Ifges(5,6); t145 * m(7) + Ifges(5,2) + Ifges(4,3) + Ifges(6,3); mrSges(4,1); mrSges(4,2); mrSges(5,1); mrSges(5,2); mrSges(5,3); m(5); m(7) * pkin(5) + mrSges(6,1); mrSges(6,2); mrSges(6,3) - t149; m(6) + m(7); Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
