% Return the minimum parameter vector for
% S6RRRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
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
% Datum: 2019-03-09 15:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RRRPPR4_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR4_convert_par2_MPV_fixb: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR4_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR4_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR4_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t155 = 2 * pkin(9) * mrSges(7,3) + Ifges(7,2);
t139 = sin(pkin(10));
t140 = cos(pkin(10));
t154 = t139 * t140;
t144 = pkin(9) ^ 2;
t146 = (pkin(5) ^ 2);
t131 = Ifges(5,2) + Ifges(6,3) + (t144 + t146) * m(7) + t155;
t133 = m(7) * t144 + Ifges(5,1) + Ifges(6,1) + t155;
t135 = t139 ^ 2;
t136 = t140 ^ 2;
t153 = t136 * t131 + t135 * t133 + Ifges(4,2);
t152 = m(4) * pkin(8) + mrSges(4,3);
t151 = m(7) * pkin(9) + mrSges(7,3);
t142 = Ifges(5,4) - Ifges(6,5);
t150 = t142 * t154;
t149 = (2 * pkin(8) * mrSges(4,3)) + 0.2e1 * t150 + t153;
t148 = t140 * mrSges(5,1) - t139 * mrSges(5,2);
t147 = (pkin(2) ^ 2);
t145 = pkin(8) ^ 2;
t143 = (m(3) + m(4));
t141 = Ifges(5,6) - Ifges(6,6);
t134 = t151 * pkin(5) + Ifges(6,4) + Ifges(5,5);
t1 = [Ifges(2,3) + t147 * m(4) + Ifges(3,2) + 2 * pkin(7) * mrSges(3,3) + (pkin(1) ^ 2 + pkin(7) ^ 2) * t143; pkin(1) * t143 + mrSges(2,1); -pkin(7) * t143 + mrSges(2,2) - mrSges(3,3); Ifges(3,1) - Ifges(3,2) + ((t145 - t147) * m(4)) + t149; t152 * pkin(2) + Ifges(3,4); Ifges(3,5); Ifges(3,6); Ifges(3,3) + ((t145 + t147) * m(4)) + t149; m(4) * pkin(2) + mrSges(3,1); mrSges(3,2) - t152; t135 * t131 + t136 * t133 + Ifges(4,1) - 0.4e1 * t150 - t153; Ifges(4,4) + (t136 - t135) * t142 + (-t131 + t133) * t154; t140 * t134 - t139 * t141 + Ifges(4,5); t139 * t134 + t140 * t141 + Ifges(4,6); (t146 * m(7)) + 0.2e1 * pkin(3) * t148 + Ifges(6,2) + Ifges(4,3) + Ifges(5,3); mrSges(4,1) + t148; t139 * mrSges(5,1) + t140 * mrSges(5,2) + mrSges(4,2); mrSges(5,3); m(5); m(7) * pkin(5) + mrSges(6,1); mrSges(6,2) - t151; mrSges(6,3); m(6) + m(7); Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
