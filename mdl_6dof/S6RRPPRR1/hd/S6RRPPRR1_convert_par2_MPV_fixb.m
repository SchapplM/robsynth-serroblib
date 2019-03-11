% Return the minimum parameter vector for
% S6RRPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
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
% Datum: 2019-03-09 08:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RRPPRR1_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR1_convert_par2_MPV_fixb: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR1_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR1_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR1_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t139 = (m(6) + m(7));
t142 = (pkin(8) ^ 2);
t144 = (pkin(4) ^ 2);
t143 = (pkin(5) ^ 2);
t153 = (t143 * m(7) + Ifges(6,2));
t149 = 2 * pkin(8) * mrSges(6,3) + t153;
t128 = Ifges(4,2) + Ifges(5,3) + (t142 + t144) * t139 + t149;
t129 = t142 * t139 + Ifges(4,1) + Ifges(5,1) + t149;
t154 = -t128 + t129;
t152 = 2 * pkin(9) * mrSges(7,3) + Ifges(7,2);
t135 = sin(pkin(10));
t136 = cos(pkin(10));
t151 = t135 * t136;
t131 = t135 ^ 2;
t132 = t136 ^ 2;
t150 = t132 - t131;
t148 = pkin(9) * m(7) + mrSges(7,3);
t138 = Ifges(4,4) - Ifges(5,5);
t147 = t138 * t151;
t146 = pkin(8) * t139 + mrSges(6,3);
t145 = t136 * mrSges(4,1) - t135 * mrSges(4,2);
t141 = pkin(9) ^ 2;
t137 = Ifges(4,6) - Ifges(5,6);
t130 = t146 * pkin(4) + Ifges(5,4) + Ifges(4,5);
t1 = [Ifges(2,3) + Ifges(3,2) + t131 * t129 + 0.2e1 * t147 + t132 * t128 + (2 * pkin(7) * mrSges(3,3)) + ((pkin(1) ^ 2 + pkin(7) ^ 2) * m(3)); m(3) * pkin(1) + mrSges(2,1); -pkin(7) * m(3) + mrSges(2,2) - mrSges(3,3); t154 * t150 + Ifges(3,1) - Ifges(3,2) - 0.4e1 * t147; t150 * t138 + t154 * t151 + Ifges(3,4); t136 * t130 - t135 * t137 + Ifges(3,5); t135 * t130 + t136 * t137 + Ifges(3,6); 0.2e1 * pkin(2) * t145 + (t144 * t139) + Ifges(5,2) + Ifges(3,3) + Ifges(4,3); mrSges(3,1) + t145; t135 * mrSges(4,1) + t136 * mrSges(4,2) + mrSges(3,2); mrSges(4,3); m(4); pkin(4) * t139 + mrSges(5,1); mrSges(5,2) - t146; mrSges(5,3); m(5) + t139; m(7) * t141 + Ifges(6,1) + t152 - t153; t148 * pkin(5) + Ifges(6,4); Ifges(6,5); Ifges(6,6); Ifges(6,3) + (t141 + t143) * m(7) + t152; m(7) * pkin(5) + mrSges(6,1); mrSges(6,2) - t148; Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
