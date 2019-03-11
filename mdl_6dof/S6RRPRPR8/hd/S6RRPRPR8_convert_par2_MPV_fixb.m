% Return the minimum parameter vector for
% S6RRPRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
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
% Datum: 2019-03-09 10:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RRPRPR8_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR8_convert_par2_MPV_fixb: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR8_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR8_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR8_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t142 = sin(pkin(10));
t143 = cos(pkin(10));
t154 = t142 * t143;
t144 = (pkin(9) ^ 2);
t146 = (pkin(5) ^ 2);
t153 = (Ifges(5,2) + Ifges(6,3) + (t144 + t146) * m(7));
t152 = Ifges(4,4) * t154;
t151 = -pkin(8) * m(5) - mrSges(5,3);
t150 = pkin(9) * m(7) + mrSges(7,3);
t147 = (pkin(3) ^ 2);
t149 = (t147 * m(5) + Ifges(3,2) + Ifges(4,3));
t148 = 2 * pkin(8) * mrSges(5,3) + 2 * pkin(9) * mrSges(7,3) + Ifges(7,2) + t153;
t145 = pkin(8) ^ 2;
t139 = t143 ^ 2;
t138 = t142 ^ 2;
t135 = t151 * pkin(3) + Ifges(4,5);
t134 = t145 * m(5) + Ifges(4,1) + t148;
t133 = Ifges(4,2) + (t145 + t147) * m(5) + t148;
t1 = [Ifges(2,3) + 2 * pkin(7) * mrSges(3,3) + (pkin(1) ^ 2 + pkin(7) ^ 2) * m(3) + t149; m(3) * pkin(1) + mrSges(2,1); -pkin(7) * m(3) + mrSges(2,2) - mrSges(3,3); t138 * t133 + t139 * t134 + Ifges(3,1) - t149 - 0.2e1 * t152; t142 * Ifges(4,6) - t143 * t135 + Ifges(3,4); Ifges(3,5) + (t139 - t138) * Ifges(4,4) + (-t133 + t134) * t154; -t143 * Ifges(4,6) - t142 * t135 + Ifges(3,6); t139 * t133 + t138 * t134 + Ifges(3,3) + 0.2e1 * t152; mrSges(3,1); mrSges(3,2); m(5) * pkin(3) + mrSges(4,1); mrSges(4,2); mrSges(4,3) - t151; m(4) + m(5); m(7) * t144 + Ifges(5,1) + Ifges(6,1) - t153; Ifges(5,4) - Ifges(6,5); t150 * pkin(5) + Ifges(6,4) + Ifges(5,5); Ifges(5,6) - Ifges(6,6); t146 * m(7) + Ifges(6,2) + Ifges(5,3); mrSges(5,1); mrSges(5,2); m(7) * pkin(5) + mrSges(6,1); mrSges(6,2) - t150; mrSges(6,3); m(6) + m(7); Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
