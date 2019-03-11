% Return the minimum parameter vector for
% S6RRPRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
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
% Datum: 2019-03-09 11:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RRPRPR13_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR13_convert_par2_MPV_fixb: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR13_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR13_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR13_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t169 = (m(3) * pkin(8));
t168 = (pkin(9) * mrSges(5,3));
t167 = 2 * pkin(10) * mrSges(7,3) + Ifges(7,2);
t153 = sin(pkin(11));
t155 = cos(pkin(11));
t166 = t153 * t155;
t159 = (pkin(5) ^ 2);
t165 = (t159 * m(7) + Ifges(5,2) + Ifges(6,3));
t164 = Ifges(6,4) * t166;
t163 = -pkin(9) * m(5) - mrSges(5,3);
t162 = -pkin(10) * m(7) - mrSges(7,3);
t161 = (Ifges(3,2) + Ifges(4,3) + t165);
t160 = (pkin(3) ^ 2);
t158 = (pkin(9) ^ 2);
t157 = pkin(10) ^ 2;
t154 = sin(pkin(6));
t152 = 2 * t168;
t150 = t155 ^ 2;
t148 = t153 ^ 2;
t147 = (t158 + t160);
t146 = t162 * pkin(5) + Ifges(6,5);
t145 = m(7) * t157 + Ifges(6,1) + t167;
t144 = Ifges(6,2) + (t157 + t159) * m(7) + t167;
t1 = [pkin(1) ^ 2 * m(3) + Ifges(2,3) + (t147 * m(5) + t152 + (2 * mrSges(3,3) + t169) * pkin(8) + t161) * t154 ^ 2; m(3) * pkin(1) + mrSges(2,1); mrSges(2,2) + (-mrSges(3,3) - t169) * t154; -2 * t168 + Ifges(3,1) + Ifges(4,2) + (-t147 + t160) * m(5) - t161; Ifges(3,4) + Ifges(4,6); t163 * pkin(3) - Ifges(4,4) + Ifges(3,5); Ifges(3,6) - Ifges(4,5); t158 * m(5) + Ifges(4,1) + Ifges(3,3) + t152 + t165; mrSges(3,1); mrSges(3,2); m(5) * pkin(3) + mrSges(4,1); mrSges(4,2) + t163; mrSges(4,3); m(4) + m(5); t148 * t144 + t150 * t145 + Ifges(5,1) - 0.2e1 * t164 - t165; t153 * Ifges(6,6) - t155 * t146 + Ifges(5,4); Ifges(5,5) + (t150 - t148) * Ifges(6,4) + (-t144 + t145) * t166; -t155 * Ifges(6,6) - t153 * t146 + Ifges(5,6); t150 * t144 + t148 * t145 + Ifges(5,3) + 0.2e1 * t164; mrSges(5,1); mrSges(5,2); m(7) * pkin(5) + mrSges(6,1); mrSges(6,2); mrSges(6,3) - t162; m(6) + m(7); Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
