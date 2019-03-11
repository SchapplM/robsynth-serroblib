% Return the minimum parameter vector for
% S6PRRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1,theta5]';
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
% MPV [24x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6PRRRPP1_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPP1_convert_par2_MPV_fixb: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPP1_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPP1_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPP1_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t152 = (m(4) + m(5));
t145 = sin(pkin(11));
t146 = cos(pkin(11));
t160 = t145 * t146;
t142 = t145 ^ 2;
t143 = t146 ^ 2;
t148 = Ifges(6,2) + Ifges(7,3);
t151 = Ifges(6,1) + Ifges(7,1);
t159 = t142 * t151 + t143 * t148 + Ifges(5,2);
t158 = pkin(9) * m(5) + mrSges(5,3);
t150 = Ifges(6,4) - Ifges(7,5);
t157 = t150 * t160;
t156 = (2 * pkin(9) * mrSges(5,3)) + 0.2e1 * t157 + t159;
t155 = t146 * mrSges(6,1) - t145 * mrSges(6,2);
t154 = (pkin(3) ^ 2);
t153 = pkin(9) ^ 2;
t149 = Ifges(6,5) + Ifges(7,4);
t147 = Ifges(6,6) - Ifges(7,6);
t1 = [m(2) + m(3) + t152; Ifges(3,3) + t154 * m(5) + Ifges(4,2) + 2 * pkin(8) * mrSges(4,3) + (pkin(2) ^ 2 + pkin(8) ^ 2) * t152; pkin(2) * t152 + mrSges(3,1); -pkin(8) * t152 + mrSges(3,2) - mrSges(4,3); Ifges(4,1) - Ifges(4,2) + ((t153 - t154) * m(5)) + t156; t158 * pkin(3) + Ifges(4,4); Ifges(4,5); Ifges(4,6); Ifges(4,3) + ((t153 + t154) * m(5)) + t156; m(5) * pkin(3) + mrSges(4,1); mrSges(4,2) - t158; t142 * t148 + t143 * t151 + Ifges(5,1) - 0.4e1 * t157 - t159; Ifges(5,4) + (t143 - t142) * t150 + (-t148 + t151) * t160; -t145 * t147 + t146 * t149 + Ifges(5,5); t145 * t149 + t146 * t147 + Ifges(5,6); 0.2e1 * pkin(4) * t155 + Ifges(7,2) + Ifges(5,3) + Ifges(6,3); mrSges(5,1) + t155; t145 * mrSges(6,1) + t146 * mrSges(6,2) + mrSges(5,2); mrSges(6,3); m(6); mrSges(7,1); mrSges(7,2); mrSges(7,3); m(7);];
MPV  = t1;
