% Return the minimum parameter vector for
% S6RRRPPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta5]';
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
% Datum: 2019-03-09 16:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RRRPPR10_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR10_convert_par2_MPV_fixb: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR10_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR10_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR10_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t161 = m(3) + m(4);
t173 = (t161 * pkin(8));
t172 = 2 * pkin(10) * mrSges(7,3) + Ifges(7,2);
t158 = sin(pkin(11));
t160 = cos(pkin(11));
t171 = t158 * t160;
t170 = Ifges(6,4) * t171;
t169 = pkin(9) * m(4) + mrSges(4,3);
t168 = -pkin(10) * m(7) - mrSges(7,3);
t162 = pkin(10) ^ 2;
t164 = (pkin(5) ^ 2);
t148 = Ifges(6,2) + (t162 + t164) * m(7) + t172;
t151 = m(7) * t162 + Ifges(6,1) + t172;
t153 = t158 ^ 2;
t155 = t160 ^ 2;
t167 = t155 * t148 + t153 * t151 + Ifges(4,2) + Ifges(5,3);
t166 = (2 * pkin(9) * mrSges(4,3)) + t167 + 0.2e1 * t170;
t165 = (pkin(2) ^ 2);
t163 = pkin(9) ^ 2;
t159 = sin(pkin(6));
t152 = t168 * pkin(5) + Ifges(6,5);
t150 = -0.2e1 * t170;
t1 = [pkin(1) ^ 2 * t161 + Ifges(2,3) + (t165 * m(4) + Ifges(3,2) + (2 * mrSges(3,3) + t173) * pkin(8)) * t159 ^ 2; pkin(1) * t161 + mrSges(2,1); mrSges(2,2) + (-mrSges(3,3) - t173) * t159; Ifges(3,1) - Ifges(3,2) + ((t163 - t165) * m(4)) + t166; t169 * pkin(2) + Ifges(3,4); Ifges(3,5); Ifges(3,6); Ifges(3,3) + ((t163 + t165) * m(4)) + t166; m(4) * pkin(2) + mrSges(3,1); mrSges(3,2) - t169; (t164 * m(7)) + Ifges(4,1) + Ifges(5,2) + Ifges(6,3) + t150 - t167; -t160 * Ifges(6,6) - t158 * t152 + Ifges(4,4) + Ifges(5,6); -t158 * Ifges(6,6) + t160 * t152 - Ifges(5,4) + Ifges(4,5); Ifges(4,6) - Ifges(5,5) - (t155 - t153) * Ifges(6,4) + (t148 - t151) * t171; t153 * t148 + t155 * t151 + Ifges(5,1) + Ifges(4,3) + t150; mrSges(4,1); mrSges(4,2); mrSges(5,1); mrSges(5,2); mrSges(5,3); m(5); m(7) * pkin(5) + mrSges(6,1); mrSges(6,2); mrSges(6,3) - t168; m(6) + m(7); Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
