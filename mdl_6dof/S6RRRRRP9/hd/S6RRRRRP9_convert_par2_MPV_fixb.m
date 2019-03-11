% Return the minimum parameter vector for
% S6RRRRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
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
% MPV [33x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 02:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RRRRRP9_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP9_convert_par2_MPV_fixb: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP9_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP9_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP9_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t151 = (m(5) + m(6));
t147 = (m(4) + t151);
t143 = m(3) + t147;
t166 = (t143 * pkin(8));
t165 = (-Ifges(6,2) - Ifges(7,2));
t156 = (pkin(3) ^ 2);
t164 = (t156 * t151 + Ifges(4,2));
t152 = (pkin(11) ^ 2);
t155 = (pkin(4) ^ 2);
t163 = (Ifges(5,2) + (t152 + t155) * m(6));
t162 = 2 * pkin(9) * mrSges(4,3) + t164;
t161 = -pkin(11) * m(6) - mrSges(6,3);
t160 = pkin(9) * t147 + mrSges(4,3);
t142 = (mrSges(5,3) - t161);
t159 = pkin(10) * t151 + t142;
t158 = 2 * pkin(11) * mrSges(6,3) + 2 * pkin(10) * t142 + t163 - t165;
t157 = (pkin(2) ^ 2);
t154 = pkin(9) ^ 2;
t153 = pkin(10) ^ 2;
t150 = sin(pkin(6));
t1 = [pkin(1) ^ 2 * t143 + Ifges(2,3) + (t157 * t147 + Ifges(3,2) + (2 * mrSges(3,3) + t166) * pkin(8)) * t150 ^ 2; pkin(1) * t143 + mrSges(2,1); mrSges(2,2) + (-mrSges(3,3) - t166) * t150; Ifges(3,1) - Ifges(3,2) + (t154 - t157) * t147 + t162; t160 * pkin(2) + Ifges(3,4); Ifges(3,5); Ifges(3,6); Ifges(3,3) + (t154 + t157) * t147 + t162; pkin(2) * t147 + mrSges(3,1); mrSges(3,2) - t160; t153 * t151 + Ifges(4,1) + t158 - t164; t159 * pkin(3) + Ifges(4,4); Ifges(4,5); Ifges(4,6); Ifges(4,3) + (t153 + t156) * t151 + t158; pkin(3) * t151 + mrSges(4,1); mrSges(4,2) - t159; t152 * m(6) + Ifges(5,1) - t163; Ifges(5,4); t161 * pkin(4) + Ifges(5,5); Ifges(5,6); t155 * m(6) + Ifges(5,3); m(6) * pkin(4) + mrSges(5,1); mrSges(5,2); Ifges(6,1) + Ifges(7,1) + t165; Ifges(6,4) + Ifges(7,4); Ifges(6,5) + Ifges(7,5); Ifges(6,6) + Ifges(7,6); 2 * pkin(5) * mrSges(7,1) + Ifges(6,3) + Ifges(7,3); mrSges(6,1) + mrSges(7,1); mrSges(6,2) + mrSges(7,2); mrSges(7,3); m(7);];
MPV  = t1;
