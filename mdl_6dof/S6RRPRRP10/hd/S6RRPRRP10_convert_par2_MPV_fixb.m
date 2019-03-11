% Return the minimum parameter vector for
% S6RRPRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-03-09 12:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RRPRRP10_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP10_convert_par2_MPV_fixb: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP10_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP10_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP10_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t177 = (m(3) * pkin(8));
t162 = (m(5) + m(6));
t176 = (-Ifges(6,2) - Ifges(7,3));
t166 = (pkin(4) ^ 2);
t175 = (t166 * m(6) + Ifges(5,2));
t159 = sin(pkin(11));
t161 = cos(pkin(11));
t174 = t159 * t161;
t173 = 2 * pkin(10) * mrSges(6,3) - t176;
t172 = Ifges(4,4) * t174;
t171 = 2 * pkin(9) * mrSges(5,3) + t175;
t170 = pkin(10) * m(6) + mrSges(6,3);
t169 = -pkin(9) * t162 - mrSges(5,3);
t167 = (pkin(3) ^ 2);
t168 = t167 * t162 + Ifges(3,2) + Ifges(4,3);
t165 = pkin(9) ^ 2;
t164 = pkin(10) ^ 2;
t160 = sin(pkin(6));
t156 = t161 ^ 2;
t154 = t159 ^ 2;
t153 = t169 * pkin(3) + Ifges(4,5);
t152 = t165 * t162 + Ifges(4,1) + t171;
t151 = Ifges(4,2) + (t165 + t167) * t162 + t171;
t1 = [pkin(1) ^ 2 * m(3) + Ifges(2,3) + ((2 * mrSges(3,3) + t177) * pkin(8) + t168) * t160 ^ 2; m(3) * pkin(1) + mrSges(2,1); mrSges(2,2) + (-mrSges(3,3) - t177) * t160; t154 * t151 + t156 * t152 + Ifges(3,1) - t168 - 0.2e1 * t172; t159 * Ifges(4,6) - t161 * t153 + Ifges(3,4); Ifges(3,5) + (t156 - t154) * Ifges(4,4) + (-t151 + t152) * t174; -t161 * Ifges(4,6) - t159 * t153 + Ifges(3,6); t156 * t151 + t154 * t152 + Ifges(3,3) + 0.2e1 * t172; mrSges(3,1); mrSges(3,2); pkin(3) * t162 + mrSges(4,1); mrSges(4,2); mrSges(4,3) - t169; m(4) + t162; t164 * m(6) + Ifges(5,1) + t173 - t175; t170 * pkin(4) + Ifges(5,4); Ifges(5,5); Ifges(5,6); Ifges(5,3) + (t164 + t166) * m(6) + t173; m(6) * pkin(4) + mrSges(5,1); mrSges(5,2) - t170; Ifges(6,1) + Ifges(7,1) + t176; Ifges(6,4) - Ifges(7,5); Ifges(6,5) + Ifges(7,4); Ifges(6,6) - Ifges(7,6); Ifges(6,3) + Ifges(7,2); mrSges(6,1); mrSges(6,2); mrSges(7,1); mrSges(7,2); mrSges(7,3); m(7);];
MPV  = t1;
