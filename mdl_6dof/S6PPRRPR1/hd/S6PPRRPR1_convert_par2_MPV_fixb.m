% Return the minimum parameter vector for
% S6PPRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2,theta5]';
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
% MPV [23x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6PPRRPR1_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRPR1_convert_par2_MPV_fixb: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRPR1_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRPR1_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRRPR1_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t171 = 2 * pkin(10) * mrSges(7,3) + Ifges(7,2);
t163 = sin(pkin(13));
t164 = cos(pkin(13));
t170 = t163 * t164;
t169 = Ifges(6,4) * t170;
t168 = -pkin(10) * m(7) - mrSges(7,3);
t166 = (pkin(5) ^ 2);
t167 = t166 * m(7) + Ifges(5,2) + Ifges(6,3);
t165 = pkin(10) ^ 2;
t161 = t164 ^ 2;
t160 = t163 ^ 2;
t159 = t168 * pkin(5) + Ifges(6,5);
t158 = m(7) * t165 + Ifges(6,1) + t171;
t157 = Ifges(6,2) + (t165 + t166) * m(7) + t171;
t1 = [m(2); m(3) + m(4) + m(5); Ifges(4,3) + 2 * pkin(9) * mrSges(5,3) + (pkin(3) ^ 2 + pkin(9) ^ 2) * m(5) + t167; m(5) * pkin(3) + mrSges(4,1); -pkin(9) * m(5) + mrSges(4,2) - mrSges(5,3); t160 * t157 + t161 * t158 + Ifges(5,1) - t167 - 0.2e1 * t169; t163 * Ifges(6,6) - t164 * t159 + Ifges(5,4); Ifges(5,5) + (t161 - t160) * Ifges(6,4) + (-t157 + t158) * t170; -t164 * Ifges(6,6) - t163 * t159 + Ifges(5,6); t161 * t157 + t160 * t158 + Ifges(5,3) + 0.2e1 * t169; mrSges(5,1); mrSges(5,2); m(7) * pkin(5) + mrSges(6,1); mrSges(6,2); mrSges(6,3) - t168; m(6) + m(7); Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
