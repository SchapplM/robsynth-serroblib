% Return the minimum parameter vector for
% S6PRRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-03-08 21:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6PRRPRP2_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP2_convert_par2_MPV_fixb: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP2_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRP2_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRP2_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t150 = (pkin(9) ^ 2);
t158 = (-Ifges(6,2) - Ifges(7,3));
t155 = 2 * pkin(9) * mrSges(6,3) - t158;
t140 = t150 * m(6) + Ifges(5,1) + t155;
t151 = (pkin(4) ^ 2);
t143 = t151 * m(6) + Ifges(5,2);
t159 = t140 - t143;
t148 = sin(pkin(11));
t149 = cos(pkin(11));
t157 = t148 * t149;
t145 = t148 ^ 2;
t146 = t149 ^ 2;
t156 = t146 - t145;
t154 = pkin(9) * m(6) + mrSges(6,3);
t141 = t154 * pkin(4) + Ifges(5,4);
t153 = t141 * t157;
t142 = mrSges(5,2) - t154;
t144 = m(6) * pkin(4) + mrSges(5,1);
t152 = -t148 * t142 + t149 * t144;
t1 = [m(2) + m(3) + m(4); Ifges(3,3) + Ifges(4,2) + t145 * t140 + 0.2e1 * t153 + t146 * t143 + (2 * pkin(8) * mrSges(4,3)) + ((pkin(2) ^ 2 + pkin(8) ^ 2) * m(4)); m(4) * pkin(2) + mrSges(3,1); -pkin(8) * m(4) + mrSges(3,2) - mrSges(4,3); t159 * t156 + Ifges(4,1) - Ifges(4,2) - 0.4e1 * t153; t156 * t141 + t159 * t157 + Ifges(4,4); t149 * Ifges(5,5) - t148 * Ifges(5,6) + Ifges(4,5); t148 * Ifges(5,5) + t149 * Ifges(5,6) + Ifges(4,6); Ifges(4,3) + Ifges(5,3) + ((t150 + t151) * m(6)) + 0.2e1 * t152 * pkin(3) + t155; mrSges(4,1) + t152; t149 * t142 + t148 * t144 + mrSges(4,2); mrSges(5,3); m(5) + m(6); Ifges(6,1) + Ifges(7,1) + t158; Ifges(6,4) - Ifges(7,5); Ifges(6,5) + Ifges(7,4); Ifges(6,6) - Ifges(7,6); Ifges(6,3) + Ifges(7,2); mrSges(6,1); mrSges(6,2); mrSges(7,1); mrSges(7,2); mrSges(7,3); m(7);];
MPV  = t1;
