% Return the minimum parameter vector for
% S6PRRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
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
% MPV [27x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6PRRPRR1_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR1_convert_par2_MPV_fixb: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR1_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR1_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR1_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t145 = (m(6) + m(7));
t148 = (pkin(9) ^ 2);
t150 = (pkin(4) ^ 2);
t149 = (pkin(5) ^ 2);
t159 = (t149 * m(7) + Ifges(6,2));
t154 = 2 * pkin(9) * mrSges(6,3) + t159;
t135 = Ifges(5,2) + (t148 + t150) * t145 + t154;
t136 = t148 * t145 + Ifges(5,1) + t154;
t160 = -t135 + t136;
t158 = 2 * pkin(10) * mrSges(7,3) + Ifges(7,2);
t143 = sin(pkin(12));
t144 = cos(pkin(12));
t157 = t143 * t144;
t139 = t143 ^ 2;
t140 = t144 ^ 2;
t156 = t140 - t139;
t155 = Ifges(5,4) * t157;
t153 = pkin(10) * m(7) + mrSges(7,3);
t152 = -pkin(9) * t145 - mrSges(6,3);
t138 = pkin(4) * t145 + mrSges(5,1);
t151 = -t143 * mrSges(5,2) + t144 * t138;
t147 = pkin(10) ^ 2;
t137 = t152 * pkin(4) + Ifges(5,5);
t1 = [m(2) + m(3) + m(4); Ifges(3,3) + Ifges(4,2) + t139 * t136 + 0.2e1 * t155 + t140 * t135 + (2 * pkin(8) * mrSges(4,3)) + ((pkin(2) ^ 2 + pkin(8) ^ 2) * m(4)); m(4) * pkin(2) + mrSges(3,1); -pkin(8) * m(4) + mrSges(3,2) - mrSges(4,3); t160 * t156 + Ifges(4,1) - Ifges(4,2) - 0.4e1 * t155; t156 * Ifges(5,4) + t160 * t157 + Ifges(4,4); -t143 * Ifges(5,6) + t144 * t137 + Ifges(4,5); t144 * Ifges(5,6) + t143 * t137 + Ifges(4,6); 0.2e1 * pkin(3) * t151 + (t150 * t145) + Ifges(4,3) + Ifges(5,3); mrSges(4,1) + t151; t144 * mrSges(5,2) + t143 * t138 + mrSges(4,2); mrSges(5,3) - t152; m(5) + t145; m(7) * t147 + Ifges(6,1) + t158 - t159; t153 * pkin(5) + Ifges(6,4); Ifges(6,5); Ifges(6,6); Ifges(6,3) + (t147 + t149) * m(7) + t158; m(7) * pkin(5) + mrSges(6,1); mrSges(6,2) - t153; Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
