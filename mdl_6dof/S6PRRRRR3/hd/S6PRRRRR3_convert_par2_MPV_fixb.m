% Return the minimum parameter vector for
% S6PRRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2019-03-09 00:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6PRRRRR3_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR3_convert_par2_MPV_fixb: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR3_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRR3_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRR3_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t160 = -pkin(11) * m(7) - mrSges(7,3);
t140 = (mrSges(6,3) - t160);
t146 = (m(6) + m(7));
t159 = -pkin(10) * t146 - t140;
t149 = (pkin(10) ^ 2);
t152 = (pkin(4) ^ 2);
t157 = (Ifges(5,2) + (t149 + t152) * t146);
t148 = (pkin(11) ^ 2);
t151 = (pkin(5) ^ 2);
t156 = (Ifges(6,2) + (t148 + t151) * m(7));
t144 = (m(5) + t146);
t141 = (m(4) + t144);
t136 = (mrSges(5,3) - t159);
t155 = pkin(9) * t144 + t136;
t154 = 2 * pkin(11) * mrSges(7,3) + 2 * pkin(9) * t136 + 2 * pkin(10) * t140 + Ifges(7,2) + t156 + t157;
t153 = (pkin(3) ^ 2);
t150 = pkin(9) ^ 2;
t1 = [m(2) + m(3) + t141; Ifges(3,3) + Ifges(4,2) + t153 * t144 + 2 * pkin(8) * mrSges(4,3) + (pkin(2) ^ 2 + pkin(8) ^ 2) * t141; pkin(2) * t141 + mrSges(3,1); -pkin(8) * t141 + mrSges(3,2) - mrSges(4,3); Ifges(4,1) - Ifges(4,2) + (t150 - t153) * t144 + t154; pkin(3) * t155 + Ifges(4,4); Ifges(4,5); Ifges(4,6); Ifges(4,3) + (t150 + t153) * t144 + t154; pkin(3) * t144 + mrSges(4,1); mrSges(4,2) - t155; t146 * t149 + Ifges(5,1) - t157; Ifges(5,4); pkin(4) * t159 + Ifges(5,5); Ifges(5,6); t146 * t152 + Ifges(5,3); pkin(4) * t146 + mrSges(5,1); mrSges(5,2); m(7) * t148 + Ifges(6,1) - t156; Ifges(6,4); pkin(5) * t160 + Ifges(6,5); Ifges(6,6); m(7) * t151 + Ifges(6,3); m(7) * pkin(5) + mrSges(6,1); mrSges(6,2); Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
