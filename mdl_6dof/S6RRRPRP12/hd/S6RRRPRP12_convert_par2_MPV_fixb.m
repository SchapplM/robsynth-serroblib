% Return the minimum parameter vector for
% S6RRRPRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5]';
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
% Datum: 2019-03-09 18:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RRRPRP12_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP12_convert_par2_MPV_fixb: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP12_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP12_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP12_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t162 = (pkin(10) * mrSges(6,3));
t151 = m(3) + m(4);
t161 = (t151 * pkin(8));
t160 = (Ifges(6,2) + Ifges(7,3));
t159 = pkin(9) * m(4) + mrSges(4,3);
t158 = -pkin(10) * m(6) - mrSges(6,3);
t152 = (pkin(10) ^ 2);
t154 = (pkin(4) ^ 2);
t157 = (Ifges(4,2) + Ifges(5,3) + (t152 + t154) * m(6) + t160);
t148 = 2 * t162;
t156 = 2 * pkin(9) * mrSges(4,3) + t148 + t157;
t155 = (pkin(2) ^ 2);
t153 = pkin(9) ^ 2;
t150 = sin(pkin(6));
t1 = [pkin(1) ^ 2 * t151 + Ifges(2,3) + (t155 * m(4) + Ifges(3,2) + (2 * mrSges(3,3) + t161) * pkin(8)) * t150 ^ 2; pkin(1) * t151 + mrSges(2,1); mrSges(2,2) + (-mrSges(3,3) - t161) * t150; Ifges(3,1) - Ifges(3,2) + (t153 - t155) * m(4) + t156; t159 * pkin(2) + Ifges(3,4); Ifges(3,5); Ifges(3,6); Ifges(3,3) + (t153 + t155) * m(4) + t156; m(4) * pkin(2) + mrSges(3,1); mrSges(3,2) - t159; t154 * m(6) + Ifges(4,1) + Ifges(5,2) - t157 - 2 * t162; Ifges(4,4) + Ifges(5,6); t158 * pkin(4) - Ifges(5,4) + Ifges(4,5); Ifges(4,6) - Ifges(5,5); t152 * m(6) + Ifges(5,1) + Ifges(4,3) + t148 + t160; mrSges(4,1); mrSges(4,2); m(6) * pkin(4) + mrSges(5,1); mrSges(5,2) + t158; mrSges(5,3); m(5) + m(6); Ifges(6,1) + Ifges(7,1) - t160; Ifges(6,4) - Ifges(7,5); Ifges(6,5) + Ifges(7,4); Ifges(6,6) - Ifges(7,6); Ifges(6,3) + Ifges(7,2); mrSges(6,1); mrSges(6,2); mrSges(7,1); mrSges(7,2); mrSges(7,3); m(7);];
MPV  = t1;
