% Return the minimum parameter vector for
% S6PPRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d5,d6,theta1,theta2,theta4]';
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
% MPV [20x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6PPRPRR1_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRPRR1_convert_par2_MPV_fixb: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRPRR1_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRPRR1_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRPRR1_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t152 = (m(6) + m(7));
t157 = 2 * pkin(10) * mrSges(7,3) + Ifges(7,2);
t156 = pkin(10) * m(7) + mrSges(7,3);
t147 = -pkin(9) * t152 + mrSges(5,2) - mrSges(6,3);
t148 = pkin(4) * t152 + mrSges(5,1);
t150 = sin(pkin(13));
t151 = cos(pkin(13));
t155 = -t150 * t147 + t151 * t148;
t154 = (pkin(5) ^ 2);
t153 = pkin(10) ^ 2;
t1 = [m(2); m(3) + m(4); Ifges(4,3) + Ifges(5,3) + (t154 * m(7)) + Ifges(6,2) + (2 * pkin(9) * mrSges(6,3)) + ((pkin(4) ^ 2 + pkin(9) ^ 2) * t152) + 0.2e1 * t155 * pkin(3); mrSges(4,1) + t155; t147 * t151 + t148 * t150 + mrSges(4,2); m(5) + t152; Ifges(6,1) - Ifges(6,2) + ((t153 - t154) * m(7)) + t157; t156 * pkin(5) + Ifges(6,4); Ifges(6,5); Ifges(6,6); Ifges(6,3) + (t153 + t154) * m(7) + t157; m(7) * pkin(5) + mrSges(6,1); mrSges(6,2) - t156; Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
