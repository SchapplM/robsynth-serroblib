% Return the minimum parameter vector for
% S6PRPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
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
% MPV [26x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6PRPRRR1_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR1_convert_par2_MPV_fixb: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR1_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR1_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRR1_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t143 = m(6) + m(7);
t152 = 2 * pkin(10) * mrSges(7,3) + Ifges(7,2);
t139 = (m(5) + t143);
t151 = pkin(10) * m(7) + mrSges(7,3);
t150 = -pkin(9) * t143 - mrSges(6,3);
t149 = (mrSges(5,3) - t150);
t136 = -pkin(8) * t139 + mrSges(4,2) - t149;
t137 = pkin(3) * t139 + mrSges(4,1);
t141 = sin(pkin(12));
t142 = cos(pkin(12));
t148 = -t141 * t136 + t142 * t137;
t147 = pkin(4) ^ 2;
t146 = (pkin(5) ^ 2);
t145 = pkin(9) ^ 2;
t144 = pkin(10) ^ 2;
t138 = t145 + t147;
t1 = [m(2) + m(3); Ifges(3,3) + Ifges(4,3) + Ifges(5,2) + (t146 * m(7)) + Ifges(6,2) + (2 * pkin(9) * mrSges(6,3)) + t138 * t143 + (2 * pkin(8) * t149) + ((pkin(3) ^ 2 + pkin(8) ^ 2) * t139) + 0.2e1 * t148 * pkin(2); mrSges(3,1) + t148; t142 * t136 + t141 * t137 + mrSges(3,2); m(4) + t139; Ifges(5,1) - Ifges(5,2) + (-t138 + t145) * t143; Ifges(5,4); t150 * pkin(4) + Ifges(5,5); Ifges(5,6); t147 * t143 + Ifges(5,3); pkin(4) * t143 + mrSges(5,1); mrSges(5,2); Ifges(6,1) - Ifges(6,2) + ((t144 - t146) * m(7)) + t152; t151 * pkin(5) + Ifges(6,4); Ifges(6,5); Ifges(6,6); Ifges(6,3) + (t144 + t146) * m(7) + t152; m(7) * pkin(5) + mrSges(6,1); mrSges(6,2) - t151; Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
