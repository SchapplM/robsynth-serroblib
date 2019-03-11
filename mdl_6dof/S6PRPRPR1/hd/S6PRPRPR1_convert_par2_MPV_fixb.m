% Return the minimum parameter vector for
% S6PRPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
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
% MPV [21x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6PRPRPR1_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR1_convert_par2_MPV_fixb: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR1_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR1_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRPR1_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t152 = (pkin(9) ^ 2);
t160 = 2 * pkin(9) * mrSges(7,3) + Ifges(7,2);
t138 = m(7) * t152 + Ifges(6,1) + t160;
t153 = (pkin(5) ^ 2);
t142 = m(7) * t153 + Ifges(6,2);
t161 = t138 - t142;
t148 = sin(pkin(12));
t150 = cos(pkin(12));
t159 = t148 * t150;
t145 = t148 ^ 2;
t146 = t150 ^ 2;
t158 = t146 - t145;
t157 = m(7) * pkin(9) + mrSges(7,3);
t139 = pkin(5) * t157 + Ifges(6,4);
t156 = t139 * t159;
t140 = mrSges(6,2) - t157;
t143 = m(7) * pkin(5) + mrSges(6,1);
t155 = -t140 * t148 + t143 * t150;
t141 = -pkin(8) * m(5) + mrSges(4,2) - mrSges(5,3);
t144 = m(5) * pkin(3) + mrSges(4,1);
t149 = sin(pkin(11));
t151 = cos(pkin(11));
t154 = -t141 * t149 + t144 * t151;
t1 = [m(2) + m(3); Ifges(3,3) + Ifges(4,3) + Ifges(5,2) + t145 * t138 + 0.2e1 * t156 + t146 * t142 + (2 * pkin(8) * mrSges(5,3)) + ((pkin(3) ^ 2 + pkin(8) ^ 2) * m(5)) + 0.2e1 * t154 * pkin(2); mrSges(3,1) + t154; t141 * t151 + t144 * t149 + mrSges(3,2); m(4) + m(5); t158 * t161 + Ifges(5,1) - Ifges(5,2) - 0.4e1 * t156; t139 * t158 + t159 * t161 + Ifges(5,4); Ifges(6,5) * t150 - Ifges(6,6) * t148 + Ifges(5,5); Ifges(6,5) * t148 + Ifges(6,6) * t150 + Ifges(5,6); Ifges(5,3) + Ifges(6,3) + ((t152 + t153) * m(7)) + 0.2e1 * t155 * pkin(4) + t160; mrSges(5,1) + t155; t140 * t150 + t143 * t148 + mrSges(5,2); mrSges(6,3); m(6) + m(7); Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
