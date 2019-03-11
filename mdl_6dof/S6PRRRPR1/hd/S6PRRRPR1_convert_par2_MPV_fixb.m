% Return the minimum parameter vector for
% S6PRRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
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
% Datum: 2019-03-08 23:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6PRRRPR1_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR1_convert_par2_MPV_fixb: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR1_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR1_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR1_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t146 = (pkin(10) ^ 2);
t157 = 2 * pkin(10) * mrSges(7,3) + Ifges(7,2);
t134 = m(7) * t146 + Ifges(6,1) + t157;
t148 = (pkin(5) ^ 2);
t137 = t148 * m(7) + Ifges(6,2);
t158 = t134 - t137;
t145 = (m(4) + m(5));
t143 = sin(pkin(12));
t144 = cos(pkin(12));
t156 = t143 * t144;
t140 = t143 ^ 2;
t141 = t144 ^ 2;
t155 = t141 - t140;
t154 = -pkin(9) * m(5) - mrSges(5,3);
t153 = pkin(10) * m(7) + mrSges(7,3);
t135 = t153 * pkin(5) + Ifges(6,4);
t152 = t135 * t156;
t151 = (mrSges(4,3) - t154);
t136 = mrSges(6,2) - t153;
t138 = m(7) * pkin(5) + mrSges(6,1);
t150 = -t143 * t136 + t144 * t138;
t149 = pkin(3) ^ 2;
t147 = pkin(9) ^ 2;
t139 = t147 + t149;
t1 = [m(2) + m(3) + t145; Ifges(3,3) + Ifges(4,2) + Ifges(5,2) + t140 * t134 + 0.2e1 * t152 + t141 * t137 + (2 * pkin(9) * mrSges(5,3)) + t139 * m(5) + (2 * pkin(8) * t151) + ((pkin(2) ^ 2 + pkin(8) ^ 2) * t145); pkin(2) * t145 + mrSges(3,1); -pkin(8) * t145 + mrSges(3,2) - t151; Ifges(4,1) - Ifges(4,2) + (-t139 + t147) * m(5); Ifges(4,4); t154 * pkin(3) + Ifges(4,5); Ifges(4,6); t149 * m(5) + Ifges(4,3); m(5) * pkin(3) + mrSges(4,1); mrSges(4,2); t158 * t155 + Ifges(5,1) - Ifges(5,2) - 0.4e1 * t152; t155 * t135 + t158 * t156 + Ifges(5,4); t144 * Ifges(6,5) - t143 * Ifges(6,6) + Ifges(5,5); t143 * Ifges(6,5) + t144 * Ifges(6,6) + Ifges(5,6); Ifges(5,3) + Ifges(6,3) + ((t146 + t148) * m(7)) + 0.2e1 * t150 * pkin(4) + t157; mrSges(5,1) + t150; t144 * t136 + t143 * t138 + mrSges(5,2); mrSges(6,3); m(6) + m(7); Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
