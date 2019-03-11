% Return the minimum parameter vector for
% S6RPRRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2]';
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
% Datum: 2019-03-09 05:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RPRRPR12_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRPR12_convert_par2_MPV_fixb: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR12_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR12_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR12_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t197 = m(4) + m(5);
t217 = (pkin(9) * t197);
t192 = sin(pkin(7));
t209 = mrSges(4,3) + t217;
t216 = t209 * t192;
t195 = cos(pkin(7));
t215 = t209 * t195;
t203 = (pkin(3) ^ 2);
t184 = (t203 * m(5) + Ifges(4,2));
t205 = t184 + (2 * mrSges(4,3) + t217) * pkin(9);
t214 = (pkin(11) * mrSges(7,3));
t213 = pkin(2) ^ 2 * t197;
t211 = pkin(10) * m(5) + mrSges(5,3);
t210 = -pkin(11) * m(7) - mrSges(7,3);
t199 = (pkin(11) ^ 2);
t202 = (pkin(5) ^ 2);
t208 = (Ifges(5,2) + Ifges(7,2) + Ifges(6,3) + (t199 + t202) * m(7));
t189 = 2 * t214;
t206 = 2 * pkin(10) * mrSges(5,3) + t189 + t208;
t200 = pkin(10) ^ 2;
t196 = cos(pkin(6));
t194 = cos(pkin(12));
t193 = sin(pkin(6));
t191 = sin(pkin(12));
t1 = [Ifges(2,3) + t196 ^ 2 * (t205 * t192 ^ 2 + Ifges(3,3) + t213) + (0.2e1 * (t191 * (-pkin(2) * t215 + Ifges(3,5)) + t194 * (t205 * t195 * t192 + Ifges(3,6))) * t196 + (t194 ^ 2 * (t205 * t195 ^ 2 + Ifges(3,2) + t213) + (0.2e1 * t194 * (pkin(2) * t216 + Ifges(3,4)) + (Ifges(3,1) + t205) * t191) * t191) * t193) * t193; mrSges(2,1); mrSges(2,2); pkin(2) * t197 + mrSges(3,1); mrSges(3,2) - t216; mrSges(3,3) + t215; m(3) + t197; t200 * m(5) + Ifges(4,1) - t184 + t206; t211 * pkin(3) + Ifges(4,4); Ifges(4,5); Ifges(4,6); Ifges(4,3) + (t200 + t203) * m(5) + t206; m(5) * pkin(3) + mrSges(4,1); mrSges(4,2) - t211; t202 * m(7) + Ifges(5,1) + Ifges(6,2) - t208 - 2 * t214; Ifges(5,4) + Ifges(6,6); t210 * pkin(5) - Ifges(6,4) + Ifges(5,5); Ifges(5,6) - Ifges(6,5); m(7) * t199 + Ifges(6,1) + Ifges(7,2) + Ifges(5,3) + t189; mrSges(5,1); mrSges(5,2); m(7) * pkin(5) + mrSges(6,1); mrSges(6,2) + t210; mrSges(6,3); m(6) + m(7); Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
