% Return the minimum parameter vector for
% S6RRRRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5]';
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
% MPV [35x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 03:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RRRRRP12_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP12_convert_par2_MPV_fixb: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP12_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP12_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP12_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t197 = sin(pkin(7));
t200 = (m(5) + m(6));
t194 = (m(4) + t200);
t210 = pkin(10) * t194 + mrSges(4,3);
t220 = t210 * t197;
t199 = cos(pkin(7));
t219 = t210 * t199;
t218 = (pkin(10) * mrSges(4,3));
t189 = m(3) + t194;
t217 = t189 * pkin(9);
t216 = (-Ifges(6,2) - Ifges(7,3));
t206 = (pkin(3) ^ 2);
t188 = (t206 * t200 + Ifges(4,2));
t205 = (pkin(4) ^ 2);
t215 = (t205 * m(6) + Ifges(5,2));
t214 = 2 * pkin(12) * mrSges(6,3) - t216;
t213 = 2 * pkin(11) * mrSges(5,3) + t215;
t212 = pkin(12) * m(6) + mrSges(6,3);
t211 = t188 + 2 * t218;
t209 = pkin(11) * t200 + mrSges(5,3);
t204 = (pkin(10) ^ 2);
t208 = t204 * t194 + t211;
t207 = pkin(2) ^ 2;
t203 = pkin(11) ^ 2;
t202 = pkin(12) ^ 2;
t198 = sin(pkin(6));
t193 = t199 ^ 2;
t187 = t204 * t193 + t207;
t186 = mrSges(3,3) + t219;
t1 = [pkin(1) ^ 2 * t189 + Ifges(2,3) + (t187 * t194 + Ifges(3,2) + t211 * t193 + (0.2e1 * t186 + t217) * pkin(9)) * t198 ^ 2; pkin(1) * t189 + mrSges(2,1); mrSges(2,2) + (-t186 - t217) * t198; -t193 * t188 + Ifges(3,1) - Ifges(3,2) + (-t187 + t204) * t194 + (-0.2e1 * t193 + 0.2e1) * t218 + t188; pkin(2) * t220 + Ifges(3,4); -pkin(2) * t219 + Ifges(3,5); t208 * t199 * t197 + Ifges(3,6); t208 * t197 ^ 2 + t207 * t194 + Ifges(3,3); pkin(2) * t194 + mrSges(3,1); mrSges(3,2) - t220; t203 * t200 + Ifges(4,1) - t188 + t213; t209 * pkin(3) + Ifges(4,4); Ifges(4,5); Ifges(4,6); Ifges(4,3) + (t203 + t206) * t200 + t213; pkin(3) * t200 + mrSges(4,1); mrSges(4,2) - t209; t202 * m(6) + Ifges(5,1) + t214 - t215; t212 * pkin(4) + Ifges(5,4); Ifges(5,5); Ifges(5,6); Ifges(5,3) + (t202 + t205) * m(6) + t214; m(6) * pkin(4) + mrSges(5,1); mrSges(5,2) - t212; Ifges(6,1) + Ifges(7,1) + t216; Ifges(6,4) - Ifges(7,5); Ifges(6,5) + Ifges(7,4); Ifges(6,6) - Ifges(7,6); Ifges(6,3) + Ifges(7,2); mrSges(6,1); mrSges(6,2); mrSges(7,1); mrSges(7,2); mrSges(7,3); m(7);];
MPV  = t1;
