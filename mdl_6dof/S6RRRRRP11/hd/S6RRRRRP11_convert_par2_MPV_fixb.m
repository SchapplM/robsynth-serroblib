% Return the minimum parameter vector for
% S6RRRRRP11
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
% MPV [33x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 02:58
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RRRRRP11_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP11_convert_par2_MPV_fixb: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP11_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP11_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP11_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t189 = sin(pkin(7));
t192 = (m(5) + m(6));
t186 = (m(4) + t192);
t202 = pkin(10) * t186 + mrSges(4,3);
t212 = t202 * t189;
t191 = cos(pkin(7));
t211 = t202 * t191;
t210 = (pkin(10) * mrSges(4,3));
t181 = m(3) + t186;
t209 = t181 * pkin(9);
t208 = (-Ifges(6,2) - Ifges(7,2));
t198 = (pkin(3) ^ 2);
t180 = (t198 * t192 + Ifges(4,2));
t197 = (pkin(4) ^ 2);
t207 = (t197 * m(6) + Ifges(5,2));
t206 = 2 * pkin(12) * mrSges(6,3) - t208;
t205 = 2 * pkin(11) * mrSges(5,3) + t207;
t204 = pkin(12) * m(6) + mrSges(6,3);
t203 = t180 + 2 * t210;
t201 = pkin(11) * t192 + mrSges(5,3);
t196 = (pkin(10) ^ 2);
t200 = t196 * t186 + t203;
t199 = pkin(2) ^ 2;
t195 = pkin(11) ^ 2;
t194 = pkin(12) ^ 2;
t190 = sin(pkin(6));
t185 = t191 ^ 2;
t179 = t196 * t185 + t199;
t178 = mrSges(3,3) + t211;
t1 = [pkin(1) ^ 2 * t181 + Ifges(2,3) + (t179 * t186 + Ifges(3,2) + t203 * t185 + (0.2e1 * t178 + t209) * pkin(9)) * t190 ^ 2; pkin(1) * t181 + mrSges(2,1); mrSges(2,2) + (-t178 - t209) * t190; -t185 * t180 + Ifges(3,1) - Ifges(3,2) + (-t179 + t196) * t186 + (-0.2e1 * t185 + 0.2e1) * t210 + t180; pkin(2) * t212 + Ifges(3,4); -pkin(2) * t211 + Ifges(3,5); t200 * t191 * t189 + Ifges(3,6); t200 * t189 ^ 2 + t199 * t186 + Ifges(3,3); pkin(2) * t186 + mrSges(3,1); mrSges(3,2) - t212; t195 * t192 + Ifges(4,1) - t180 + t205; t201 * pkin(3) + Ifges(4,4); Ifges(4,5); Ifges(4,6); Ifges(4,3) + (t195 + t198) * t192 + t205; pkin(3) * t192 + mrSges(4,1); mrSges(4,2) - t201; t194 * m(6) + Ifges(5,1) + t206 - t207; t204 * pkin(4) + Ifges(5,4); Ifges(5,5); Ifges(5,6); Ifges(5,3) + (t194 + t197) * m(6) + t206; m(6) * pkin(4) + mrSges(5,1); mrSges(5,2) - t204; Ifges(6,1) + Ifges(7,1) + t208; Ifges(6,4) + Ifges(7,4); Ifges(6,5) + Ifges(7,5); Ifges(6,6) + Ifges(7,6); 2 * pkin(5) * mrSges(7,1) + Ifges(6,3) + Ifges(7,3); mrSges(6,1) + mrSges(7,1); mrSges(6,2) + mrSges(7,2); mrSges(7,3); m(7);];
MPV  = t1;
