% Return the minimum parameter vector for
% S6RPRPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
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
% Datum: 2019-03-09 04:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RPRPRR11_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR11_convert_par2_MPV_fixb: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR11_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR11_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR11_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t232 = (m(4) * pkin(9));
t204 = sin(pkin(7));
t222 = mrSges(4,3) + t232;
t231 = t222 * t204;
t208 = cos(pkin(7));
t230 = t222 * t208;
t210 = (m(6) + m(7));
t216 = (pkin(4) ^ 2);
t193 = (t216 * t210 + Ifges(4,2) + Ifges(5,3));
t218 = t193 + (2 * mrSges(4,3) + t232) * pkin(9);
t229 = pkin(2) ^ 2 * m(4);
t215 = (pkin(5) ^ 2);
t228 = (t215 * m(7) + Ifges(6,2));
t227 = 2 * pkin(11) * mrSges(7,3) + Ifges(7,2);
t202 = sin(pkin(13));
t206 = cos(pkin(13));
t226 = t202 * t206;
t224 = Ifges(5,4) * t226;
t223 = 2 * pkin(10) * mrSges(6,3) + t228;
t221 = pkin(11) * m(7) + mrSges(7,3);
t220 = -pkin(10) * t210 - mrSges(6,3);
t213 = pkin(10) ^ 2;
t212 = pkin(11) ^ 2;
t209 = cos(pkin(6));
t207 = cos(pkin(12));
t205 = sin(pkin(6));
t203 = sin(pkin(12));
t198 = t206 ^ 2;
t195 = t202 ^ 2;
t192 = t220 * pkin(4) + Ifges(5,5);
t191 = t213 * t210 + Ifges(5,1) + t223;
t190 = Ifges(5,2) + (t213 + t216) * t210 + t223;
t1 = [Ifges(2,3) + t209 ^ 2 * (t218 * t204 ^ 2 + Ifges(3,3) + t229) + (0.2e1 * (t203 * (-pkin(2) * t230 + Ifges(3,5)) + t207 * (t218 * t208 * t204 + Ifges(3,6))) * t209 + (t207 ^ 2 * (t218 * t208 ^ 2 + Ifges(3,2) + t229) + (0.2e1 * t207 * (pkin(2) * t231 + Ifges(3,4)) + (Ifges(3,1) + t218) * t203) * t203) * t205) * t205; mrSges(2,1); mrSges(2,2); m(4) * pkin(2) + mrSges(3,1); mrSges(3,2) - t231; mrSges(3,3) + t230; m(3) + m(4); t195 * t190 + t198 * t191 + Ifges(4,1) - t193 - 0.2e1 * t224; t202 * Ifges(5,6) - t206 * t192 + Ifges(4,4); Ifges(4,5) + (t198 - t195) * Ifges(5,4) + (-t190 + t191) * t226; -t206 * Ifges(5,6) - t202 * t192 + Ifges(4,6); t198 * t190 + t195 * t191 + Ifges(4,3) + 0.2e1 * t224; mrSges(4,1); mrSges(4,2); pkin(4) * t210 + mrSges(5,1); mrSges(5,2); mrSges(5,3) - t220; m(5) + t210; m(7) * t212 + Ifges(6,1) + t227 - t228; t221 * pkin(5) + Ifges(6,4); Ifges(6,5); Ifges(6,6); Ifges(6,3) + (t212 + t215) * m(7) + t227; m(7) * pkin(5) + mrSges(6,1); mrSges(6,2) - t221; Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
