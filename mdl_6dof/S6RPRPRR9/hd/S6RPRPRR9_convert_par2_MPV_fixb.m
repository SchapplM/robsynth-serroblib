% Return the minimum parameter vector for
% S6RPRPRR9
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
% MPV [30x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RPRPRR9_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR9_convert_par2_MPV_fixb: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR9_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR9_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR9_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t233 = (m(4) * pkin(9));
t203 = sin(pkin(7));
t223 = mrSges(4,3) + t233;
t232 = t223 * t203;
t207 = cos(pkin(7));
t231 = t223 * t207;
t209 = (m(6) + m(7));
t220 = pkin(10) * t209 + mrSges(6,3);
t190 = t220 * pkin(4) + Ifges(5,4);
t201 = sin(pkin(13));
t205 = cos(pkin(13));
t227 = t201 * t205;
t221 = t190 * t227;
t212 = (pkin(10) ^ 2);
t214 = (pkin(5) ^ 2);
t229 = (t214 * m(7) + Ifges(6,2));
t224 = 2 * pkin(10) * mrSges(6,3) + t229;
t188 = t212 * t209 + Ifges(5,1) + t224;
t215 = (pkin(4) ^ 2);
t193 = t215 * t209 + Ifges(5,2);
t194 = t201 ^ 2;
t197 = t205 ^ 2;
t225 = t194 * t188 + t197 * t193 + Ifges(4,2);
t217 = 0.2e1 * t221 + t225 + ((2 * mrSges(4,3) + t233) * pkin(9));
t230 = pkin(2) ^ 2 * m(4);
t228 = 2 * pkin(11) * mrSges(7,3) + Ifges(7,2);
t222 = pkin(11) * m(7) + mrSges(7,3);
t191 = mrSges(5,2) - t220;
t192 = pkin(4) * t209 + mrSges(5,1);
t218 = -t201 * t191 + t205 * t192;
t211 = pkin(11) ^ 2;
t208 = cos(pkin(6));
t206 = cos(pkin(12));
t204 = sin(pkin(6));
t202 = sin(pkin(12));
t1 = [Ifges(2,3) + t208 ^ 2 * (t217 * t203 ^ 2 + Ifges(3,3) + t230) + (0.2e1 * (t202 * (-pkin(2) * t231 + Ifges(3,5)) + t206 * (t217 * t207 * t203 + Ifges(3,6))) * t208 + (t206 ^ 2 * (t217 * t207 ^ 2 + Ifges(3,2) + t230) + (0.2e1 * t206 * (pkin(2) * t232 + Ifges(3,4)) + (Ifges(3,1) + t217) * t202) * t202) * t204) * t204; mrSges(2,1); mrSges(2,2); m(4) * pkin(2) + mrSges(3,1); mrSges(3,2) - t232; mrSges(3,3) + t231; m(3) + m(4); t197 * t188 + t194 * t193 + Ifges(4,1) - 0.4e1 * t221 - t225; Ifges(4,4) + (t197 - t194) * t190 + (t188 - t193) * t227; t205 * Ifges(5,5) - t201 * Ifges(5,6) + Ifges(4,5); t201 * Ifges(5,5) + t205 * Ifges(5,6) + Ifges(4,6); Ifges(4,3) + Ifges(5,3) + ((t212 + t215) * t209) + 0.2e1 * t218 * pkin(3) + t224; mrSges(4,1) + t218; t205 * t191 + t201 * t192 + mrSges(4,2); mrSges(5,3); m(5) + t209; m(7) * t211 + Ifges(6,1) + t228 - t229; t222 * pkin(5) + Ifges(6,4); Ifges(6,5); Ifges(6,6); Ifges(6,3) + (t211 + t214) * m(7) + t228; m(7) * pkin(5) + mrSges(6,1); mrSges(6,2) - t222; Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
