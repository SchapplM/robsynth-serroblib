% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PPRRRR2
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 09:26
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PPRRRR2_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR2_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR2_jacobiR_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 09:26:44
% EndTime: 2019-02-22 09:26:44
% DurationCPUTime: 0.14s
% Computational Cost: add. (279->44), mult. (694->91), div. (0->0), fcn. (958->14), ass. (0->53)
t232 = cos(qJ(3));
t206 = qJ(5) + qJ(6);
t204 = sin(t206);
t217 = cos(qJ(4));
t231 = t204 * t217;
t205 = cos(t206);
t230 = t205 * t217;
t208 = sin(pkin(12));
t214 = cos(pkin(6));
t229 = t208 * t214;
t209 = sin(pkin(7));
t228 = t209 * t214;
t210 = sin(pkin(6));
t227 = t210 * t209;
t213 = cos(pkin(7));
t226 = t210 * t213;
t211 = cos(pkin(13));
t225 = t211 * t213;
t212 = cos(pkin(12));
t224 = t212 * t214;
t223 = t210 * t232;
t222 = t209 * t223;
t207 = sin(pkin(13));
t221 = -t208 * t207 + t211 * t224;
t220 = t212 * t207 + t211 * t229;
t219 = t221 * t213;
t218 = t220 * t213;
t216 = sin(qJ(3));
t215 = sin(qJ(4));
t200 = -t207 * t229 + t212 * t211;
t199 = t207 * t224 + t208 * t211;
t198 = -t211 * t227 + t214 * t213;
t195 = t208 * t226 + t209 * t220;
t194 = -t209 * t221 - t212 * t226;
t193 = t216 * t228 + (t232 * t207 + t216 * t225) * t210;
t192 = t210 * t207 * t216 - t223 * t225 - t232 * t228;
t191 = t193 * t217 + t198 * t215;
t190 = -t193 * t215 + t198 * t217;
t189 = t200 * t232 + (t208 * t227 - t218) * t216;
t188 = t200 * t216 - t208 * t222 + t232 * t218;
t187 = t199 * t232 + (-t212 * t227 + t219) * t216;
t186 = t199 * t216 + t212 * t222 - t232 * t219;
t185 = t189 * t217 + t195 * t215;
t184 = -t189 * t215 + t195 * t217;
t183 = t187 * t217 + t194 * t215;
t182 = -t187 * t215 + t194 * t217;
t181 = -t191 * t205 - t192 * t204;
t180 = -t191 * t204 + t192 * t205;
t179 = -t185 * t205 - t188 * t204;
t178 = -t185 * t204 + t188 * t205;
t177 = -t183 * t205 - t186 * t204;
t176 = -t183 * t204 + t186 * t205;
t1 = [0, 0, -t188 * t230 + t189 * t204, t184 * t205, t178, t178; 0, 0, -t186 * t230 + t187 * t204, t182 * t205, t176, t176; 0, 0, -t192 * t230 + t193 * t204, t190 * t205, t180, t180; 0, 0, t188 * t231 + t189 * t205, -t184 * t204, t179, t179; 0, 0, t186 * t231 + t187 * t205, -t182 * t204, t177, t177; 0, 0, t192 * t231 + t193 * t205, -t190 * t204, t181, t181; 0, 0, -t188 * t215, t185, 0, 0; 0, 0, -t186 * t215, t183, 0, 0; 0, 0, -t192 * t215, t191, 0, 0;];
JR_rot  = t1;
