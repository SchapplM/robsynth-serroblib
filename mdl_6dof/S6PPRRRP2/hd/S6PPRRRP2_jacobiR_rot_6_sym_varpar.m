% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PPRRRP2
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 09:25
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PPRRRP2_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP2_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP2_jacobiR_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 09:25:25
% EndTime: 2019-02-22 09:25:25
% DurationCPUTime: 0.15s
% Computational Cost: add. (172->40), mult. (516->91), div. (0->0), fcn. (712->14), ass. (0->46)
t230 = cos(qJ(3));
t204 = sin(pkin(11));
t210 = cos(pkin(6));
t229 = t204 * t210;
t205 = sin(pkin(7));
t228 = t205 * t210;
t206 = sin(pkin(6));
t227 = t206 * t205;
t209 = cos(pkin(7));
t226 = t206 * t209;
t207 = cos(pkin(12));
t225 = t207 * t209;
t208 = cos(pkin(11));
t224 = t208 * t210;
t211 = sin(qJ(5));
t215 = cos(qJ(4));
t223 = t211 * t215;
t214 = cos(qJ(5));
t222 = t214 * t215;
t221 = t206 * t230;
t220 = t205 * t221;
t203 = sin(pkin(12));
t219 = -t204 * t203 + t207 * t224;
t218 = t208 * t203 + t207 * t229;
t217 = t219 * t209;
t216 = t218 * t209;
t213 = sin(qJ(3));
t212 = sin(qJ(4));
t199 = -t203 * t229 + t208 * t207;
t198 = t203 * t224 + t204 * t207;
t197 = -t207 * t227 + t210 * t209;
t194 = t204 * t226 + t218 * t205;
t193 = -t219 * t205 - t208 * t226;
t192 = t213 * t228 + (t230 * t203 + t213 * t225) * t206;
t191 = t206 * t203 * t213 - t221 * t225 - t230 * t228;
t190 = t192 * t215 + t197 * t212;
t189 = -t192 * t212 + t197 * t215;
t188 = t199 * t230 + (t204 * t227 - t216) * t213;
t187 = t199 * t213 - t204 * t220 + t230 * t216;
t186 = t198 * t230 + (-t208 * t227 + t217) * t213;
t185 = t198 * t213 + t208 * t220 - t230 * t217;
t184 = t188 * t215 + t194 * t212;
t183 = -t188 * t212 + t194 * t215;
t182 = t186 * t215 + t193 * t212;
t181 = -t186 * t212 + t193 * t215;
t1 = [0, 0, -t187 * t222 + t188 * t211, t183 * t214, -t184 * t211 + t187 * t214, 0; 0, 0, -t185 * t222 + t186 * t211, t181 * t214, -t182 * t211 + t185 * t214, 0; 0, 0, -t191 * t222 + t192 * t211, t189 * t214, -t190 * t211 + t191 * t214, 0; 0, 0, -t187 * t212, t184, 0, 0; 0, 0, -t185 * t212, t182, 0, 0; 0, 0, -t191 * t212, t190, 0, 0; 0, 0, -t187 * t223 - t188 * t214, t183 * t211, t184 * t214 + t187 * t211, 0; 0, 0, -t185 * t223 - t186 * t214, t181 * t211, t182 * t214 + t185 * t211, 0; 0, 0, -t191 * t223 - t192 * t214, t189 * t211, t190 * t214 + t191 * t211, 0;];
JR_rot  = t1;
