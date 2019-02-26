% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRRP11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,theta2]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:13
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRRRP11_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP11_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP11_jacobiR_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:13:42
% EndTime: 2019-02-26 21:13:42
% DurationCPUTime: 0.26s
% Computational Cost: add. (237->47), mult. (692->91), div. (0->0), fcn. (956->14), ass. (0->54)
t218 = cos(pkin(6));
t213 = sin(pkin(12));
t226 = cos(qJ(1));
t234 = t226 * t213;
t216 = cos(pkin(12));
t222 = sin(qJ(1));
t236 = t222 * t216;
t208 = t218 * t234 + t236;
t221 = sin(qJ(3));
t225 = cos(qJ(3));
t233 = t226 * t216;
t237 = t222 * t213;
t207 = -t218 * t233 + t237;
t217 = cos(pkin(7));
t214 = sin(pkin(7));
t215 = sin(pkin(6));
t240 = t215 * t226;
t232 = t214 * t240;
t230 = t207 * t217 + t232;
t196 = -t208 * t225 + t230 * t221;
t202 = -t207 * t214 + t217 * t240;
t220 = sin(qJ(4));
t224 = cos(qJ(4));
t188 = t196 * t224 + t202 * t220;
t219 = sin(qJ(5));
t250 = t188 * t219;
t223 = cos(qJ(5));
t249 = t188 * t223;
t186 = t196 * t220 - t202 * t224;
t229 = t218 * t236 + t234;
t241 = t215 * t222;
t246 = -t214 * t241 + t229 * t217;
t245 = t208 * t221;
t243 = t214 * t218;
t242 = t215 * t216;
t239 = t217 * t225;
t238 = t219 * t224;
t235 = t223 * t224;
t227 = t229 * t214 + t217 * t241;
t209 = -t218 * t237 + t233;
t206 = -t214 * t242 + t218 * t217;
t201 = t221 * t243 + (t216 * t217 * t221 + t213 * t225) * t215;
t200 = t215 * t213 * t221 - t225 * t243 - t239 * t242;
t198 = t209 * t225 - t246 * t221;
t197 = t209 * t221 + t246 * t225;
t195 = -t230 * t225 - t245;
t193 = t207 * t239 + t225 * t232 + t245;
t192 = t201 * t224 + t206 * t220;
t191 = -t201 * t220 + t206 * t224;
t190 = t198 * t224 + t227 * t220;
t189 = t198 * t220 - t227 * t224;
t185 = t190 * t223 + t197 * t219;
t184 = -t190 * t219 + t197 * t223;
t1 = [t195 * t219 + t249, 0, -t197 * t235 + t198 * t219, -t189 * t223, t184, 0; t185, 0, -t193 * t235 - t196 * t219, t186 * t223, t193 * t223 + t250, 0; 0, 0, -t200 * t235 + t201 * t219, t191 * t223, -t192 * t219 + t200 * t223, 0; t195 * t223 - t250, 0, t197 * t238 + t198 * t223, t189 * t219, -t185, 0; t184, 0, t193 * t238 - t196 * t223, -t186 * t219, -t193 * t219 + t249, 0; 0, 0, t200 * t238 + t201 * t223, -t191 * t219, -t192 * t223 - t200 * t219, 0; t186, 0, -t197 * t220, t190, 0, 0; t189, 0, -t193 * t220, -t188, 0, 0; 0, 0, -t200 * t220, t192, 0, 0;];
JR_rot  = t1;
