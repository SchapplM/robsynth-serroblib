% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRR11
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 12:09
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRPRR11_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR11_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR11_jacobiR_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 12:09:31
% EndTime: 2019-02-22 12:09:31
% DurationCPUTime: 0.20s
% Computational Cost: add. (189->42), mult. (534->79), div. (0->0), fcn. (756->12), ass. (0->50)
t231 = sin(qJ(2));
t232 = sin(qJ(1));
t236 = cos(qJ(2));
t237 = cos(qJ(1));
t249 = cos(pkin(6));
t241 = t237 * t249;
t220 = t231 * t241 + t232 * t236;
t230 = sin(qJ(3));
t235 = cos(qJ(3));
t227 = sin(pkin(6));
t243 = t227 * t237;
t210 = t220 * t230 + t235 * t243;
t213 = -t220 * t235 + t230 * t243;
t229 = sin(qJ(5));
t234 = cos(qJ(5));
t199 = t210 * t229 - t213 * t234;
t219 = t232 * t231 - t236 * t241;
t228 = sin(qJ(6));
t233 = cos(qJ(6));
t257 = t199 * t228 + t219 * t233;
t256 = -t199 * t233 + t219 * t228;
t198 = t210 * t234 + t213 * t229;
t255 = t198 * t228;
t254 = t198 * t233;
t246 = t227 * t231;
t217 = t230 * t246 - t249 * t235;
t245 = t227 * t235;
t218 = t249 * t230 + t231 * t245;
t207 = t217 * t234 - t218 * t229;
t253 = t207 * t228;
t252 = t207 * t233;
t242 = t232 * t249;
t222 = -t231 * t242 + t237 * t236;
t214 = t222 * t230 - t232 * t245;
t215 = t232 * t227 * t230 + t222 * t235;
t240 = -t214 * t234 + t215 * t229;
t251 = t240 * t228;
t250 = t240 * t233;
t244 = t227 * t236;
t203 = t214 * t229 + t215 * t234;
t208 = t217 * t229 + t218 * t234;
t239 = t229 * t235 - t230 * t234;
t238 = t229 * t230 + t234 * t235;
t221 = -t237 * t231 - t236 * t242;
t216 = t238 * t244;
t205 = t238 * t221;
t204 = t238 * t219;
t196 = t203 * t233 + t221 * t228;
t195 = -t203 * t228 + t221 * t233;
t1 = [t256, t205 * t233 - t222 * t228, t250, 0, -t250, t195; t196, -t204 * t233 - t220 * t228, -t254, 0, t254, -t257; 0, t216 * t233 - t228 * t246, -t252, 0, t252, -t208 * t228 + t233 * t244; t257, -t205 * t228 - t222 * t233, -t251, 0, t251, -t196; t195, t204 * t228 - t220 * t233, t255, 0, -t255, t256; 0, -t216 * t228 - t233 * t246, t253, 0, -t253, -t208 * t233 - t228 * t244; t198, t239 * t221, -t203, 0, t203, 0; t240, -t239 * t219, -t199, 0, t199, 0; 0, t239 * t244, -t208, 0, t208, 0;];
JR_rot  = t1;
