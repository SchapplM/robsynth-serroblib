% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:53
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRPRR9_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR9_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR9_jacobiR_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:53:24
% EndTime: 2019-02-26 20:53:25
% DurationCPUTime: 0.29s
% Computational Cost: add. (375->45), mult. (1060->92), div. (0->0), fcn. (1462->16), ass. (0->54)
t239 = sin(pkin(7));
t237 = sin(pkin(13));
t241 = cos(pkin(13));
t247 = sin(qJ(3));
t251 = cos(qJ(3));
t257 = t251 * t237 + t247 * t241;
t225 = t257 * t239;
t243 = cos(pkin(7));
t227 = t257 * t243;
t244 = cos(pkin(6));
t242 = cos(pkin(12));
t252 = cos(qJ(1));
t258 = t252 * t242;
t238 = sin(pkin(12));
t248 = sin(qJ(1));
t262 = t248 * t238;
t229 = -t244 * t258 + t262;
t259 = t252 * t238;
t261 = t248 * t242;
t230 = t244 * t259 + t261;
t233 = t247 * t237 - t251 * t241;
t240 = sin(pkin(6));
t264 = t240 * t252;
t213 = t225 * t264 + t229 * t227 + t230 * t233;
t221 = -t229 * t239 + t243 * t264;
t246 = sin(qJ(5));
t250 = cos(qJ(5));
t204 = t213 * t250 + t221 * t246;
t245 = sin(qJ(6));
t249 = cos(qJ(6));
t224 = t233 * t239;
t226 = t233 * t243;
t255 = t224 * t264 + t229 * t226 - t230 * t257;
t269 = t204 * t245 - t249 * t255;
t268 = t204 * t249 + t245 * t255;
t202 = t213 * t246 - t221 * t250;
t265 = t240 * t248;
t263 = t245 * t250;
t260 = t249 * t250;
t231 = -t244 * t261 - t259;
t256 = -t231 * t239 + t243 * t265;
t232 = -t244 * t262 + t258;
t254 = t225 * t265 + t231 * t227 - t232 * t233;
t253 = t244 * t225 + (t227 * t242 - t233 * t238) * t240;
t228 = -t240 * t242 * t239 + t244 * t243;
t218 = -t244 * t224 + (-t226 * t242 - t238 * t257) * t240;
t215 = -t224 * t265 - t231 * t226 - t232 * t257;
t208 = t228 * t246 + t250 * t253;
t207 = t228 * t250 - t246 * t253;
t206 = t246 * t256 + t250 * t254;
t205 = t246 * t254 - t250 * t256;
t201 = t206 * t249 - t215 * t245;
t200 = -t206 * t245 - t215 * t249;
t1 = [t268, 0, t215 * t260 + t245 * t254, 0, -t205 * t249, t200; t201, 0, -t213 * t245 + t255 * t260, 0, t202 * t249, t269; 0, 0, t218 * t260 + t245 * t253, 0, t207 * t249, -t208 * t245 - t218 * t249; -t269, 0, -t215 * t263 + t249 * t254, 0, t205 * t245, -t201; t200, 0, -t213 * t249 - t255 * t263, 0, -t202 * t245, t268; 0, 0, -t218 * t263 + t249 * t253, 0, -t207 * t245, -t208 * t249 + t218 * t245; t202, 0, t215 * t246, 0, t206, 0; t205, 0, t255 * t246, 0, -t204, 0; 0, 0, t218 * t246, 0, t208, 0;];
JR_rot  = t1;
