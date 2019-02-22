% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRPR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 10:50
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRRPR12_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR12_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRPR12_jacobiR_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 10:50:30
% EndTime: 2019-02-22 10:50:31
% DurationCPUTime: 0.25s
% Computational Cost: add. (234->47), mult. (692->91), div. (0->0), fcn. (956->14), ass. (0->54)
t226 = cos(pkin(6));
t221 = sin(pkin(12));
t234 = cos(qJ(1));
t241 = t234 * t221;
t224 = cos(pkin(12));
t230 = sin(qJ(1));
t242 = t230 * t224;
t216 = t226 * t241 + t242;
t229 = sin(qJ(3));
t233 = cos(qJ(3));
t240 = t234 * t224;
t243 = t230 * t221;
t215 = -t226 * t240 + t243;
t225 = cos(pkin(7));
t222 = sin(pkin(7));
t223 = sin(pkin(6));
t247 = t223 * t234;
t239 = t222 * t247;
t237 = t215 * t225 + t239;
t204 = -t216 * t233 + t237 * t229;
t209 = -t215 * t222 + t225 * t247;
t228 = sin(qJ(4));
t232 = cos(qJ(4));
t196 = t204 * t228 - t209 * t232;
t227 = sin(qJ(6));
t258 = t196 * t227;
t231 = cos(qJ(6));
t257 = t196 * t231;
t256 = t204 * t232 + t209 * t228;
t236 = t226 * t242 + t241;
t248 = t223 * t230;
t253 = -t222 * t248 + t236 * t225;
t252 = t216 * t229;
t250 = t222 * t226;
t249 = t223 * t224;
t246 = t225 * t233;
t245 = t227 * t228;
t244 = t228 * t231;
t217 = -t226 * t243 + t240;
t214 = -t222 * t249 + t226 * t225;
t211 = t236 * t222 + t225 * t248;
t208 = t229 * t250 + (t224 * t225 * t229 + t221 * t233) * t223;
t207 = t223 * t221 * t229 - t233 * t250 - t246 * t249;
t206 = t217 * t233 - t253 * t229;
t205 = t217 * t229 + t253 * t233;
t203 = -t237 * t233 - t252;
t201 = t215 * t246 + t233 * t239 + t252;
t200 = t208 * t232 + t214 * t228;
t199 = t208 * t228 - t214 * t232;
t198 = t206 * t232 + t211 * t228;
t197 = t206 * t228 - t211 * t232;
t193 = t197 * t227 + t205 * t231;
t192 = t197 * t231 - t205 * t227;
t1 = [t203 * t231 + t258, 0, -t205 * t245 + t206 * t231, t198 * t227, 0, t192; t193, 0, -t201 * t245 - t204 * t231, -t256 * t227, 0, -t201 * t227 - t257; 0, 0, -t207 * t245 + t208 * t231, t200 * t227, 0, t199 * t231 - t207 * t227; -t203 * t227 + t257, 0, -t205 * t244 - t206 * t227, t198 * t231, 0, -t193; t192, 0, -t201 * t244 + t204 * t227, -t256 * t231, 0, -t201 * t231 + t258; 0, 0, -t207 * t244 - t208 * t227, t200 * t231, 0, -t199 * t227 - t207 * t231; t256, 0, -t205 * t232, -t197, 0, 0; t198, 0, -t201 * t232, t196, 0, 0; 0, 0, -t207 * t232, -t199, 0, 0;];
JR_rot  = t1;
